/*
 * Created on Apr 10, 2017
 *
 */
package org.reactome.cancer.driver;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math.stat.correlation.SpearmansCorrelation;
import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.persistence.MySQLAdaptor;
import org.jgrapht.Graph;
import org.jgrapht.GraphPath;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.alg.shortestpath.AStarShortestPath;
import org.jgrapht.alg.shortestpath.GraphMeasurer;
import org.jgrapht.graph.AsSubgraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.Pseudograph;
import org.junit.Test;
import org.reactome.annotate.AnnotationHelper;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.cancer.*;
import org.reactome.r3.Interactome3dAnalyzer;
import org.reactome.r3.ReactionMapGenerator;
import org.reactome.r3.ReactomeAnalyzer;
import org.reactome.r3.graph.BreadthFirstSearch;
import org.reactome.r3.util.*;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author gwu
 */
public class MechismoAnalyzer {
    private final static Logger logger = Logger.getLogger(MechismoOutputLoader.class);
    private String dirName = "datasets/Mechismo/";
    private String outputFileName = dirName + "COSMICv74_somatic_noSNPs_GWS_mechismo_output.tsv";
    private String pciContactFile = dirName + "human_pci_contact_hits.tsv";
    private FileUtility fu = new FileUtility();

    /**
     * Default constructor.
     */
    public MechismoAnalyzer() {
    }

    public static void main(String[] args) throws Exception {
        MechismoAnalyzer analyzer = new MechismoAnalyzer();
        analyzer.analyzeDistribution();
    }

    @Test
    public void generateSampleToReactionMatrix() throws IOException {
        String resultDirName = dirName + "FrancescoResults/";
        // Based a new analysis result
        String reactioneFileName = resultDirName + "tcga_mechismo_stat_cancer_wise.tsv";
        Map<String, Set<String>> fiToReactions = new HashMap<>();
        Map<String, Set<String>> cancerToReactions = new HashMap<>();
        Files.lines(Paths.get(reactioneFileName))
                .skip(1)
                .map(line -> line.split("\t"))
                .filter(tokens -> Double.parseDouble(tokens[11]) < 0.01d) // FDR cutoff = 0.01
                .filter(tokens -> !tokens[13].equals("-")) // Make sure we have DB_IDs
                .forEach(tokens -> {
                    String fi = InteractionUtilities.generateFIFromGene(tokens[1], tokens[2]);
                    String[] reactionIds = tokens[13].split(",");
                    for (String reactionId : reactionIds) {
                        fiToReactions.compute(fi, (key, set) -> {
                            if (set == null)
                                set = new HashSet<>();
                            set.add(reactionId);
                            return set;
                        });
                        cancerToReactions.compute(tokens[0], (key, set) -> {
                            if (set == null)
                                set = new HashSet<>();
                            set.add(reactionId);
                            return set;
                        });
                    }
                });

        System.out.println("Total reactions: " + fiToReactions.size());
        // Load TCGA samples
        String mechResultFileName = dirName + "TCGA/TCGA_mech_output.tsv";
        Map<String, Set<String>> sampleToFIs = new HashMap<>(); // Sample to mutated FIs
        Map<String, String> sampleToCancerType = new HashMap<>();
        Files.lines(Paths.get(mechResultFileName))
                .skip(1) // Skip the first line
                .map(line -> line.split("\t"))
                .filter(tokens -> !tokens[0].equals("name_a1")) // Escape the head lines inside the file
                .filter(tokens -> !tokens[4].equals(tokens[5])) // Escape synonymous mutations
                .filter(tokens -> tokens.length > 18) // Need partner's information
                .filter(tokens -> !tokens[18].equals("[PROT]")) // Escape self-interaction
                .forEach(tokens -> {
                    String fi = InteractionUtilities.generateFIFromGene(tokens[0], tokens[18]);
                    String[] tokensTmp = tokens[6].split(" ");
                    String[] sampleTokens = tokensTmp[tokensTmp.length - 1].split(";");
                    for (String sampleToken : sampleTokens) {
                        tokensTmp = sampleToken.split(":");
                        String cancer = tokensTmp[0];
                        tokensTmp = tokensTmp[1].split(",");
                        for (String sample : tokensTmp) {
                            sampleToCancerType.put(sample, cancer);
                            sampleToFIs.compute(sample, (key, set) -> {
                                if (set == null)
                                    set = new HashSet<>();
                                set.add(fi);
                                return set;
                            });
                        }
                    }
                });
        System.out.println("Total samples: " + sampleToCancerType.size());
        // Just a quick check
        Map<String, Set<String>> cancerTypeToSamples = new HashMap<>();
        sampleToCancerType.forEach((sample, cancer) -> cancerTypeToSamples.compute(cancer, (key, set) -> {
            if (set == null)
                set = new HashSet<>();
            set.add(sample);
            return set;
        }));
        cancerTypeToSamples.forEach((cancer, samples) -> System.out.println(cancer + "\t" + samples.size()));
        // Generate a map from samples to hit reactions
        Map<String, Set<String>> sampleToReactions = new HashMap<>();
        sampleToFIs.forEach((sample, fis) -> {
            String cancer = sampleToCancerType.get(sample);
            Set<String> cancerReactions = cancerToReactions.get(cancer);
            if (cancerReactions == null)
                return; // Don't do anything
            Set<String> sampleFIs = sampleToFIs.get(sample);
            Set<String> fiReactions = new HashSet<>();
            for (String fi : sampleFIs) {
                Set<String> tmp = fiToReactions.get(fi);
                if (tmp != null)
                    fiReactions.addAll(tmp);
            }
            fiReactions.retainAll(cancerReactions);
            sampleToReactions.put(sample, fiReactions);
        });
        sampleToReactions.forEach((sample, reactions) -> System.out.println(sample + "\t" + reactions.size()));

        // Output into a local file
        String outFileName = "results/MechismoSamplesToReactions_103017.txt";
        fu.setOutput(outFileName);
        Set<String> reactions = sampleToReactions.values().stream().flatMap(set -> set.stream()).collect(Collectors.toSet());
        List<String> reactionList = new ArrayList<>(reactions);
        Collections.sort(reactionList);
        StringBuilder builder = new StringBuilder();
        builder.append("Cancer\tSample");
        reactionList.forEach(reaction -> builder.append("\t").append(reaction));
        fu.printLine(builder.toString());
        builder.setLength(0);
        for (String sample : sampleToReactions.keySet()) {
            Set<String> sampleReactions = sampleToReactions.get(sample);
            String cancer = sampleToCancerType.get(sample);
            builder.append(cancer).append("\t").append(sample);
            for (String reaction : reactionList) {
                int out = sampleReactions.contains(reaction) ? 1 : 0;
                builder.append("\t").append(out);
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }

    @Test
    public void generateReactionToTCGACancerMatrix() throws IOException {
        String resultDirName = dirName + "FrancescoResults/";

        Map<String, Set<String>> reactionToCancers = new HashMap<>();
        Set<String> cancers = new HashSet<>();

//        String reactionFileName = resultDirName + "TCGA_mechismo_reactome_091417_reaction.txt";
//        String output = "results/MechismoReactionToTCGACancer_092817.txt";
//        Files.lines(Paths.get(reactionFileName))
//             .skip(1)
//             .map(line -> line.split("\t"))
//             .forEach(tokens -> {
//                 cancers.add(tokens[0]);
//                 reactionToCancers.compute(tokens[1], (reaction, set) -> {
//                    if (set == null)
//                        set = new HashSet<>();
//                    set.add(tokens[0]);
//                    return set;
//                 });
//             });

        // Based a new analysis result
        String reactioneFileName = resultDirName + "tcga_mechismo_stat_cancer_wise.tsv";
        String output = "results/MechismoReactionToTCGACancer_100217.txt";
        Files.lines(Paths.get(reactioneFileName))
                .skip(1)
                .map(line -> line.split("\t"))
                .filter(tokens -> Double.parseDouble(tokens[11]) < 0.01d) // FDR cutoff = 0.01
                .filter(tokens -> !tokens[13].equals("-")) // Make sure we have DB_IDs
                .forEach(tokens -> {
                    cancers.add(tokens[0]);
                    String[] reactionIds = tokens[13].split(",");
                    for (String reactionId : reactionIds) {
                        reactionToCancers.compute(reactionId, (reaction, set) -> {
                            if (set == null)
                                set = new HashSet<>();
                            set.add(tokens[0]);
                            return set;
                        });
                    }
                });

        System.out.println("Total cancers: " + cancers.size());
        System.out.println("Total reactions: " + reactionToCancers.size());

        // Output the file
        List<String> cancerList = new ArrayList<>(cancers);
        Collections.sort(cancerList);

        if (true)
            return;

        fu.setOutput(output);
        StringBuilder builder = new StringBuilder();
        builder.append("Reaction");
        cancerList.forEach(cancer -> builder.append("\t").append(cancer));
        fu.printLine(builder.toString());
        builder.setLength(0);
        reactionToCancers.forEach((reaction, set) -> {
            builder.append(reaction);
            cancerList.forEach(cancer -> {
                builder.append("\t");
                builder.append(set.contains(cancer));
            });
            try {
                fu.printLine(builder.toString());
            } catch (IOException e) {
                e.printStackTrace();
            }
            builder.setLength(0);
        });
        fu.close();
    }

    @Test
    public void checkSampleToHitReactions() throws IOException {
        String srcFileName = "results/MechismoSamplesToReactions_103017.txt";
        Map<String, Set<String>> cancerToSamples = new HashMap<>();
        Map<String, Set<String>> cancerToHitSamples = new HashMap<>();
        fu.setInput(srcFileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            cancerToSamples.compute(tokens[0], (key, set) -> {
                if (set == null)
                    set = new HashSet<>();
                set.add(tokens[1]);
                return set;
            });
            for (int i = 2; i < tokens.length; i++) {
                if (tokens[i].equals("1")) {
                    cancerToHitSamples.compute(tokens[0], (key, set) -> {
                        if (set == null)
                            set = new HashSet<>();
                        set.add(tokens[1]);
                        return set;
                    });
                    break;
                }
            }
        }
        fu.close();

        System.out.println("Cancer\tTotalSamples\tHitSamples\tPercentage");
        for (String cancer : cancerToSamples.keySet()) {
            Set<String> samples = cancerToSamples.get(cancer);
            Set<String> hitSamples = cancerToHitSamples.get(cancer);
            Double percent = (double) hitSamples.size() / samples.size();
            System.out.println(cancer + "\t" +
                    samples.size() + "\t" +
                    hitSamples.size() + "\t" +
                    percent);
        }
    }

    @Test
    public void analyzeSampleClusterReactions() throws Exception {
        String cancer = "LGG";
        cancer = "UCEC";
        cancer = "HNSC";

        String fileName = "results/" + cancer + "TwoClusters_112917.txt";

        String srcFileName = "datasets/Mechismo/FrancescoResults/110317/tcga_mechismo_stat_cancer_wise_significant.tsv";
        double fdrCutoff = 0.01d;
        boolean hasHeader = false;

        Map<String, Set<String>> sampleToReactions = loadSampleToReactions(srcFileName,
                cancer,
                fdrCutoff,
                hasHeader);

        // Generate reaction names
        CancerDriverReactomeAnalyzer reactomeAnalyzer = new CancerDriverReactomeAnalyzer();
        MySQLAdaptor dba = reactomeAnalyzer.getDBA();


        Map<String, Set<String>> clusterToSamples = loadSampleCluseters(fileName);
        Map<String, List<String>> clusterToReactions = new HashMap<>();
        clusterToSamples.forEach((cluster, samples) -> {
            Set<String> reactions = new HashSet<>();
            samples.forEach(sample -> reactions.addAll(sampleToReactions.get(sample)));
            List<String> idList = new ArrayList<>(reactions);
            idList.sort(Comparator.naturalOrder());
            System.out.println(cluster + ": " + samples.size() + " samples, " + reactions.size() + " reactions");
            try {
                for (String id : idList) {
                    GKInstance reaction = dba.fetchInstance(new Long(id));
                    System.out.println(id + "\t" + reaction.getDisplayName());
                }
            } catch (Exception e) {
            }
            System.out.println();
            clusterToReactions.put(cluster, idList);
        });

        // Get shared reactions
        List<String> shared = new ArrayList<>();
        clusterToReactions.forEach((cluster, reactions) -> {
            if (shared.size() == 0)
                shared.addAll(reactions);
            else
                shared.retainAll(reactions);
        });
        System.out.println("\nShared reactions: " + shared.size());
        for (String id : shared) {
            GKInstance reaction = dba.fetchInstance(new Long(id));
            System.out.println(id + "\t" + reaction.getDisplayName());
        }

        // Get different reactions
        System.out.println("\nUnique reactions:");
        for (String cluster : clusterToReactions.keySet()) {
            List<String> reactions = clusterToReactions.get(cluster);
            reactions.removeAll(shared);
            System.out.println("\n" + cluster + ": " + reactions.size() + " reactions");
            for (String id : reactions) {
                GKInstance reaction = dba.fetchInstance(new Long(id));
                System.out.println(id + "\t" + reaction.getDisplayName());
            }
        }
    }

    /**
     * Load sample cluster results output from R script.
     *
     * @param fileName
     * @return
     * @throws IOException
     */
    private Map<String, Set<String>> loadSampleCluseters(String fileName) throws IOException {
        Map<String, Set<String>> clusterToSamples = new HashMap<>();
        fu.setInput(fileName);
        String line = null;
        String cluster = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("Total clusters"))
                continue;
            if (line.startsWith("Cluster ")) {
                int index = line.indexOf(":");
                cluster = line.substring(0, index).trim();
                continue;
            }
            if (cluster != null) {
                String[] tokens = line.split("\t");
                Set<String> set = Arrays.asList(tokens).stream().collect(Collectors.toSet());
                clusterToSamples.put(cluster, set);
                cluster = null;
            }
        }

        fu.close();
        return clusterToSamples;
    }

    @Test
    public void calculateSampleDistancesOnReactionNetwork() throws IOException {

        String cancer = "GBM";
        cancer = "PAAD";
        cancer = "SKCM";
//        cancer = "LGG";
//        cancer = "UCEC";
        cancer = "COADREAD";
//        cancer = "THCA";
//        cancer = "STAD";
//        cancer = "LUAD";
//        cancer = "HNSC";

        // Generate distances based on significant reactions only
//        String srcFileName = "datasets/Mechismo/FrancescoResults/110317/tcga_mechismo_stat_cancer_wise_significant.tsv";
//        double fdrCutoff = 0.01d;
//        boolean hasHeader = false;

        // Generate distances based on all reactions
        String srcFileName = "datasets/Mechismo/FrancescoResults/110317/tcga_mechismo_stat_cancer_wise.tsv";
        double fdrCutoff = 1.01d;
        boolean hasHeader = true;

        Map<String, Set<String>> sampleToReactions = loadSampleToReactions(srcFileName,
                cancer,
                fdrCutoff,
                hasHeader);

        // A quick quality check
        sampleToReactions.forEach((sample, set) -> System.out.println(sample + "\t" + set.size()));

        filterSamplesWithoutReactions(sampleToReactions);

        // Some test code for two UCEC samples have many hit reactions
//        String sample1 = "TCGA-A5-A0G3-01";
//        String sample2 = "TCGA-A5-A0G9-01";
//        calculateNetworkDistance(sampleToReactions, sample1, sample2);
//        
//        if (true)
//            return;

        String date = "103117";
        date = "111317";
        date = "111617"; // Should be a very quick way to calculate distances in BFS
        date = "112817"; // Use all reactions

        String fileName = "results/MechismoSamplePairWiseReactionNetworkDist_" + cancer + "_" + date + ".txt";

        generatePairWiseNetworkDistance(sampleToReactions, fileName);
    }

    private Map<String, Set<String>> loadSampleToReactions(String cancer) throws IOException {
        String srcFileName = "results/MechismoSamplesToReactions_103017.txt";
        Map<String, Set<String>> sampleToReactions = new HashMap<>();
        fu.setInput(srcFileName);
        String line = fu.readLine();
        String[] reactions = line.split("\t");
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (!tokens[0].equals(cancer))
                continue;
            Set<String> set = new HashSet<>();
            for (int i = 2; i < tokens.length; i++) {
                if (tokens[i].equals("1"))
                    set.add(reactions[i]);
            }
            sampleToReactions.put(tokens[1], set);
        }
        fu.close();
        return sampleToReactions;
    }

    /**
     * Use this method to load sample to hit reactions using updated files from Francesco, which contain samples
     * in the analysis output.
     *
     * @param cancer
     * @param fdrCutoff
     * @return
     * @throws IOException
     */
    private Map<String, Set<String>> loadSampleToReactions(String srcFileName,
                                                           String cancer,
                                                           double fdrCutoff,
                                                           boolean hasHeader) throws IOException {
        Map<String, Set<String>> sampleToReactions = new HashMap<>();
        try (Stream<String> stream = Files.lines(Paths.get(srcFileName))) {
            stream.skip(hasHeader ? 1 : 0)
                    .map(line -> line.split("\t"))
                    .filter(tokens -> tokens[0].equals(cancer)) // Filter to cancer type
                    .filter(tokens -> !tokens[13].equals("-")) // Need to have reactions
                    .filter(tokens -> Double.parseDouble(tokens[11]) < fdrCutoff) // Less than the predefined fdr cutoff
                    .forEach(tokens -> {
                        String[] reactionIds = tokens[13].split(",");
                        // Get samples
                        int index1 = tokens[15].indexOf("\'");
                        int index2 = tokens[15].lastIndexOf("\'");
                        String tmp = tokens[15].substring(index1, index2 + 1);
                        String[] samples = tmp.split(", ");
                        Arrays.asList(samples).forEach(sample -> {
                            // Need to remove single quotes
                            sample = sample.substring(1, sample.length() - 1);
                            sampleToReactions.compute(sample, (key, set) -> {
                                if (set == null)
                                    set = new HashSet<>();
                                set.addAll(Arrays.asList(reactionIds));
                                return set;
                            });
                        });
                    });
        }
        return sampleToReactions;
    }

    private void filterSamplesWithoutReactions(Map<String, Set<String>> sampleToReactions) {
        // Filter out samples don't have any hit pathways
        for (Iterator<String> it = sampleToReactions.keySet().iterator(); it.hasNext(); ) {
            String sample = it.next();
            Set<String> set = sampleToReactions.get(sample);
            if (set.size() == 0)
                it.remove();
        }
    }

    /**
     * Calculate a pair-wise distance matrix for TCGA cancer types based
     * on a list of significant reactions generated by Francesco.
     *
     * @throws IOException
     */
    @Test
    public void calculateCancerTypeDistancesOnReactionNetwork() throws IOException {
        // Load cancer type to significant reactions
        Map<String, Set<String>> cancerToReactions = new HashMap<>();
        String resultDirName = dirName + "FrancescoResults/";

//        String reactionFileName = resultDirName + "TCGA_mechismo_reactome_091417_reaction.txt";
//        Files.lines(Paths.get(reactionFileName))
//             .skip(1)
//             .map(line -> line.split("\t"))
//             .filter(tokens -> tokens.length > 1)
//             .forEach(tokens -> {
//                 cancerToReactions.compute(tokens[0], (key, set) -> {
//                     if (set == null)
//                         set = new HashSet<>();
//                     set.add(tokens[1]);
//                     return set;
//                 });
//             });

        // Based a new analysis result
        String reactioneFileName = resultDirName + "tcga_mechismo_stat_cancer_wise.tsv";
        Files.lines(Paths.get(reactioneFileName))
                .skip(1)
                .map(line -> line.split("\t"))
                .filter(tokens -> Double.parseDouble(tokens[11]) < 0.01d) // FDR cutoff = 0.01
                .filter(tokens -> !tokens[13].equals("-")) // Make sure we have DB_IDs
                .forEach(tokens -> {
                    String[] reactionIds = tokens[13].split(",");
                    for (String reactionId : reactionIds) {
                        cancerToReactions.compute(tokens[0], (key, set) -> {
                            if (set == null)
                                set = new HashSet<>();
                            set.add(reactionId);
                            return set;
                        });
                    }
                });

        String fileName = "results/TCGACancerPairWiseReactionNetworkDist_100217.txt";
        generatePairWiseNetworkDistance(cancerToReactions, fileName);
    }

    void calculateNetworkDistance(Map<String, Set<String>> sampleToReactions,
                                  String sample1,
                                  String sample2) throws IOException {
        // Calculate pair-wise average distances between two cancer types
        ReactionMapGenerator mapGenerator = new ReactionMapGenerator();
        Set<String> network = mapGenerator.loadLargestComponents();
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> idToPartners = bfs.generateIdToPartnersMap(network);
        // Need a little bit cleanup
        sampleToReactions.forEach((cancer, reactions) -> reactions.retainAll(idToPartners.keySet()));

        Set<String> set1 = sampleToReactions.get(sample1);
        Set<String> set2 = sampleToReactions.get(sample2);
        long time1 = System.currentTimeMillis();
        double dist = bfs.calculateMinShortestPath(set1, set2, idToPartners);
        long time2 = System.currentTimeMillis();
        System.out.println(dist);
        System.out.println("Total time: " + (time2 - time1));
    }

    /**
     * Use this method to generate the pair-wise distance for samples listed in the passed parameter.
     *
     * @param sampleToReactions
     * @throws IOException
     */
    private void generatePairWiseNetworkDistance(Map<String, Set<String>> sampleToReactions,
                                                 String fileName) throws IOException {
        // Calculate pair-wise average distances between two cancer types
        ReactionMapGenerator mapGenerator = new ReactionMapGenerator();
        Set<String> network = mapGenerator.loadLargestComponents();
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> idToPartners = bfs.generateIdToPartnersMap(network);

        // Need a little bit cleanup
        sampleToReactions.forEach((cancer, reactions) -> reactions.retainAll(idToPartners.keySet()));
        filterSamplesWithoutReactions(sampleToReactions);
        List<String> samples = sampleToReactions.keySet().stream().sorted(Comparator.naturalOrder()).collect(Collectors.toList());

        // Generate output
        StringBuilder builder = new StringBuilder();
        builder.append("Cancer");
        samples.forEach(cancer -> builder.append("\t").append(cancer));
        fu.setOutput(fileName);
        fu.printLine(builder.toString());
        builder.setLength(0);

        // Use a helper object
        MechismoBFSHelper helper = new MechismoBFSHelper();
        Map<String, Integer> pairToDist = helper.calculateShortestPath(sampleToReactions,
                idToPartners,
                bfs);
        for (int i = 0; i < samples.size(); i++) {
            String sample1 = samples.get(i);
            builder.append(sample1);
            Set<String> set1 = sampleToReactions.get(sample1);
            for (int j = 0; j < samples.size(); j++) {
                String sample2 = samples.get(j);
                builder.append("\t");
                if (i == j)
                    builder.append(0.0d);
                else {
                    Set<String> set2 = sampleToReactions.get(sample2);
                    double dist = helper.calculateMinShortestPath(set1,
                            set2,
                            pairToDist);
                    builder.append(dist);
                }
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }

    public Map<String, Double> loadReactionToMaxMechismoScore() throws IOException {
        String mechismoFile = "results/ReactionsInMechismo_051017.txt";
        Map<String, Double> reactionToMechismo = Files.lines(Paths.get(mechismoFile))
                .skip(1)
                .map(line -> line.split("\t"))
                .filter(tokens -> tokens[tokens.length - 1].length() > 0)
                .collect(Collectors.toMap(tokens -> tokens[1],
                        tokens -> new Double(tokens[tokens.length - 1]),
                        (v1, v2) -> (Math.abs(v1) > Math.abs(v2) ? v1 : v2)));
        System.out.println("Total reactions in mechismo scores: " + reactionToMechismo.size());
        return reactionToMechismo;
    }

    /**
     * Analyze correlations between different measurements for Reactome reactions.
     *
     * @throws IOException
     */
    @Test
    public void analyzeCorrelations() throws Exception {
        String resultDir = "results/";

        // Check correlation between cancer driver gene enrichment scores and maximum Mechismo scores
        Map<String, Double> reactionToDriverEnrichment = new CancerDriverReactomeAnalyzer().loadReactionToCancerGeneEnrichment();
        Map<String, Double> reactionToMechismo = loadReactionToMaxMechismoScore();
        Map<String, Double> reactionToInteractome3dScore = new Interactome3dDriverAnalyzer().loadReactionTo3dScore();

        System.out.println("\nCorrelation between driver enrichment and mechismo score:");
        analyzeCorrelation(reactionToDriverEnrichment, reactionToMechismo);

        System.out.println("\nCorrelation between driver enrichment and interaction3d score:");
        analyzeCorrelation(reactionToDriverEnrichment, reactionToInteractome3dScore);

        System.out.println("\nCorrelation between interaction3d score and driver enrichment:");
        analyzeCorrelation(reactionToInteractome3dScore, reactionToDriverEnrichment);

        System.out.println("\nCorrelation between interaction3d and mechismo score:");
        analyzeCorrelation(reactionToInteractome3dScore, reactionToMechismo);
    }

    protected void analyzeCorrelation(Map<String, Double> reactionToScore1,
                                      Map<String, Double> reactionToScore2) throws MathException {
        List<Double> values1 = new ArrayList<>();
        List<Double> values2 = new ArrayList<>();
//        double scoreCutoff = 1.3d;
        double scoreCutoff = 2.0d;
        reactionToScore1.forEach((reaction, enrichment) -> {
            if (enrichment < scoreCutoff) // Choose cutoff < 0.01
                return;
            if (reactionToScore2.containsKey(reaction)) {
                values1.add(enrichment);
                values2.add(Math.abs(reactionToScore2.get(reaction)));
            }
        });
        System.out.println("Score cutoff: " + scoreCutoff);
        System.out.println("Total selected reactions for correlation: " + values1.size());

        PearsonsCorrelation correlation = MathUtilities.constructPearsonCorrelation(values1, values2);
        System.out.println(correlation.getCorrelationMatrix().getEntry(0, 1) + ", " +
                correlation.getCorrelationPValues().getEntry(0, 1));

        SpearmansCorrelation spearman = MathUtilities.constructSpearmansCorrelation(values1, values2);
        correlation = spearman.getRankCorrelation();
        System.out.printf("Correlation (rank): %f, correlation (via Pearson): %f, p-value: %f\n",
                spearman.getCorrelationMatrix().getEntry(0, 1),
                correlation.getCorrelationMatrix().getEntry(0, 1),
                correlation.getCorrelationPValues().getEntry(0, 1));

        performFisherTest(reactionToScore1,
                reactionToScore2);
    }

    protected void performFisherTest(Map<String, Double> reactionToScore1,
                                     Map<String, Double> reactionToScore2) {
        // Perform a Fisher exact test for sharing top reactions
        List<String> list1 = new ArrayList<>(reactionToScore1.keySet());
        List<String> list2 = new ArrayList<>(reactionToScore2.keySet());
        // Keep only shared reactions
        list1.retainAll(list2);
        list2.retainAll(list1);
        System.out.println("Total shared reactions: " + list1.size() + ", " + list2.size());
        // Sort these two lists
        list1.sort((r1, r2) -> {
            Double s1 = reactionToScore1.get(r1);
            Double s2 = reactionToScore1.get(r2);
            return s2.compareTo(s1);
        });
        list2.sort((r1, r2) -> {
            Double s1 = reactionToScore2.get(r1);
            Double s2 = reactionToScore2.get(r2);
            return s2.compareTo(s1);
        });
        // Take top 100
        int top = 100;
        Set<String> topShared = InteractionUtilities.getShared(list1.subList(0, top),
                list2.subList(0, top));
        System.out.println("Chosen top reactions: " + top);
        System.out.println("Top shared: " + topShared.size());
        FisherExact fisherExact = new FisherExact(list1.size());
        double pvalue = fisherExact.getRightTailedP(topShared.size(),
                top - topShared.size(),
                top - topShared.size(),
                list1.size() - 2 * top + topShared.size());
        System.out.println("p-value from Fisher Exact test: " + pvalue);
    }

    @Test
    public void analyzeDistribution() throws IOException {
        // Get scores involved known cancer genes only
        Set<String> cancerGenes = new CancerDriverAnalyzer().getDriverGenes(null);
        System.out.println("Total cancer genes: " + cancerGenes.size());

        Path filePath = Paths.get(outputFileName);

        List<Double> values = new ArrayList<>();
        List<Double> cancerGeneValues = new ArrayList<>();
        Map<String, Double> geneToMaxValue = new HashMap<>();

        // Get scores
        Files.lines(filePath)
                .map(line -> line.split("\t"))
                .filter(tokens -> tokens.length >= 28)
//        .filter(tokens -> tokens[19].trim().length() > 0)
                .forEach(tokens -> {
                    if (tokens[27].trim().length() == 0)
                        return;
                    Double value = new Double(tokens[27]);
                    values.add(value);
                    if (cancerGenes.contains(tokens[0]))
                        cancerGeneValues.add(value);
                    BiFunction<String, Double, Double> func = (k, v) ->
                            (v == null || Math.abs(v) < Math.abs(value)) ? value : v;
                    geneToMaxValue.compute(tokens[0], func);
                });

        System.out.println("Total values: " + values.size());
        System.out.println("\tMinimum value: " + values.stream().mapToDouble(Double::doubleValue).min());
        System.out.println("Total cancer gene values: " + cancerGeneValues.size());
        System.out.println("\tMinimum value: " + cancerGeneValues.stream().mapToDouble(Double::doubleValue).min());
        System.out.println("Max values: " + geneToMaxValue.size());
        System.out.println("\tMinimum value: " + geneToMaxValue.values().stream().mapToDouble(Double::doubleValue).min());

        // Two plots
        // All values
        Plotter plotter = new Plotter("Plot");
        Map<String, List<Double>> datasetToName = new HashMap<>();
        datasetToName.put("All Scores", values);
        datasetToName.put("Cancer Gene Scores", cancerGeneValues);
        plotter.plotHistograpm(datasetToName,
                11,
                "Histogram of Mechismo Score",
                "Mechismo Score",
                "Frequency");
        // Maximum values
        plotter = new Plotter("Plot");
        datasetToName = new HashMap<>();
        datasetToName.put("All Genes", new ArrayList<Double>(geneToMaxValue.values()));
        Map<String, Double> cancerGeneToMax = new HashMap<>();
        cancerGeneToMax.putAll(geneToMaxValue);
        cancerGeneToMax.keySet().retainAll(cancerGenes);
        datasetToName.put("Cancer Genes", new ArrayList<Double>(cancerGeneToMax.values()));
        plotter.plotHistograpm(datasetToName,
                11,
                "Histogram of Maximum Mechismo Score",
                "Mechismo Maximum Score",
                "Frequency");
    }

    @Test
    public void performPathwayEnrichmentAnalysisForReactions() throws Exception {
        String reactionsFileName = "results/ReactionsInMechismo_051017.txt";

        AnnotationHelper annotationHelper = new AnnotationHelper();
        annotationHelper.setReactionIdToPathwayFile("resources/ReactomeReactionsToPathways_051017.txt");
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        annotator.setUseBenjaminiHochbergForFDR(true);
        annotator.setAnnotationHelper(annotationHelper);

        // Choose a threshold to get reactions to be analyzed
        double[] thresholds = new double[]{5.0d, 4.0d, 3.0d};
        for (double threshold : thresholds) {
            List<String> reactionIds = getReactionIds(reactionsFileName, threshold);
            List<GeneSetAnnotation> annotations = annotator.annotateReactionsWithReactomePathways(reactionIds);
            System.out.println("Threshold: " + threshold + " with total selected reactions: " + reactionIds.size());
            System.out.println("Pathway\tNumberInPathway\tRatioOfPathway\tHitNumber\tpValue\tFDR\tHitIds");
            for (GeneSetAnnotation annotation : annotations) {
                System.out.println(annotation.getTopic() + "\t" +
                        annotation.getNumberInTopic() + "\t" +
                        annotation.getRatioOfTopic() + "\t" +
                        annotation.getHitNumber() + "\t" +
                        annotation.getPValue() + "\t" +
                        annotation.getFdr() + "\t" +
                        annotation.getHitIds());
            }
            System.out.println();
        }
    }

    private List<String> getReactionIds(String fileName,
                                        double threshold) throws IOException {
        fu.setInput(fileName);
        String line = fu.readLine();
        List<String> rtn = new ArrayList<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 5)
                continue;
            double score = new Double(tokens[4]);
            if (Math.abs(score) >= threshold)
                rtn.add(tokens[0]);
        }
        fu.close();
        return rtn;
    }

    /**
     * Load ppi to maximum interaction effect score using the default output filename.
     *
     * @return
     * @throws IOException
     */
    public Map<String, Double> loadPPIToMaxScore() throws IOException {
        return loadPPIToMaxScore(outputFileName);
    }

    @Test
    public void checkAllHumanReactions() throws Exception {
        String output = "results/ReactionsInMechismo_051017.txt";
        output = "results/ReactionsInMechismo_052217.txt";
        CancerDriverReactomeAnalyzer reactomeAnalyzer = new CancerDriverReactomeAnalyzer();

        Long targetReactionId = 5658435L;
        targetReactionId = 5672950L;

        checkAllHumanReactions(reactomeAnalyzer,
                targetReactionId,
                output);
    }

    /**
     * Check all human reactions based on mechismo data.
     *
     * @throws Exception
     */
    public void checkAllHumanReactions(CancerDriverReactomeAnalyzer reactomeAnalyzer,
                                       Long targetReactionId,
                                       String outputFileName) throws Exception {
        // Load all the non-disease reactions from Reactome
        List<GKInstance> reactions = reactomeAnalyzer.loadHumanReactions();
        System.out.println("Total human reactions in Reactome: " + reactions.size());

        Map<String, Double> ppiToScore = loadPPIToMaxScore();
        System.out.println("Total PPIs in mechismo: " + ppiToScore.size());

        ReactomeAnalyzer reactomeDataAnalyzer = new ReactomeAnalyzer();
        if (targetReactionId == null) {
            fu.setOutput(outputFileName);
            fu.printLine("DB_ID\tName\tTotal_FIs\tMechismo_FIs\tMax_Score");
        }
        int total = 0;
        int checkedReactions = 0;
        for (GKInstance reaction : reactions) {
            if (targetReactionId != null && !reaction.getDBID().equals(targetReactionId))
                continue;
            System.out.println("Checking " + reaction + "...");
            Set<String> fis = reactomeDataAnalyzer.generateTentativePPIsForReaction(reaction,
                    false);
            System.out.println("\tExtracted FIs: " + fis.size());
            if (fis.size() == 0)
                continue;
            // Check if there is a structure available from interactome3d
            int overlap = 0;
            Double score = null;
            for (String fi : fis) {
                if (ppiToScore.containsKey(fi)) {
                    overlap++;
                    Double currentScore = ppiToScore.get(fi);
                    if (score == null || Math.abs(currentScore) > Math.abs(score))
                        score = currentScore;
                    System.out.println(fi + "\t" + currentScore);
                }
            }
            if (targetReactionId == null) {
                fu.printLine(reaction.getDBID() + "\t" +
                        reaction.getDisplayName() + "\t" +
                        fis.size() + "\t" +
                        overlap + "\t" +
                        (score == null ? "" : score));
            }
            if (overlap > 0)
                total++;
            checkedReactions++;
        }
        if (targetReactionId == null)
            fu.close();
        System.out.println("Total checked reactions: " + checkedReactions);
        System.out.println("Total reactions having FIs in mechismo: " + total);
        // Output ran on May 10, 2017
//        Total human reactions in Reactome: 9453
//        Total PPIs in mechismo: 13323
//        Total checked reactions: 4589
//        Total reactions having FIs in mechsimo: 1265
    }

    @Test
    public void testLoadPPIToMaxScore() throws IOException {
        Map<String, Double> ppiToScore = loadPPIToMaxScore();
        System.out.println("Total PPIs: " + ppiToScore.size());
        // The above output should be: Total PPIs: 13323
    }

    @Test
    public void checkKnownCancerDriverGenes() throws IOException {
        CancerDriverAnalyzer driverAnalyzer = new CancerDriverAnalyzer();
        Set<String> knownCancerGenes = driverAnalyzer.getDriverGenes(null);
        System.out.println("Total known driver genes: " + knownCancerGenes.size());
    }

    @Test
    public void checkSignificantReactions() throws IOException {
        // Get the list of significant reactions
        String fileName = "results/MechismoReactionList_092517.txt";
        Set<String> significantIds = Files.lines(Paths.get(fileName))
                .skip(1)
                .map(line -> line.split("\t"))
                .map(tokens -> tokens[0])
                .collect(Collectors.toSet());
        System.out.println("Total significant ids: " + significantIds.size());

        // Get ids in the network
        fileName = "results/MechismoAllReactions_092517.txt";
        Set<String> networkIds = Files.lines(Paths.get(fileName))
                .map(line -> line.split(" "))
                .flatMap(tokens -> Arrays.asList(tokens[0], tokens[2]).stream())
                .collect(Collectors.toSet());
        System.out.println("Total reaction ids in network: " + networkIds.size());
        significantIds.removeAll(networkIds);
        System.out.println("Significant ids not in network: " + significantIds.size());
        significantIds.forEach(System.out::println);
    }

    @Test
    public void checkReactionsInAnalysisResults() throws IOException {
        String resultDirName = dirName + "FrancescoResults/";
        // Based on an old analysis results
//        String reactionFileName = resultDirName + "TCGA_mechismo_reactome_091417_reaction.txt";
//        Set<String> reactionIds = Files.lines(Paths.get(reactionFileName))
//                .skip(1)
//                .map(line -> line.split("\t"))
//                .filter(tokens -> tokens.length > 1)
//                .map(tokens -> tokens[1])
//                .collect(Collectors.toSet());
        // Based a new analys result
        String reactioneFileName = resultDirName + "tcga_mechismo_stat_cancer_wise.tsv";
        Set<String> reactionIds = Files.lines(Paths.get(reactioneFileName))
                .skip(1)
                .map(line -> line.split("\t"))
                .filter(tokens -> Double.parseDouble(tokens[11]) < 0.01d) // FDR cutoff = 0.01
                .filter(tokens -> !tokens[13].equals("-")) // Make sure we have DB_IDs
                .map(tokens -> tokens[13].split(","))
                .flatMap(tokens -> Arrays.asList(tokens).stream())
                .collect(Collectors.toSet());

        System.out.println("Total reaction ids: " + reactionIds.size());
        List<String> reactionIdList = new ArrayList<>(reactionIds);
        Collections.sort(reactionIdList);
        System.out.println("ReactionId\tIsSignificant");
        reactionIdList.forEach(id -> System.out.println(id + "\ttrue"));
//        if (true)
//            return;
        // Generate the reaction sub-network
        System.out.println("\nGenerate the reaction network:");
        ReactionMapGenerator mapGenerator = new ReactionMapGenerator();
//        mapGenerator.generateSubNetwork(reactionIds);
        mapGenerator.generateSubNetworkForAll(reactionIds, new HashMap<>());
    }

    /**
     * This is a quick test for results generated by Francesco.
     *
     * @throws IOException
     */
    @Test
    public void checkAnalysisResults() throws IOException {
        String resultDirName = dirName + "FrancescoResults/";
        String reactionFileName = resultDirName + "TCGA_mechismo_reactome_091417_reaction.txt";
        Set<String> reactionIds = Files.lines(Paths.get(reactionFileName))
                .skip(1)
                .map(line -> line.split("\t"))
                .filter(tokens -> tokens.length > 1)
                .map(tokens -> tokens[1])
                .collect(Collectors.toSet());
        System.out.println("Total reaction ids: " + reactionIds.size());
        List<String> reactionIdList = new ArrayList<>(reactionIds);
        Collections.sort(reactionIdList);
        System.out.println("ReactionId\tIsSignificant");
        reactionIdList.forEach(id -> System.out.println(id + "\ttrue"));
        if (true)
            return;
        // Check pathways
        String pathwayFileName = resultDirName + "TCGA_mechismo_reactome_091417_pathway.txt";
        Set<String> pathwayIds = Files.lines(Paths.get(pathwayFileName))
                .skip(1)
                .map(line -> line.split("\t"))
                .filter(tokens -> tokens.length > 1)
                .map(tokens -> tokens[1])
                .collect(Collectors.toSet());
        System.out.println("Total pathway ids: " + pathwayIds.size());

        String interfaceFileName = resultDirName + "TCGA_mechismo_reactome_091417_interface.txt";
        Set<String> fis = new HashSet<>();
        Set<String> cancerTypes = new HashSet<>();
        Files.lines(Paths.get(interfaceFileName))
                .skip(1)
                .forEach(line -> {
                    String[] tokens = line.split("\t");
                    cancerTypes.add(tokens[0]);
                    fis.add(InteractionUtilities.generateFIFromGene(tokens[1], tokens[2]));
                });
        System.out.println("Total interfaces: " + fis.size());
        System.out.println("Total cancer types: " + cancerTypes.size());

        // Generate the reaction sub-network
        System.out.println("\nGenerate the reaction network:");
        ReactionMapGenerator mapGenerator = new ReactionMapGenerator();
//        mapGenerator.generateSubNetwork(reactionIds);
        mapGenerator.generateSubNetworkForAll(reactionIds, new HashMap<>());
    }

    /**
     * Load the largest functional impact scores for mutations listed in cosmic between two proteins
     * interacting each other. The functional impact scores for mutations are reported in column 27
     * ( ie: interaction effect for all contacts of the site with this interactor).
     *
     * @return
     * @throws IOException
     */
    public Map<String, Double> loadPPIToMaxScore(String mechismoFileName) throws IOException {
        Path path = Paths.get(mechismoFileName);
        Comparator<Double> comparator = (v1, v2) -> Double.compare(Math.abs(v1), Math.abs(v2));
        // Use new Java 8 stream API
        Map<String, Optional<Double>> ppiToScore = Files.lines(path)
                .map(line -> line.split("\t"))
                .filter(tokens -> tokens.length >= 28)
                .filter(tokens -> tokens[19].trim().length() > 0)
                .collect(Collectors.groupingBy(tokens -> InteractionUtilities.generateFIFromGene(tokens[1], tokens[19]), // Key by PPI
                        Collectors.mapping(tokens -> new Double(tokens[27]), Collectors.maxBy(comparator)))); // Value should be maximum
        Map<String, Double> rtn = new HashMap<>();
        ppiToScore.forEach((ppi, score) -> {
            if (score.isPresent())
                rtn.put(ppi, score.get());
        });
        return rtn;

//        while ((line = fu.readLine()) != null) {
//            String[] tokens = line.split("\t");
//            if (tokens.length < 28)
//                continue;
//            // If no second protein is reported, escape it.
//            if (tokens[19].trim().length() == 0)
//                continue;
//            String ppi = InteractionUtilities.generateFIFromGene(tokens[1], 
//                                                                 tokens[19]);
//            Double score = ppiToScore.get(ppi);
//            Double currentScore = new Double(tokens[27]);
//            if (score == null || Math.abs(currentScore) > Math.abs(score))
//                ppiToScore.put(ppi, currentScore);
//        }
    }

    @Test
    public void testGetMechismoInteractions() throws IOException {
        Set<String> interactions = getInteractions(outputFileName);
        System.out.println("Total interactions: " + interactions.size());
        // Output is: 13323
    }

    @Test
    public void checkOverlapWithInteractome3d() throws IOException {
        Set<String> mechismoPPIs = getInteractions(outputFileName);
        System.out.println("Total PPIs in mechismo output: " + mechismoPPIs.size());
        fu.saveInteractions(mechismoPPIs, "results/MechismoPPIs_051017.txt");

        Interactome3dAnalyzer interactomeAnalyser = new Interactome3dAnalyzer();
        String interactomeDirName = "datasets/interactome3d/2016_06/prebuilt/representative/";
        Map<String, File> fiToPDB = interactomeAnalyser.loadPPIToPDBFile(interactomeDirName,
                false);
        System.out.println("Total PPIs in interactome3d: " + fiToPDB.size());
        fu.saveInteractions(fiToPDB.keySet(), "results/Interactome3dPPIs_051017.txt");

        Set<String> shared = InteractionUtilities.getShared(mechismoPPIs, fiToPDB.keySet());
        System.out.println("Shared: " + shared.size());
        // Output:
        //        Total PPIs in mechismo output: 13,323
        //        Total PPIs in interactome3d: 9,875
        //        Shared: 4011
        // TODO: The 50% coverage of interactome3d by mechismo may result from no cosmic mutations are found in other
        // PPIs. We should check the original cosmic mutation data and then compare the outputs.
    }

    /**
     * Get the interactions reported in the mechismo output file.
     *
     * @return
     * @throws IOException
     */
    public Set<String> getInteractions(String fileName) throws IOException {
        Set<String> interactions = new HashSet<>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 20)
                continue;
            if (tokens[19].trim().length() == 0)
                continue;
            String ppi = InteractionUtilities.generateFIFromGene(tokens[1], tokens[19]);
            interactions.add(ppi);
        }
        fu.close();
        return interactions;
    }

    /**
     * The mechismo scores are calculated for each mutation in a specific coordinate in a specific protein. It is
     * collected from scores for this specific mutation across all interactions invovled in the protein having
     * the mutation. So they are unique regard mutation, coordinate, and protein. However, scores used for cancer
     * analysis are not: the interacting partners need to be considered. In other words,
     * they are specific to mutations, coordinates, and two proteins involved in the interaction.
     *
     * @throws IOException
     */
    @Test
    public void checkMechismoScores() throws IOException {
        String fileName = dirName + "SOS1_RAS.txt";
        fileName = outputFileName;

        fu.setInput(fileName);
        String line = null;
        Map<String, String> mutationToScore = new HashMap<>();
        Set<String> scores = new HashSet<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");

            //            for (int i = 0; i < tokens.length; i++) {
            //                System.out.println(i + "\t" + tokens[i]);
            //            }
            //            if (true)
            //                break;

            if (tokens.length < 20) {
                //                System.out.println(line);
                continue;
            }
            // Make sure there are AA mutations
            if (tokens[6].length() > 0 && tokens[19].length() > 0) {
                //                mutationToScore.put(tokens[6] + "\t" + tokens[19], tokens[14]);
                mutationToScore.put(tokens[6], tokens[14]);
                scores.add(tokens[14]);
            }
        }
        System.out.println("Total mutations: " + mutationToScore.size()); // Should be 20 * 19 = 380 
        System.out.println("Total scores: " + new HashSet<>(mutationToScore.values()).size());
        System.out.println("Actual scores: " + scores.size());
        for (String mutation : mutationToScore.keySet())
            System.out.println(mutation + "\t" + mutationToScore.get(mutation));
        System.out.println("\nActual extra scores: ");
        scores.removeAll(mutationToScore.values());
        for (String score : scores)
            System.out.println(score);
    }

    @Test
    public void checkOutputFile() throws IOException {
        try (Stream<String> stream = Files.lines(Paths.get(outputFileName))) {
            String id_seq_a1 = "316";
            String id_hit = "293865";
            stream.forEach(line -> {
                String[] tokens = line.split("\t");
                if (tokens[2].equals(id_seq_a1))
                    System.out.println(line);
//                if (tokens.length > 31 && tokens[31].equals(id_hit))
//                    System.out.println(line);
            });
        }
    }

    @Test
    public void checkPCIContactFile() throws IOException {
//        String targetGene = "AKT1";
        String targetGene = "PTEN";
        String pdbId = null;
//        String pdbId = "2uvm";

        Map<String, String> idSeqA1ToGene = loadIdSeqA1ToGene();

        try (Stream<String> stream = Files.lines(Paths.get(pciContactFile))) {
            stream.skip(1)
                    .filter(line -> {
                        String[] tokens = line.split("\t");
                        String gene = idSeqA1ToGene.get(tokens[0]);
                        if (pdbId == null)
                            return targetGene.equals(gene);
                        else
                            return targetGene.equals(gene) && tokens[2].equals(pdbId);
                    })
                    .forEach(System.out::println);
        }
    }

    private Map<String, String> loadIdSeqA1ToGene() throws IOException {
        try (Stream<String> stream = Files.lines(Paths.get(outputFileName))) {
            Map<String, String> idToGene = new HashMap<>();
            stream.map(line -> line.split("\t"))
                    .forEach(tokens -> idToGene.put(tokens[2], tokens[0]));
            return idToGene;
        }
    }

    @Test
    public void pickRows() throws IOException {
        boolean needChemical = false;
        String[] targetGenes = new String[]{
                "HRAS",
//                "NRAS", 
//                "KRAS"
//                , "SOS1"
//                "PTEN"
//        		"AKT1"
        };
        Set<String> geneSet = new HashSet<>();
        for (String gene : targetGenes)
            geneSet.add(gene);

        String fileName = outputFileName;
        //        fileName = dirName + "HRAS_mechismo_op.tsv";
        try (Stream<String> stream = Files.lines(Paths.get(fileName))) {
            stream.map(line -> line.split("\t"))
                    .filter(tokens -> tokens.length >= 19)
                    .forEach(tokens -> {
                        if (geneSet.contains(tokens[0])) {
                            if (!needChemical ||
                                    needChemical && tokens[18].contains("CHEM"))
                                System.out.println(String.join("\t", tokens));
                        }
                    });
        }
//        fu.setInput(fileName);
//        String line = null;
//        while ((line = fu.readLine()) != null) {
//            String[] tokens = line.split("\t");
//            if (tokens.length < 19)
//                continue; // Cannot get interactors
//            //            System.out.println(tokens[0] + "\t" + tokens[18]);
////            if (geneSet.contains(tokens[0]) && geneSet.contains(tokens[18]))
////                System.out.println(line);
//            if (geneSet.contains(tokens[0]))
//            		System.out.println();
//        }
//        fu.close();
    }

    public void mapReactomeReactions(CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer,
                                     String mechismoOutputFilePath,
                                     String reactomeReactionNetworkFilePath,
                                     String mechismoSamples2SigReactionsFilePath,
                                     String mechismoFIFilterFilePath,
                                     String outputDir,
                                     String outputFilePrefix,
                                     String tcgaCancerType,
                                     int depth,
                                     int numPermutations,
                                     double mechScoreLowerBoundInclusive,
                                     double eThresh,
                                     boolean ignoreDependentUpstreamReactions,
                                     boolean ignoreIndependentUpstreamReactions,
                                     boolean excludeMultipleImmediateUpstreamReactions,
                                     boolean rewireLargestComponentOnly,
                                     boolean useRxnDist,
                                     String rxnFilter,
                                     Integer minNumTargetRxnPatients) throws Exception {
        ResourceMonitor resourceMonitor = new ResourceMonitor();
        resourceMonitor.StartMethodTimer();

        System.out.println("Loading reactome data...");
        ReactomeReactionGraphLoader reactomeReactionGraphLoader =
                new ReactomeReactionGraphLoader(
                        reactomeReactionNetworkFilePath,
                        mechismoSamples2SigReactionsFilePath,
                        rxnFilter);

        System.out.println("Loading mechismo data...");
        MechismoOutputLoader mechismoOutputLoader =
                new MechismoOutputLoader(mechismoOutputFilePath,
                        mechismoFIFilterFilePath,
                        tcgaCancerType,
                        mechScoreLowerBoundInclusive,
                        eThresh);

        System.out.println("Mapping reactome & mechismo data...");
        ReactomeMechismoDataMap reactomeMechismoDataMap =
                new ReactomeMechismoDataMap(
                        cancerDriverReactomeAnalyzer,
                        mechismoOutputLoader,
                        reactomeReactionGraphLoader);

        List<CooccurrenceResult> rewiredNetworkResults = new ArrayList<>();

        System.out.println("Creating reaction graph...");
        DefaultDirectedGraph<Long, DefaultEdge> reactionGraph =
                reactomeReactionGraphLoader.getReactionGraph();

        ReactionGraphAnalyzer reactionGraphAnalyzer = new ReactionGraphAnalyzer(
                depth,
                reactomeMechismoDataMap,
                ignoreDependentUpstreamReactions,
                ignoreIndependentUpstreamReactions,
                excludeMultipleImmediateUpstreamReactions);

        RandomGraphGenerator randomGraphGenerator = new RandomGraphGenerator(reactionGraph);

        for (int i = 0; i < numPermutations; i++) {
            resourceMonitor.StartLoopTimer();

            Graph<Long, DefaultEdge> rewiredReactionGraph =
                    randomGraphGenerator.GenerateRandomGraph(rewireLargestComponentOnly);

            ConnectivityInspector<Long, DefaultEdge> connectivityInspector =
                    new ConnectivityInspector<>(rewiredReactionGraph);

            System.out.println(String.format("Created reaction graph permutation with %d vertices, %d edges, %d components",
                    rewiredReactionGraph.vertexSet().size(),
                    rewiredReactionGraph.edgeSet().size(),
                    connectivityInspector.connectedSets().size()));

            CooccurrenceResult rewiredNetworkResult =
                    reactionGraphAnalyzer.SearchRxnNetworkAndCalculateCooccurrencePValues(
                            rewiredReactionGraph,
                            useRxnDist,
                            minNumTargetRxnPatients);

            rewiredNetworkResult.CalculateBHAdjustedPValues();
            //rewiredNetworkResult.writeToFile(outputDir, "RandomRewiring_" + (i + 1) + "_");

            resourceMonitor.CalculateMemUsed();

            rewiredNetworkResult.MagicallyShrinkMemoryFootprint();
            rewiredNetworkResults.add(rewiredNetworkResult);

            resourceMonitor.EndLoopTimer(i + 1 + "");//start with iteration '1'
        }

        CooccurrenceResult realResult =
                reactionGraphAnalyzer.SearchRxnNetworkAndCalculateCooccurrencePValues(
                        reactionGraph,
                        useRxnDist,
                        minNumTargetRxnPatients);

        realResult.CalculateBHAdjustedPValues();
        realResult.CalculateEmpiricalPValues(rewiredNetworkResults);

        reactomeMechismoDataMap.WriteRxn2SamplesToFile(outputDir, "");

        realResult.writeToFile(outputDir, outputFilePrefix);
        realResult.writePatientGroupingsToFile(outputDir, outputFilePrefix);
        realResult.writePatientDistancesToFile(outputDir,
                outputFilePrefix,
                reactomeMechismoDataMap);

        resourceMonitor.CalculateMemUsed();
        resourceMonitor.EndMethodTimer();
    }

    private Map<String, Set<String>> transformNetworkToLineGraph(Map<String, Set<String>> network) {
        Map<String, Set<String>> lineGraph = new HashMap<>();

        for (String commonGene : network.keySet()) {
            Set<String> partnerGenes = network.get(commonGene);
            Set<String> fis = new HashSet<>();
            for (String partnerGene : partnerGenes) {
                fis.add(FI.convertGeneNamePairToFIName("\t", commonGene, partnerGene));
            }
            for (String fi : fis) {
                Set<String> partnerFIs = new HashSet<>(fis);
                partnerFIs.remove(fi);
                Set<String> existingPartners;
                if (lineGraph.containsKey(fi)) {
                    existingPartners = lineGraph.get(fi);
                } else {
                    existingPartners = new HashSet<>();
                }
                partnerFIs.addAll(existingPartners);
                lineGraph.put(fi, existingPartners);
            }
        }

        return lineGraph;
    }

    private Map<String, Map<String,Set<String>>> transformNetworkToLineGraphNamedEdges(Map<String, Set<String>> network) {
        Map<String, Map<String,Set<String>>> lineGraph = new HashMap<>();

        for (String commonGene : network.keySet()) {
            Set<String> partnerGenes = network.get(commonGene);
            Set<String> fis = new HashSet<>();
            for (String partnerGene : partnerGenes) {
                fis.add(FI.convertGeneNamePairToFIName("-", commonGene, partnerGene));
            }
            for (String fi : fis) {
                Set<String> partnerFIs = new HashSet<>(fis);
                partnerFIs.remove(fi);
                Map<String,Set<String>> existingPartners;
                if (lineGraph.containsKey(fi)) {
                    existingPartners = lineGraph.get(fi);
                } else {
                    existingPartners = new HashMap<>();
                }
                if(existingPartners.containsKey(commonGene)) {
                    partnerFIs.addAll(existingPartners.get(commonGene));
                }
                existingPartners.put(commonGene,partnerFIs);
                lineGraph.put(fi, existingPartners);
            }
        }

        return lineGraph;
    }

    public Map<FI,Double> parseSignificantFIs(String tcgaCancerType, String significantFIsFilePath, MySQLAdaptor dba) throws Exception {
        Map<String,String> geneToUniprotMap = new ReactomeAnalyzer().getGeneToUniprotMap(dba);
        Map<FI,Double> significantFIsToFDR = new HashMap<>();
        FileUtility fileUtility = new FileUtility();
        try {
            fileUtility.setInput(significantFIsFilePath);
            String line;
            String[] tokens;
            while ((line = fileUtility.readLine()) != null) {
                tokens = line.split("\t");
                if (tokens[0].trim().toUpperCase().equals(tcgaCancerType.trim().toUpperCase())) {
                    try {
                        Gene gene1 = new Gene(tokens[1], geneToUniprotMap.get(tokens[1]));
                        Gene gene2 = new Gene(tokens[2], geneToUniprotMap.get(tokens[2]));
                        Double fdr = new Double(tokens[11]);
                        significantFIsToFDR.put(new FI(gene1,gene2),fdr);
                    }catch(NullPointerException npe){
                        int debug = 1;
                        //ignore for now.. these are protein-chemical interactions
                    }
                }
            }
            fileUtility.close();
        } catch (IOException ioe) {
            logger.error(String.format("Couldn't use %s, %s: %s",
                    significantFIsFilePath,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
        return significantFIsToFDR;
    }

    public Set<String> parseSignificantFIs(String tcgaCancerType, String significantFIsFilePath) {
        Set<String> significantFIs = new HashSet<>();
        FileUtility fileUtility = new FileUtility();
        try {
            fileUtility.setInput(significantFIsFilePath);
            String line;
            String[] tokens;
            while ((line = fileUtility.readLine()) != null) {
                tokens = line.split("\t");
                if (tokens[0].trim().toUpperCase().equals(tcgaCancerType.trim().toUpperCase())) {
                    String geneName1 = tokens[1];
                    String geneName2 = tokens[2];
                    Set<String> fis;
                    significantFIs.add(FI.convertGeneNamePairToFIName("\t", geneName1, geneName2));
                }
            }
            fileUtility.close();
        } catch (IOException ioe) {
            logger.error(String.format("Couldn't use %s, %s: %s",
                    significantFIsFilePath,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
        return significantFIs;
    }

    public Set<String> findLargestComponent(Map<String, Set<String>> vtxToConnectedVtxs) {
        Set<String> vtxsInLargestComponent = new HashSet<>();
        Graph<String, DefaultEdge> pseudograph = new Pseudograph(DefaultEdge.class);

        int count = 0;
        for (String vtx : vtxToConnectedVtxs.keySet()) {
            if (!pseudograph.containsVertex(vtx)) {
                pseudograph.addVertex(vtx);
            }
            for (String connectedVtx : vtxToConnectedVtxs.get(vtx)) {
                if (!pseudograph.containsVertex(connectedVtx)) {
                    pseudograph.addVertex(connectedVtx);
                }
                if (!pseudograph.containsEdge(vtx, connectedVtx)) {
                    pseudograph.addEdge(vtx, connectedVtx);
                }
            }
            count++;
            if (count % 1000 == 0) {
                System.out.println(String.format("Searched %d of %d vertices...", count, vtxToConnectedVtxs.keySet().size()));
            }
        }

        ConnectivityInspector connectivityInspector = new ConnectivityInspector(pseudograph);
        List<Set<String>> components = connectivityInspector.connectedSets();
        for (Set<String> component : components) {
            if (component.size() > vtxsInLargestComponent.size()) {
                vtxsInLargestComponent = new HashSet<>(component);
            }
        }
        return vtxsInLargestComponent;
    }

    private boolean stringsMostlyMatch(String s1, String s2, boolean emptyAlwaysMatches) {
        if(emptyAlwaysMatches){
            if(s1.isEmpty() || s2.isEmpty()){
                return true;
            }
        }
        return s1.trim().toUpperCase().equals(s2.trim().toUpperCase());
    }

    private Map<String, Set<String>> loadPatientToFIs(String srcFileName,
                                                      String tcgaCancerType,
                                                      double fdrCutoff) throws IOException {

        int cancerTypeIdx = 0;
        int gene1Idx = 1;
        int gene2Idx = 2;
        int fdrIdx = 11;
        int patientSetIdx = 15;

        Map<String, Set<String>> patientToReactions = new HashMap<>();
        try (Stream<String> stream = Files.lines(Paths.get(srcFileName))) {
            stream
                    .map(line -> line.split("\t"))
                    .filter(tokens -> stringsMostlyMatch(tokens[cancerTypeIdx], tcgaCancerType, true)) // Filter to cancer type
                    .filter(tokens -> !tokens[gene2Idx].contains("[")) // exclude non-transcribed chemicals
                    .filter(tokens -> Double.parseDouble(tokens[fdrIdx]) < fdrCutoff) // Less than the predefined fdr cutoff
                    .forEach(tokens -> {
                        String fiName = FI.convertGeneNamePairToFIName("\t",
                                tokens[gene1Idx],
                                tokens[gene2Idx]);
                        int index1 = tokens[patientSetIdx].indexOf("\'");
                        int index2 = tokens[patientSetIdx].lastIndexOf("\'");
                        String tmp = tokens[patientSetIdx].substring(index1, index2 + 1);
                        String[] patients = tmp.split(", ");
                        Arrays.asList(patients).forEach(patient -> {
                            // Need to remove single quotes
                            patient = patient.substring(1, patient.length() - 1);
                            patientToReactions.compute(patient, (key, set) -> {
                                if (set == null)
                                    set = new HashSet<>();
                                set.add(fiName);
                                return set;
                            });
                        });
                    });
        }
        return patientToReactions;
    }

    private void generatePairwiseFINetworkDistance(Map<String, Set<String>> patientToSigFIsForCancerType,
                                                   Map<String, Set<String>> fiToFisSharingGene,
                                                   String outputFileName) throws IOException {
        // Calculate pair-wise average distances between two cancer types
        ReactionMapGenerator mapGenerator = new ReactionMapGenerator();
        BreadthFirstSearch bfs = new BreadthFirstSearch();

        // Need a little bit cleanup
        // ensure patient FIs are in FI network
        patientToSigFIsForCancerType.forEach((patient, sigFIsForCancerType) ->  sigFIsForCancerType.retainAll(fiToFisSharingGene.keySet()));
        filterSamplesWithoutReactions(patientToSigFIsForCancerType);
        List<String> patients = patientToSigFIsForCancerType.keySet().stream().sorted(Comparator.naturalOrder()).collect(Collectors.toList());

        // Generate output
        StringBuilder builder = new StringBuilder();
        builder.append("Cancer");
        patients.forEach(patient -> builder.append("\t").append(patient));
        fu.setOutput(outputFileName);
        fu.printLine(builder.toString());
        builder.setLength(0);

        // Use a helper object
        MechismoBFSHelper bfsHelper = new MechismoBFSHelper();
        //what distance is this?
        Map<String, Integer> fiPairToDist = bfsHelper.calculateShortestFIPath(patientToSigFIsForCancerType,
                fiToFisSharingGene,
                bfs);
        for (int i = 0; i < patients.size(); i++) {
            String patient1 = patients.get(i);
            builder.append(patient1);
            Set<String> patient1SigFIsForCancerType = patientToSigFIsForCancerType.get(patient1);
            for (int j = 0; j < patients.size(); j++) {
                String patient2 = patients.get(j);
                builder.append("\t");
                if (i == j)
                    builder.append(0.0d);
                else {
                    Set<String> patient2SigFIsForCancerType = patientToSigFIsForCancerType.get(patient2);
                    //what distance is this?
                    double dist = bfsHelper.calculateAvgShortestFIPath(
                            patient1SigFIsForCancerType,
                            patient2SigFIsForCancerType,
                            fiPairToDist);
                    builder.append(dist);
                }
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }

    public void calculateSampleDistancesOnFINetwork(String tcgaCancerType,
                                                    String fiNetworkFilePath,
                                                    String significantMechismoFIsFilePath,
                                                    double fdrCutoff,
                                                    String outputDir) throws IOException {

        Map<String, Set<String>> patientToFIs = loadPatientToFIs(
                significantMechismoFIsFilePath,
                tcgaCancerType,
                fdrCutoff);

        // A quick quality check
        patientToFIs.forEach((patient, set) -> System.out.println(patient + "\t" + set.size()));

        //this is a general method that should work for FI sets too
        filterSamplesWithoutReactions(patientToFIs);

        Set<String> geneNetwork =
                new ReactionMapGenerator().loadNetwork(fiNetworkFilePath, "\t", 1);

        Map<String, Set<String>> fiToFisSharingGeneMap = transformNetworkToLineGraph(
                new BreadthFirstSearch().generateIdToPartnersMap(geneNetwork));

        generatePairwiseFINetworkDistance(
                patientToFIs,
                fiToFisSharingGeneMap,
                String.format("%s/MechismoSamplePairwiseFINetworkDist_%s.tsv",
                        outputDir,
                        tcgaCancerType.isEmpty()
                                ? "Pancancer"
                                : tcgaCancerType));
    }

    private Set<Patient> parsePatientsFromClusterFile(String patientClusterFileName, ReactomeMechismoDataMap dataMap) throws IOException {
        Set<Patient> patients = new HashSet<>();
        FileUtility fileUtility0 = new FileUtility();
        fileUtility0.setInput(patientClusterFileName);
        String line;
        while((line = fileUtility0.readLine()) != null){
            Patient patient = dataMap.getPatientBarcodeToPatientMap(line.trim().toUpperCase());
            if(patient != null) {
                patients.add(patient);
            }else{
                int debug = 1;
                System.out.println(String.format("Can't find patient '%s'",line));
            }
        }
        fileUtility0.close();
        return patients;
    }

    private Map<Gene,Integer> generateGenePatientCountMap(Set<Patient> patients,
                                                          Map<FI,Double> significantFIsToFDR,
                                                          ReactomeMechismoDataMap dataMap){
        Map<Gene,Integer> genePatientCountMap = new HashMap<>();
        for(Patient patient : patients){
            if(patient == null){
                int debug = 1;
            }
            Set<FI> patientSigFIs = dataMap.getFIs(patient);
            patientSigFIs.retainAll(significantFIsToFDR.keySet());
            for(FI fi : patientSigFIs){
                Set<Gene> geneSet = new HashSet(fi.getGenes());
                for(Gene gene : geneSet){
                    Integer count = 1;
                    if(genePatientCountMap.containsKey(gene)){
                        count = 1 + genePatientCountMap.get(gene);
                    }
                    genePatientCountMap.put(gene,count);
                }
            }
        }
        return genePatientCountMap;
    }

    private Map<String,Integer> generateGeneStringPatientCountMap(Set<Patient> patients,
                                                          Map<FI,Double> significantFIsToFDR,
                                                          ReactomeMechismoDataMap dataMap){
        Map<String,Integer> genePatientCountMap = new HashMap<>();
        for(Patient patient : patients){
            if(patient == null){
                int debug = 1;
            }
            Set<FI> patientSigFIs = dataMap.getFIs(patient);
            patientSigFIs.retainAll(significantFIsToFDR.keySet());
            for(FI fi : patientSigFIs){
                Set<Gene> geneSet = new HashSet(fi.getGenes());
                for(Gene gene : geneSet){
                    Integer count = 1;
                    if(genePatientCountMap.containsKey(gene.getHgncName())){
                        count = 1 + genePatientCountMap.get(gene.getHgncName());
                    }
                    genePatientCountMap.put(gene.getHgncName(),count);
                }
            }
        }
        return genePatientCountMap;
    }

    private void writeMapToFile(Map<Gene,Integer> geneSampleCountMap, String filePath) throws IOException {
        FileUtility fileUtility = new FileUtility();
        fileUtility.setOutput(filePath);
        fileUtility.printLine("Gene\tSampleNumber");
        for(Gene gene : geneSampleCountMap.keySet()){
            fileUtility.printLine(String.format("%s\t%d",
                    gene.getHgncName(),
                    geneSampleCountMap.get(gene)));
        }
        fileUtility.close();
    }

    public void getFIsForPatientClusters(CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer,
                                         String mechismoOutputFilePath,
                                         String mechismoFIFilterFilePath,
                                         String tcgaCancerType,
                                         double mechScoreLowerBoundInclusive,
                                         double eThresh,
                                         String fiNetworkFilePath,
                                         String reactomeReactionNetworkFilePath,
                                         String significantMechismoFIsFilePath) throws Exception{
        System.out.println("Loading reactome data...");
        ReactomeReactionGraphLoader reactomeReactionGraphLoader =
                new ReactomeReactionGraphLoader(
                        reactomeReactionNetworkFilePath);

        System.out.println("Loading mechismo data...");
        MechismoOutputLoader mechismoOutputLoader =
                new MechismoOutputLoader(mechismoOutputFilePath,
                        mechismoFIFilterFilePath,
                        tcgaCancerType,
                        mechScoreLowerBoundInclusive,
                        eThresh);

        System.out.println("Mapping reactome & mechismo data...");
        ReactomeMechismoDataMap reactomeMechismoDataMap =
                new ReactomeMechismoDataMap(
                        cancerDriverReactomeAnalyzer,
                        mechismoOutputLoader,
                        reactomeReactionGraphLoader);

        System.out.println("Generating FI network...");
        Set<String> fiNetwork =
                new ReactionMapGenerator().loadNetwork(fiNetworkFilePath, "\t", 1);
        BreadthFirstSearch breadthFirstSearch = new BreadthFirstSearch();
        Map<String, Map<String,Set<String>>> fiToFisSharingGeneMap = transformNetworkToLineGraphNamedEdges(
                breadthFirstSearch.generateIdToPartnersMap(fiNetwork));

        //System.out.println("Finding largest component of line graph...");
        //Set<String> fisInLargestComponent = findLargestComponent(fiToFisSharingGeneMap);

        System.out.println(String.format("Parsing significant FIs for cancertype: %s",
                tcgaCancerType));
        Map<FI,Double> significantFIsForCancerType = parseSignificantFIs(tcgaCancerType,
                significantMechismoFIsFilePath,
                cancerDriverReactomeAnalyzer.getDBA());

        Set<Patient> cluster1Patients = parsePatientsFromClusterFile("/home/burkhart/UCEC_cluster1_samples.txt",
                reactomeMechismoDataMap);
        Set<Patient> cluster2Patients = parsePatientsFromClusterFile("/home/burkhart/UCEC_cluster2_samples.txt",
                reactomeMechismoDataMap);

        Map<Gene,Integer> cluster1Gene2PatientCount = generateGenePatientCountMap(cluster1Patients,
                significantFIsForCancerType,
                reactomeMechismoDataMap);
        Map<Gene,Integer> cluster2Gene2PatientCount = generateGenePatientCountMap(cluster2Patients,
                significantFIsForCancerType,
                reactomeMechismoDataMap);

        writeMapToFile(cluster1Gene2PatientCount, "/home/burkhart/UCEC_cluster1_geneSampleCounts.txt");
        writeMapToFile(cluster2Gene2PatientCount, "/home/burkhart/UCEC_cluster2_geneSampleCounts.txt");

        Map<String,Integer> cluster1GeneString2PatientCount = generateGeneStringPatientCountMap(cluster1Patients,
                significantFIsForCancerType,
                reactomeMechismoDataMap);
        Map<String,Integer> cluster2GeneString2PatientCount = generateGeneStringPatientCountMap(cluster2Patients,
                significantFIsForCancerType,
                reactomeMechismoDataMap);

        Map<String,Double> significantFIStringssForCancerType = new HashMap<>();
        for(FI fi : significantFIsForCancerType.keySet()){
            String fiString = fi.toString("-");
            significantFIStringssForCancerType.put(fiString,significantFIsForCancerType.get(fi));
        }

        FileUtility sif_fileUtility = new FileUtility();
        FileUtility tbl_fileUtility = new FileUtility();
        sif_fileUtility.setOutput("/home/burkhart/UCEC_fis.sif");
        tbl_fileUtility.setOutput("/home/burkhart/UCEC_fis.tbl");
        tbl_fileUtility.printLine(String.format("%s\t%s","cluster","name"));

        Set<String> fiEdgeSet = new HashSet<>();
        Set<String> c1NodeSet = new HashSet<>();
        Set<String> c2NodeSet = new HashSet<>();

        for(String fi : fiToFisSharingGeneMap.keySet()){
                for (String gene : fiToFisSharingGeneMap.get(fi).keySet()) {
                    if (cluster1GeneString2PatientCount.containsKey(gene) &&
                            cluster1GeneString2PatientCount.get(gene) > 30) {
                        for (String fi2 : fiToFisSharingGeneMap.get(fi).get(gene)) {
                            String edgeLineFwd = String.format("%s\t-\t%s", fi, fi2);
                            String edgeLineRev = String.format("%s\t-\t%s", fi2, fi);
                            if(!fiEdgeSet.contains(edgeLineFwd) && !fiEdgeSet.contains(edgeLineRev)) {
                                sif_fileUtility.printLine(edgeLineFwd);
                                fiEdgeSet.add(edgeLineFwd);
                                fiEdgeSet.add(edgeLineRev);
                            }
                            String c1fi1 = String.format("%s\t%s","Cluster1",fi);
                            if(!c1NodeSet.contains(c1fi1)) {
                                tbl_fileUtility.printLine(c1fi1);
                                c1NodeSet.add(c1fi1);
                            }
                            String c1fi2 = String.format("%s\t%s","Cluster1",fi2);
                            if(!c1NodeSet.contains(c1fi2)) {
                                tbl_fileUtility.printLine(c1fi2);
                                c1NodeSet.add(c1fi2);
                            }
                        }
                    }

                    if (cluster2GeneString2PatientCount.containsKey(gene) &&
                            cluster2GeneString2PatientCount.get(gene) > 30) {
                        for (String fi2 : fiToFisSharingGeneMap.get(fi).get(gene)) {
                            String edgeLineFwd = String.format("%s\t-\t%s", fi, fi2);
                            String edgeLineRev = String.format("%s\t-\t%s", fi2, fi);
                            if(!fiEdgeSet.contains(edgeLineFwd) && !fiEdgeSet.contains(edgeLineRev)) {
                                sif_fileUtility.printLine(edgeLineFwd);
                                fiEdgeSet.add(edgeLineFwd);
                                fiEdgeSet.add(edgeLineRev);
                            }
                            String c2fi1 = String.format("%s\t%s","Cluster2",fi);
                            if(!c2NodeSet.contains(c2fi1)) {
                                tbl_fileUtility.printLine(c2fi1);
                                c2NodeSet.add(c2fi1);
                            }
                            String c2fi2 = String.format("%s\t%s","Cluster2",fi2);
                            if(!c2NodeSet.contains(c2fi2)) {
                                tbl_fileUtility.printLine(c2fi2);
                                c2NodeSet.add(c2fi2);
                            }
                        }
                    }
            }
        }
        sif_fileUtility.close();
        tbl_fileUtility.close();
    }

    public void analyzeInterfaceCooccurrence(CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer,
                                             String mechismoOutputFilePath,
                                             String mechismoFIFilterFilePath,
                                             String tcgaCancerType,
                                             double mechScoreLowerBoundInclusive,
                                             double eThresh,
                                             String fiNetworkFilePath,
                                             String reactomeReactionNetworkFilePath,
                                             String significantMechismoFIsFilePath,
                                             String outputDir) throws Exception {
        System.out.println("Loading reactome data...");
        ReactomeReactionGraphLoader reactomeReactionGraphLoader =
                new ReactomeReactionGraphLoader(
                        reactomeReactionNetworkFilePath);

        System.out.println("Loading mechismo data...");
        MechismoOutputLoader mechismoOutputLoader =
                new MechismoOutputLoader(mechismoOutputFilePath,
                        mechismoFIFilterFilePath,
                        tcgaCancerType,
                        mechScoreLowerBoundInclusive,
                        eThresh);

        System.out.println("Mapping reactome & mechismo data...");
        ReactomeMechismoDataMap reactomeMechismoDataMap =
                new ReactomeMechismoDataMap(
                        cancerDriverReactomeAnalyzer,
                        mechismoOutputLoader,
                        reactomeReactionGraphLoader);

        System.out.println("Generating FI network...");
        Set<String> fiNetwork =
                new ReactionMapGenerator().loadNetwork(fiNetworkFilePath, "\t", 1);
        BreadthFirstSearch breadthFirstSearch = new BreadthFirstSearch();
        Map<String, Set<String>> fiToFisSharingGeneMap = transformNetworkToLineGraph(
                breadthFirstSearch.generateIdToPartnersMap(fiNetwork));

        //System.out.println("Finding largest component of line graph...");
        //Set<String> fisInLargestComponent = findLargestComponent(fiToFisSharingGeneMap);

        System.out.println(String.format("Parsing significant FIs for cancertype: %s",
                tcgaCancerType));
        Set<String> significantFIsForCancerType = parseSignificantFIs(tcgaCancerType, significantMechismoFIsFilePath);

        System.out.println("Calculating shortest paths between patient pairs...");
        List<Patient> patientList = new ArrayList<>(reactomeMechismoDataMap.getPatients());
        FileUtility fileUtility0 = new FileUtility();
        String outFilePath0 = outputDir + tcgaCancerType + "_PatientPairDistances.csv";
        try {
            fileUtility0.setOutput(outFilePath0);
            fileUtility0.printLine("Patient 1," +
                    "Patient 2," +
                    "Min Shortest FI Path," +
                    "Average FI Distance");
            for (int i = 0; i < patientList.size() - 1; i++) {
                for (int j = i + 1; j < patientList.size(); j++) {
                    Patient patient1 = patientList.get(i);
                    Patient patient2 = patientList.get(j);
                    Set<String> p1FIs = reactomeMechismoDataMap.getFIStrings(patient1);
                    Set<String> p2FIs = reactomeMechismoDataMap.getFIStrings(patient2);
                    //p1FIs.retainAll(fisInLargestComponent);
                    p1FIs.retainAll(significantFIsForCancerType);
                    //p2FIs.retainAll(fisInLargestComponent);
                    p2FIs.retainAll(significantFIsForCancerType);

                    if (fiToFisSharingGeneMap != null) {
                        try {
                            double minShortestPath =
                                    breadthFirstSearch.calculateMinShortestPath(p1FIs, p2FIs, fiToFisSharingGeneMap);
                            double avgDistance =
                                    breadthFirstSearch.calculateAverageDistance(p1FIs, p2FIs, fiToFisSharingGeneMap);

                            //findLargestComponent() takes too long.. doing this instead
                            if (minShortestPath > 0 && avgDistance > 0) {
                                //write to file
                                fileUtility0.printLine(String.format("%s,%s,%.3f,%.3f",
                                        patient1,
                                        patient2,
                                        minShortestPath,
                                        avgDistance));
                            }
                        } catch (NullPointerException npe) {
                            //do nothing, sometimes thrown by breadthFirstSearch for unknown reason
                        }
                    }
                }
                System.out.println(String.format("Processed %d of %d patient pair distances...",
                        (i + 1) * patientList.size() - ((i + 1) * (i + 1)) / 2,
                        (patientList.size() * patientList.size()) / 2));
            }
            fileUtility0.close();
        } catch (IOException ioe) {
            System.out.println(String.format("Couldn't use %s, %s: %s",
                    outFilePath0,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }

        System.out.println("Accumulating FI pair cooccurrences & p-values...");
        List<Set<FI>> fiPairs = new ArrayList<>();
        Set<Set<FI>> fiPairsSet = new HashSet<>();
        List<Double> pvalues = new ArrayList<>();

        List<Integer> avalues = new ArrayList<>();
        List<Integer> bvalues = new ArrayList<>();
        List<Integer> cvalues = new ArrayList<>();
        List<Integer> dvalues = new ArrayList<>();

        Integer nPatients = reactomeMechismoDataMap.getNumPatients();

        FisherExact fisherExact = new FisherExact(nPatients);

        int x = 0;

        for (FI fi : reactomeMechismoDataMap.getFIs()) {
            Set<Patient> patientsAffectedByFirstFI = reactomeMechismoDataMap.getPatients(fi);
            Set<FI> cooccurringFIs = new HashSet<>();
            for (Patient patient : patientsAffectedByFirstFI) {
                cooccurringFIs.addAll(reactomeMechismoDataMap.getFIs(patient));
            }
            //find ABCDs -- this will find every pair multiple times... check hash first
            for (FI cooccurringFI : cooccurringFIs) {
                Set<Gene> pairGenes = new HashSet<>();
                pairGenes.addAll(fi.getGenes());
                pairGenes.addAll(cooccurringFI.getGenes());
                if (pairGenes.size() == 4) { // ensure FI pairs don't share genes
                    Set<FI> fiPair = new HashSet<>();
                    fiPair.add(fi);
                    fiPair.add(cooccurringFI);
                    if (!fiPairsSet.contains(fiPair)) {
                        Set<Patient> patientsAffectedBySecondFI =
                                reactomeMechismoDataMap.getPatients(cooccurringFI);

                        Set<Patient> patientsAffectedByBothFIs = new HashSet<>(patientsAffectedByFirstFI);
                        patientsAffectedByBothFIs.retainAll(patientsAffectedBySecondFI);
                        Integer D = patientsAffectedByBothFIs.size();

                        if (D > 1) {
                            Set<Patient> patientsAffectedBySecondFIOnly = new HashSet<>(patientsAffectedBySecondFI);
                            patientsAffectedBySecondFIOnly.removeAll(patientsAffectedByBothFIs);
                            Integer C = patientsAffectedBySecondFIOnly.size();

                            Set<Patient> patientsAffectedByFirstFIOnly = new HashSet<>(patientsAffectedByFirstFI);
                            patientsAffectedByFirstFIOnly.removeAll(patientsAffectedByBothFIs);
                            Integer B = patientsAffectedByFirstFIOnly.size();

                            Integer A = nPatients - (B + C + D);

                            Double p = fisherExact.getRightTailedP(A, B, C, D);
                            avalues.add(A);
                            bvalues.add(B);
                            cvalues.add(C);
                            dvalues.add(D);
                            fiPairs.add(fiPair);
                            fiPairsSet.add(fiPair);
                            pvalues.add(p);
                        }
                    }
                }
            }
            x++;
            if (x % 1000 == 0) {
                System.out.println("processed " + x + " of " + reactomeMechismoDataMap.getFIs().size() + " fis...");
            }
        }
        //calculate FDR
        Map<Double, Double> pvalueFDRMap = new HashMap<>();
        List<Double> pvaluesSorted = new ArrayList<>(pvalues);
        List<Double> fdrs = MathUtilities.calculateFDRWithBenjaminiHochberg(pvaluesSorted);
        for (int i = 0; i < fdrs.size(); i++) {
            pvalueFDRMap.put(pvaluesSorted.get(i), fdrs.get(i));
        }
        FileUtility fileUtility = new FileUtility();
        String outFilePath = outputDir + tcgaCancerType + "_InterfaceCooccurrence.csv";
        try {
            fileUtility.setOutput(outFilePath);
            fileUtility.printLine("FI Pair," +
                    "#Samples With Neither Interface Mutated," +
                    "#Samples With First Interface Mutated," +
                    "#Samples With Second Interface Mutated," +
                    "#Samples With Both Interfaces Mutated," +
                    "Fisher Exact P-value," +
                    "BH Adjusted P-value");
            for (int i = 0; i < fiPairs.size(); i++) {
                List<FI> fiPair = new ArrayList<>(fiPairs.get(i));
                Collections.sort(fiPair);
                fileUtility.printLine(String.format("%s | %s,%d,%d,%d,%d,%.100e,%.100e",
                        fiPair.get(0),
                        fiPair.get(1),
                        avalues.get(i),
                        bvalues.get(i),
                        cvalues.get(i),
                        dvalues.get(i),
                        pvalues.get(i),
                        pvalueFDRMap.get(pvalues.get(i))));
            }
            fileUtility.close();
        } catch (IOException ioe) {
            System.out.println(String.format("Couldn't use %s, %s: %s",
                    outFilePath,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
        System.out.println("done.");
    }

    private Map<Patient, Set<FI>> calculatePatientGraphFIs(Pseudograph<FI, DefaultEdge> edgeticFIGraph,
                                                           ReactomeMechismoDataMap reactomeMechismoDataMap) {
        Map<Patient, Set<FI>> patientGraphCenters = new HashMap<>();
        for (Patient patient : reactomeMechismoDataMap.getPatients()) {
            Set<FI> patientFIs = new HashSet<>(reactomeMechismoDataMap.getFIs(patient));
            patientFIs.retainAll(edgeticFIGraph.vertexSet());
            AsSubgraph<FI, DefaultEdge> patientSubgraph = new AsSubgraph<>(edgeticFIGraph,
                    patientFIs);
            GraphMeasurer graphMeasurer = new GraphMeasurer(patientSubgraph);
            /*if(!patientFIs.isEmpty()) {
                List<FI> fiList = new ArrayList<>(patientFIs);
                int half = (fiList.size() / 2) + 1;
                Set<FI> halfFISet = new HashSet<>();
                for (int i = 0; i < half; i++) {
                    halfFISet.add(fiList.get(i));
                }
                patientGraphCenters.put(patient, halfFISet);
            }*/
            patientGraphCenters.put(patient, graphMeasurer.getGraphCenter());
            System.out.println(String.format("Calculated graph center for %d of %d patients...",
                    patientGraphCenters.keySet().size(),
                    reactomeMechismoDataMap.getPatients().size()));
        }
        return patientGraphCenters;
    }

    private Double FindAverageShortestPathLengthBetweenPatients(Set<FI> patient1GraphCenter,
                                                                Set<FI> patient2GraphCenter,
                                                                ConnectivityInspector connectivityInspector,
                                                                AStarShortestPath aStarShortestPath) {
        double pathLengthSum = 0.0d;
        int pathCount = 0;

        for (FI patient1FI : patient1GraphCenter) {
            for (FI patient2FI : patient2GraphCenter) {
                if (!patient1FI.equals(patient2FI) &&
                        connectivityInspector.pathExists(patient1FI, patient2FI)) {
                    pathCount++;
                    try {
                        GraphPath graphPath = aStarShortestPath.getPath(
                                patient1FI,
                                patient2FI);
                        if (graphPath != null) {
                            pathLengthSum +=
                                    (double) graphPath.getLength();
                        } else {
                            int debug = 1;
                        }
                    } catch (NoSuchMethodError nme) {
                        int debug = 1;
                    }
                    //System.out.println(String.format("Calculated %d of %d possible paths between patients",
                    //        pathCount,
                    //        patient1GraphCenter.size() * patient2GraphCenter.size()));
                }
            }
        }
        pathLengthSum = pathCount > 0
                ? pathLengthSum / (double) pathCount
                : (double) (patient1GraphCenter.size() * patient2GraphCenter.size());

        return pathLengthSum;
    }
}
