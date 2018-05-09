package org.reactome.mechismo;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math.stat.correlation.SpearmansCorrelation;
import org.gk.model.GKInstance;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.annotate.AnnotationHelper;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.cancer.driver.CancerDriverReactomeAnalyzer;
import org.reactome.cancer.driver.Interactome3dDriverAnalyzer;
import org.reactome.cancer.driver.MechismoBFSHelper;
import org.reactome.r3.ReactionMapGenerator;
import org.reactome.r3.graph.BreadthFirstSearch;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.FisherExact;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;

/**
 * This class is used to perform the mechismo analysis in the perspective of Reactome reactions.
 * @author wug
 *
 */
public class ReactionNetworkAnalyzer {
    private String dirName = "datasets/Mechismo/";
    private FileUtility fu = new FileUtility();
    
    public ReactionNetworkAnalyzer() {
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
}
