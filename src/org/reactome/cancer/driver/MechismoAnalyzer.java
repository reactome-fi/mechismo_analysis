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
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.junit.Test;
import org.reactome.annotate.AnnotationHelper;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.cancer.*;
import org.reactome.px.util.InteractionUtilities;
import org.reactome.r3.Interactome3dAnalyzer;
import org.reactome.r3.ReactomeAnalyzer;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.FisherExact;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.Plotter;

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
        String output = "results/ReactionsInMechisom_051017.txt";
        output = "results/ReactionsInMechisom_052217.txt";
        CancerDriverReactomeAnalyzer reactomeAnalyzer = new CancerDriverReactomeAnalyzer();

        checkAllHumanReactions(reactomeAnalyzer, output);
    }

    /**
     * Check all human reactions based on mechismo data.
     *
     * @throws Exception
     */
    public void checkAllHumanReactions(CancerDriverReactomeAnalyzer reactomeAnalyzer,
                                       String outputFileName) throws Exception {
        // Load all the non-disease reactions from Reactome
        List<GKInstance> reactions = reactomeAnalyzer.loadHumanReactions();
        System.out.println("Total human reactions in Reactome: " + reactions.size());

        Map<String, Double> ppiToScore = loadPPIToMaxScore();
        System.out.println("Total PPIs in mechismo: " + ppiToScore.size());

        ReactomeAnalyzer reactomeDataAnalyzer = new ReactomeAnalyzer();
        fu.setOutput(outputFileName);
        fu.printLine("DB_ID\tName\tTotal_FIs\tMechismo_FIs\tMax_Score");
        int total = 0;
        int checkedReactions = 0;
        for (GKInstance reaction : reactions) {
            Set<String> fis = reactomeDataAnalyzer.generateTentativePPIsForReaction(reaction,
                    false);
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
                }
            }
            fu.printLine(reaction.getDBID() + "\t" +
                    reaction.getDisplayName() + "\t" +
                    fis.size() + "\t" +
                    overlap + "\t" +
                    (score == null ? "" : score));
            if (overlap > 0)
                total++;
            checkedReactions++;
        }
        fu.close();
        System.out.println("Total checked reactions: " + checkedReactions);
        System.out.println("Total reactions having FIs in mechsimo: " + total);
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
        String targetGene = "AKT1";
        String pdbId = "2uvm";

        Map<String, String> idSeqA1ToGene = loadIdSeqA1ToGene();

        try (Stream<String> stream = Files.lines(Paths.get(pciContactFile))) {
            stream.skip(1)
                    .filter(line -> {
                        String[] tokens = line.split("\t");
                        String gene = idSeqA1ToGene.get(tokens[0]);
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
        String[] targetGenes = new String[]{
                "HRAS", "NRAS", "KRAS", "SOS1"
        };
        Set<String> geneSet = new HashSet<>();
        for (String gene : targetGenes)
            geneSet.add(gene);

        String fileName = outputFileName;
        //        fileName = dirName + "HRAS_mechismo_op.tsv";
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 19)
                continue; // Cannot get interactors
            //            System.out.println(tokens[0] + "\t" + tokens[18]);
            if (geneSet.contains(tokens[0]) && geneSet.contains(tokens[18]))
                System.out.println(line);
        }
        fu.close();
    }

    private Set<Reaction[]> GenerateReactionSetPairs(Set<Reaction> reactionSet) {
        Set<Reaction[]> reactionSetPairs = new HashSet<>();
        List<Reaction> reactionList = new ArrayList<>(reactionSet);
        for (int i = 0; i < reactionList.size() - 1; i++) { //first to second-last
            for (int j = i + 1; j < reactionList.size(); j++) { //second to last
                if (j == i) {
                    continue;
                }
                Reaction[] pair = {reactionList.get(i), reactionList.get(j)};
                reactionSetPairs.add(pair);
            }
        }
        return reactionSetPairs;
    }

    private DirectedGraph<Long, DefaultEdge> RewireReactionGraph(
            DefaultDirectedGraph<Long, DefaultEdge> reactionGraph,
            Random prng) {
        DefaultDirectedGraph<Long, DefaultEdge> rewiredReactionGraph = new DefaultDirectedGraph<>(DefaultEdge.class);
        Graphs.addGraph(rewiredReactionGraph, reactionGraph);

        RewireComponents(rewiredReactionGraph, prng);

        return rewiredReactionGraph;
    }

    private void RewireComponents(DefaultDirectedGraph<Long, DefaultEdge> rewiredReactionGraph,
                                  Random prng) {
        ConnectivityInspector connectivityInspector = new ConnectivityInspector(rewiredReactionGraph);
        List<Set<Long>> connectedSets = connectivityInspector.connectedSets();

        int initialComponentCount = connectedSets.size();
        int initialEdgeCount = rewiredReactionGraph.edgeSet().size();
        int initialVtxCount = rewiredReactionGraph.vertexSet().size();

        for (Set<Long> connectedSet : connectedSets) {
            if (connectedSet.size() >= 3) {

                if (connectedSet.size() >= initialVtxCount) {
                    throw new IllegalStateException(String.format(
                            "component vtx count (%d) should be less than initial (%d)",
                            connectedSet.size(),
                            initialVtxCount));
                }

                Set<DefaultEdge> componentEdges = new HashSet<>();
                for (Long vertex : connectedSet) {
                    componentEdges.addAll(rewiredReactionGraph.incomingEdgesOf(vertex));
                    componentEdges.addAll(rewiredReactionGraph.outgoingEdgesOf(vertex));
                }

                if (componentEdges.size() >= initialEdgeCount) {
                    throw new IllegalStateException(String.format(
                            "component edge count (%d) should be less than initial (%d)",
                            componentEdges.size(),
                            initialEdgeCount));
                }


                int E = componentEdges.size();
                int V = connectedSet.size();
                DefaultEdge randEdge1, randEdge2;
                ArrayList<DefaultEdge> componentEdgesAry = new ArrayList<>(componentEdges);
                ArrayList<Long> componentVtxAry = new ArrayList<>(connectedSet);
                Long sourceVtx1, sourceVtx2, targetVtx1, targetVtx2;

                //from https://en.wikipedia.org/wiki/Degree-preserving_randomization
                int rewires = (int) ((E / 2.0) * Math.log(1 / Math.pow(10, -7)));
                for (int i = 0; i < rewires; i++) {

                    randEdge1 = componentEdgesAry.get(prng.nextInt(E));
                    randEdge2 = componentEdgesAry.get(prng.nextInt(E));

                    while (Objects.equals(randEdge1, randEdge2)) {
                        randEdge1 = componentEdgesAry.get(prng.nextInt(E));
                        randEdge2 = componentEdgesAry.get(prng.nextInt(E));
                    }

                    rewiredReactionGraph.removeEdge(randEdge1);
                    componentEdgesAry.remove(randEdge1);

                    int curEdgeCount1 = rewiredReactionGraph.edgeSet().size();
                    if (curEdgeCount1 != (initialEdgeCount - 1)) {
                        throw new IllegalStateException(String.format(
                                "current edge count (%d) should not differ from initial - 1 (%d)",
                                curEdgeCount1,
                                initialEdgeCount - 1));
                    }

                    rewiredReactionGraph.removeEdge(randEdge2);
                    componentEdgesAry.remove(randEdge2);

                    int curEdgeCount2 = rewiredReactionGraph.edgeSet().size();
                    if (curEdgeCount2 != (initialEdgeCount - 2)) {
                        throw new IllegalStateException(String.format(
                                "current edge count (%d) should not differ from initial - 2 (%d)",
                                curEdgeCount2,
                                initialEdgeCount - 2));
                    }

                    sourceVtx1 = componentVtxAry.get(prng.nextInt(V));
                    targetVtx1 = componentVtxAry.get(prng.nextInt(V));

                    sourceVtx2 = componentVtxAry.get(prng.nextInt(V));
                    targetVtx2 = componentVtxAry.get(prng.nextInt(V));

                    while (Objects.equals(sourceVtx1, sourceVtx2) ||
                            Objects.equals(targetVtx1, targetVtx2) ||
                            Objects.equals(sourceVtx1, targetVtx1) ||
                            Objects.equals(sourceVtx2, targetVtx2)) {
                        sourceVtx1 = componentVtxAry.get(prng.nextInt(V));
                        targetVtx1 = componentVtxAry.get(prng.nextInt(V));

                        sourceVtx2 = componentVtxAry.get(prng.nextInt(V));
                        targetVtx2 = componentVtxAry.get(prng.nextInt(V));
                    }

                    //randomly flip edge direction
                    //rewire edge 1
                    WireEdge(V,
                            prng,
                            sourceVtx1,
                            targetVtx2,
                            rewiredReactionGraph,
                            componentEdgesAry,
                            componentVtxAry);

                    int curEdgeCount3 = rewiredReactionGraph.edgeSet().size();
                    if (curEdgeCount3 != (initialEdgeCount - 1)) {
                        throw new IllegalStateException(String.format(
                                "current edge count (%d) should not differ from initial - 1 (%d)",
                                curEdgeCount3,
                                initialEdgeCount - 1));
                    }

                    //rewire edge 2
                    WireEdge(V,
                            prng,
                            sourceVtx2,
                            targetVtx1,
                            rewiredReactionGraph,
                            componentEdgesAry,
                            componentVtxAry);

                    int curEdgeCount4 = rewiredReactionGraph.edgeSet().size();
                    if (curEdgeCount4 != initialEdgeCount) {
                        throw new IllegalStateException(String.format(
                                "current edge count (%d) should not differ from initial (%d)",
                                curEdgeCount4,
                                initialEdgeCount));
                    }
                }
            }
        }
        int finalComponentCount = connectivityInspector.connectedSets().size();
        int finalEdgeCount = rewiredReactionGraph.edgeSet().size();
        int finalVtxCount = rewiredReactionGraph.vertexSet().size();

        if (initialComponentCount > finalComponentCount) {
            throw new IllegalStateException(String.format(
                    "Final component count (%d) should be no greater than initial (%d)",
                    finalComponentCount,
                    initialComponentCount));
        }

        if (initialEdgeCount != finalEdgeCount) {
            throw new IllegalStateException(String.format(
                    "Final edge count (%d) should not differ from initial (%d)",
                    finalEdgeCount,
                    initialEdgeCount));
        }

        if (initialVtxCount != finalVtxCount) {
            throw new IllegalStateException(String.format(
                    "Final vertex count (%d) should not differ from initial (%d)",
                    finalVtxCount,
                    initialVtxCount));
        }
    }

    private void WireEdge(int V,
                          Random prng,
                          Long sourceVtx,
                          Long targetVtx,
                          DefaultDirectedGraph<Long, DefaultEdge> graph,
                          ArrayList<DefaultEdge> edgeAry,
                          ArrayList<Long> vtxAry) {
        int rand;
        while (graph.getEdge(sourceVtx, targetVtx) != null) {
            rand = prng.nextInt(V);
            sourceVtx = vtxAry.get(rand);
            rand = prng.nextInt(V);
            targetVtx = vtxAry.get(rand);
        }
        graph.addEdge(sourceVtx, targetVtx);
        edgeAry.add(
                graph.getEdge(sourceVtx, targetVtx));
    }

    private void SearchUpstreamReactions(int depth,
                                         DirectedGraph<Long, DefaultEdge> reactionGraph,
                                         Set<DefaultEdge> incomingEdges,
                                         Reaction targetReactionCandidate,
                                         Set<TargetReactionCandidate> targetReactionCandidates,
                                         ReactomeMechismoDataMap reactomeMechismoDataMap) {
        Set<DefaultEdge> incomingEdgesCpy = new HashSet<>(incomingEdges);
        Set<Reaction> upstreamReactions = new HashSet<>();
        Set<Reaction> allUpstreamReactions = new HashSet<>();
        List<Set<Reaction>> upstreamReactionsByDepth = new ArrayList<>();
        for (int i = 1; i <= depth; i++) {
            for (DefaultEdge incomingEdge : incomingEdgesCpy) {
                upstreamReactions.add(reactomeMechismoDataMap.getReaction(reactionGraph.getEdgeSource(incomingEdge)));
            }
            upstreamReactions.removeAll(allUpstreamReactions);
            allUpstreamReactions.addAll(upstreamReactions);
            if (i < depth) {
                upstreamReactionsByDepth.add(new HashSet<>(upstreamReactions));
                upstreamReactions.clear();
                incomingEdgesCpy.clear();
                for (Reaction upstreamReaction : upstreamReactionsByDepth.get(i - 1)) {
                    incomingEdgesCpy.addAll(reactionGraph.incomingEdgesOf(upstreamReaction.getReactionID()));
                }
            }
        }
        incomingEdgesCpy.clear();
        upstreamReactions.clear();
        upstreamReactionsByDepth.clear();

        allUpstreamReactions.remove(targetReactionCandidate);

        //find intersection of upstream reactions and supported reactions
        Set<Reaction> supportedUpstreamReactions = new HashSet<>(allUpstreamReactions);
        supportedUpstreamReactions.retainAll(reactomeMechismoDataMap.getSupportedReactions());

        //remove upstream reactions directly activated by target
        supportedUpstreamReactions.removeIf(
                upstreamReaction -> !reactionGraph.getAllEdges(
                        targetReactionCandidate.getReactionID(),
                        upstreamReaction.getReactionID()).isEmpty());

        Set<FI> supportedFIs = new HashSet<>();
        for (Reaction supportedUpstreamReaction : supportedUpstreamReactions) {
            supportedFIs.addAll(reactomeMechismoDataMap.getFIs(supportedUpstreamReaction));
        }

        targetReactionCandidates.add(new TargetReactionCandidate(
                targetReactionCandidate,
                supportedUpstreamReactions,
                supportedFIs));
    }

    private void SearchRxnNetworkForTargetReactionCandidates(Set<TargetReactionCandidate> targetReactionCandidates,
                                                             DirectedGraph<Long, DefaultEdge> reactionGraph,
                                                             ReactomeMechismoDataMap reactomeMechismoDataMap,
                                                             int depth) {
        Set<DefaultEdge> outgoingEdges = new HashSet<>();
        Set<Reaction> targetReactionCandidateSet = new HashSet<>();

        for (Reaction reaction : reactomeMechismoDataMap.getSupportedReactions()) {
            //http://jgrapht.org/javadoc/org/jgrapht/graph/DefaultDirectedGraph.html
            //go one reaction downstream
            outgoingEdges.addAll(reactionGraph.outgoingEdgesOf(reaction.getReactionID()));
        }
        for (DefaultEdge outgoingEdge : outgoingEdges) {
            targetReactionCandidateSet.add(
                    reactomeMechismoDataMap.getReaction(
                            reactionGraph.getEdgeTarget(outgoingEdge)));
        }
        for (Reaction targetReactionCandidate : targetReactionCandidateSet) {
            SearchUpstreamReactions(
                    depth,
                    reactionGraph,
                    reactionGraph.incomingEdgesOf(targetReactionCandidate.getReactionID()),
                    targetReactionCandidate,
                    targetReactionCandidates,
                    reactomeMechismoDataMap);
        }
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
                                     String rxnFilter) throws Exception {
        long startMethodTime = System.currentTimeMillis();

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

        FisherExact fisherExact = new FisherExact(reactomeMechismoDataMap.getNumPatients());

        Random prng = new Random(88L);
        double maxMemUsed = CalculateJavaMemFootprintGiB();
        for (int i = 0; i < numPermutations; i++) {
            long startLoopTime = System.currentTimeMillis();

            DirectedGraph<Long, DefaultEdge> rewiredReactionGraph =
                    RewireReactionGraph(reactionGraph, prng);

            ConnectivityInspector<Long, DefaultEdge> connectivityInspector =
                    new ConnectivityInspector<>(rewiredReactionGraph);

            System.out.println(String.format("Created reaction graph permutation with %d vertices, %d edges, %d components",
                    rewiredReactionGraph.vertexSet().size(),
                    rewiredReactionGraph.edgeSet().size(),
                    connectivityInspector.connectedSets().size()));

            CooccurrenceResult rewiredNetworkResult = CalculateCooccurrencePValues(
                    fisherExact,
                    rewiredReactionGraph,
                    reactomeMechismoDataMap,
                    depth,
                    ignoreDependentUpstreamReactions,
                    ignoreIndependentUpstreamReactions,
                    excludeMultipleImmediateUpstreamReactions);

            rewiredNetworkResult.CalculateBHAdjustedPValues();

            rewiredNetworkResult.writeToFile(outputDir, i + 1 + "");

            double curMemUsed = CalculateJavaMemFootprintGiB();
            maxMemUsed = curMemUsed > maxMemUsed ?
                    curMemUsed :
                    maxMemUsed;

            rewiredNetworkResult.MagicallyShrinkMemoryFootprint();

            rewiredNetworkResults.add(rewiredNetworkResult);

            long endLoopTime = System.currentTimeMillis();
            System.out.println(String.format("Completed permutation iteration %d in %.2f minutes\n" +
                            "Total running time: %.2f minutes\n" +
                            "Max memory footprint: %.2f GiB",
                    i + 1,//start with iteration '1'
                    (endLoopTime - startLoopTime) / 60000.0,
                    (endLoopTime - startMethodTime) / 60000.0,
                    maxMemUsed));
        }

        CooccurrenceResult realResult = CalculateCooccurrencePValues(
                fisherExact,
                reactionGraph,
                reactomeMechismoDataMap,
                depth,
                ignoreDependentUpstreamReactions,
                ignoreIndependentUpstreamReactions,
                excludeMultipleImmediateUpstreamReactions);

        realResult.CalculateBHAdjustedPValues();
        realResult.CalculateEmpiricalPValues(rewiredNetworkResults);

        WriteRxn2SamplesToFile(outputDir, "", reactomeMechismoDataMap);

        realResult.writeToFile(outputDir, outputFilePrefix);
        realResult.writePatientGroupingsToFile(outputDir, outputFilePrefix);
        realResult.writePatientCluster0UnionDistancesToFile(outputDir,
                outputFilePrefix,
                reactomeMechismoDataMap);

        double curMemUsed = CalculateJavaMemFootprintGiB();
        maxMemUsed = curMemUsed > maxMemUsed ?
                curMemUsed :
                maxMemUsed;

        long endMethodTime = System.currentTimeMillis();

        System.out.println(String.format("Completed after %.2f minutes\n" +
                        "Max memory footprint: %.2f GiB",
                (endMethodTime - startMethodTime) / 60000.0,
                maxMemUsed));
    }

    private double CalculateJavaMemFootprint() {
        return Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
    }

    private double CalculateJavaMemFootprintGiB() {
        return CalculateJavaMemFootprint() / Math.pow(2.0, 30);
    }

    private CooccurrenceResult CalculateCooccurrencePValues(
            FisherExact fisherExact,
            DirectedGraph<Long, DefaultEdge> reactionGraph,
            ReactomeMechismoDataMap reactomeMechismoDataMap,
            int depth,
            boolean ignoreDependentUpstreamReactions,
            boolean ignoreIndependentUpstreamReactions,
            boolean excludeMultipleImmediateUpstreamReactions) throws MathException {

        Set<TargetReactionCandidate> targetReactionCandidates = new HashSet<>();
        SearchRxnNetworkForTargetReactionCandidates(
                targetReactionCandidates,
                reactionGraph,
                reactomeMechismoDataMap,
                depth);

        System.out.println(String.format("Found %d supported Dn/Up reactions",
                targetReactionCandidates.size()));

        return CalculateCooccurrence(
                reactionGraph,
                fisherExact,
                reactomeMechismoDataMap,
                targetReactionCandidates,
                ignoreDependentUpstreamReactions,
                ignoreIndependentUpstreamReactions,
                excludeMultipleImmediateUpstreamReactions);
    }

    private CooccurrenceResult CalculateCooccurrence(
            DirectedGraph<Long, DefaultEdge> reactionGraph,
            FisherExact fisherExact,
            ReactomeMechismoDataMap reactomeMechismoDataMap,
            Set<TargetReactionCandidate> targetReactionCandidates,
            boolean ignoreDependentUpstreamReactions,
            boolean ignoreIndependentUpstreamReactions,
            boolean excludeMultipleImmediateUpstreamReactions) throws MathException {

        List<Reaction> allTargetRxns = new ArrayList<>();
        List<Set<Reaction>> allCooccurringUpstreamRxns = new ArrayList<>();
        List<Set<FI>> allCooccurringUpstreamReactionFIs = new ArrayList<>();
        List<Set<Patient>> allSamplesWTargetRxnMutations = new ArrayList<>();
        List<Integer> allSamplesW0MutatedUpstreamRxns = new ArrayList<>();
        List<Set<Patient>> allSamplesW1MutatedUpstreamRxn = new ArrayList<>();
        List<Set<Patient>> allSamplesW2MutatedUpstreamRxns = new ArrayList<>();
        List<Set<Patient>> allSamplesW3plusMutatedUpstreamRxns = new ArrayList<>();
        List<Set<Mutation>> allSuperIndirectMutatations = new ArrayList<>();
        List<Set<Mutation>> allIndirectMutations = new ArrayList<>();
        List<Set<Mutation>> allSuperDirectMutations = new ArrayList<>();
        List<Set<Mutation>> allDirectMutations = new ArrayList<>();
        List<Double> allPValues = new ArrayList<>();

        int iterationCounter = 0;
        for (TargetReactionCandidate targetReactionCandidate : targetReactionCandidates) {
            Reaction targetRxn = targetReactionCandidate.getTargetReaction();
            List<Reaction[]> targetUpstreamReactionPairs =
                    new ArrayList<>(GenerateReactionSetPairs(targetReactionCandidate.getSupportedUpstreamRxns()));
            Set<Reaction> targetCooccurringUpstreamRxns = new HashSet<>();
            Set<FI> targetCooccurringUpstreamReactionFIs = new HashSet<>();
            Set<Patient> targetSamplesWTargetRxnMutations = new HashSet<>();
            Set<Patient> targetSamplesW0MutatedUpstreamRxns = new HashSet<>();
            Set<Patient> targetSamplesW1MutatedUpstreamRxn = new HashSet<>();
            Set<Patient> targetSamplesW2MutatedUpstreamRxns = new HashSet<>();
            Set<Patient> targetSamplesW3plusMutatedUpstreamRxns = new HashSet<>();
            Set<Mutation> targetSuperIndirectMutations = new HashSet<>();
            Set<Mutation> targetIndirectMutations = new HashSet<>();
            Set<Mutation> targetSuperDirectMutations = new HashSet<>();
            Set<Mutation> targetDirectMutations = new HashSet<>();
            List<Double> targetUpstreamReactionPValues = new ArrayList<>();

            if (!targetUpstreamReactionPairs.isEmpty()) {
                boolean targetReactionDetected = false;
                for (Reaction[] upstreamReactionPair : targetUpstreamReactionPairs) {
                    Reaction upstreamReaction1 = upstreamReactionPair[0];
                    Reaction upstreamReaction2 = upstreamReactionPair[1];

                    if (ignoreDependentUpstreamReactions &&
                            (!reactionGraph.getAllEdges(upstreamReaction1.getReactionID(),
                                    upstreamReaction2.getReactionID()).isEmpty() ||
                                    !reactionGraph.getAllEdges(upstreamReaction2.getReactionID(),
                                            upstreamReaction1.getReactionID()).isEmpty())) {
                        continue;
                    } else if (ignoreIndependentUpstreamReactions &&
                            (reactionGraph.getAllEdges(upstreamReaction1.getReactionID(),
                                    upstreamReaction2.getReactionID()).isEmpty() &&
                                    reactionGraph.getAllEdges(upstreamReaction2.getReactionID(),
                                            upstreamReaction1.getReactionID()).isEmpty())) {
                        continue;
                    } else if (excludeMultipleImmediateUpstreamReactions &&
                            (!reactionGraph.getAllEdges(upstreamReaction1.getReactionID(),
                                    targetRxn.getReactionID()).isEmpty() &&
                                    !reactionGraph.getAllEdges(upstreamReaction2.getReactionID(),
                                            targetRxn.getReactionID()).isEmpty())) {
                        continue;
                    } else {
                        targetReactionDetected = true;

                        Set<Patient> A = new HashSet<>(reactomeMechismoDataMap.getPatients());
                        Set<Patient> B = new HashSet<>(reactomeMechismoDataMap.getPatients());
                        Set<Patient> C = new HashSet<>(reactomeMechismoDataMap.getPatients());
                        Set<Patient> D = new HashSet<>(reactomeMechismoDataMap.getPatients());
                        Set<Patient> upstreamRxn1Samples = new HashSet<>(reactomeMechismoDataMap.getPatients(upstreamReaction1));
                        Set<Patient> upstreamRxn2Samples = new HashSet<>(reactomeMechismoDataMap.getPatients(upstreamReaction2));

                        D.retainAll(upstreamRxn1Samples);
                        D.retainAll(upstreamRxn2Samples);

                        C.removeAll(D);
                        C.retainAll(upstreamRxn2Samples);

                        B.removeAll(D);
                        B.retainAll(upstreamRxn1Samples);

                        A.removeAll(B);
                        A.removeAll(C);
                        A.removeAll(D);

                        if ((A.size() + B.size() + C.size() + D.size()) != reactomeMechismoDataMap.getPatients().size()) {
                            throw new IllegalStateException(String.format(
                                    "A(%d) + B(%d) + C(%d) + D(%d) sum to %d but should sum to %d",
                                    A.size(),
                                    B.size(),
                                    C.size(),
                                    D.size(),
                                    (A.size() + B.size() + C.size() + D.size()),
                                    reactomeMechismoDataMap.getPatients().size()));
                        }

                        Set<Mutation> targetRxnmuts = new HashSet<>();
                        Set<Mutation> upstreamRxn1Mutations = new HashSet<>();
                        Set<Mutation> upstreamRxn2Mutations = new HashSet<>();

                        //Search D samples for indirect, super-indirect, and super-direct mutations
                        if (D.size() > 0) {
                            targetCooccurringUpstreamRxns.add(upstreamReaction1);
                            targetCooccurringUpstreamRxns.add(upstreamReaction2);

                            targetCooccurringUpstreamReactionFIs.addAll(
                                    reactomeMechismoDataMap.getReactionPatientsFIs(upstreamReaction1, D));
                            targetCooccurringUpstreamReactionFIs.addAll(
                                    reactomeMechismoDataMap.getReactionPatientsFIs(upstreamReaction2, D));

                            upstreamRxn1Mutations.addAll(
                                    reactomeMechismoDataMap.getReactionPatientsMutations(
                                            upstreamReaction1,
                                            D));
                            upstreamRxn2Mutations.addAll(
                                    reactomeMechismoDataMap.getReactionPatientsMutations(
                                            upstreamReaction2,
                                            D));
                            targetRxnmuts.addAll(
                                    reactomeMechismoDataMap.getReactionPatientsMutations(
                                            targetRxn,
                                            D));
                        }

                        targetSamplesWTargetRxnMutations.addAll(
                                reactomeMechismoDataMap.getPatients(targetRxn));

                        //upstream union
                        Set<Mutation> patientRxnMutUnion = new HashSet<>(upstreamRxn1Mutations);
                        patientRxnMutUnion.addAll(upstreamRxn2Mutations);

                        //upstream intersection
                        Set<Mutation> sampleRxnMutIntersection = new HashSet<>(upstreamRxn1Mutations);
                        sampleRxnMutIntersection.retainAll(upstreamRxn2Mutations);

                        //upstream union & target intersection
                        Set<Mutation> sampleRxnMutUnionCpy = new HashSet<>(patientRxnMutUnion);
                        sampleRxnMutUnionCpy.retainAll(targetRxnmuts);
                        targetSuperDirectMutations.addAll(sampleRxnMutUnionCpy);

                        //target - upstream union
                        Set<Mutation> targetRxnmutsCpy = new HashSet<>(targetRxnmuts);
                        targetRxnmutsCpy.removeAll(patientRxnMutUnion);
                        targetDirectMutations.addAll(targetRxnmutsCpy);

                        //upstream intersection - target
                        Set<Mutation> sampleRxnMutIntersectionCpy = new HashSet<>(sampleRxnMutIntersection);
                        sampleRxnMutIntersectionCpy.removeAll(targetRxnmuts);
                        targetSuperIndirectMutations.addAll(sampleRxnMutIntersectionCpy);

                        //upstream union - upstream intersection - target
                        sampleRxnMutUnionCpy = new HashSet<>(patientRxnMutUnion);
                        sampleRxnMutUnionCpy.removeAll(sampleRxnMutIntersection);
                        sampleRxnMutUnionCpy.removeAll(targetRxnmuts);
                        targetIndirectMutations.addAll(sampleRxnMutUnionCpy);

                        //samplesW2 D intersection
                        Set<Patient> samplesW2Dintersection = new HashSet<>(D);
                        samplesW2Dintersection.retainAll(targetSamplesW2MutatedUpstreamRxns);
                        targetSamplesW3plusMutatedUpstreamRxns.addAll(samplesW2Dintersection);

                        targetSamplesW0MutatedUpstreamRxns.addAll(A);
                        targetSamplesW1MutatedUpstreamRxn.addAll(B);
                        targetSamplesW1MutatedUpstreamRxn.addAll(C);
                        targetSamplesW2MutatedUpstreamRxns.addAll(D);//remove samplesW3 later

                        //FisherExact uses numerical approximation and is sometimes > 1.0000000000
                        targetUpstreamReactionPValues.add(
                                MathUtilities.boundDouble01(
                                        fisherExact.getRightTailedP(
                                                A.size(),
                                                B.size(),
                                                C.size(),
                                                D.size())));
                    }
                }

                if (targetReactionDetected) {

                    if (allTargetRxns.size() != allPValues.size() ||
                            allTargetRxns.size() != allCooccurringUpstreamRxns.size() ||
                            allTargetRxns.size() != allCooccurringUpstreamReactionFIs.size()) {
                        throw new IllegalStateException(String.format("These should all have size == %d",
                                targetUpstreamReactionPValues.size()));
                    }

                    //we lose track of sample support in upstream reaction pair context
                    targetSamplesW0MutatedUpstreamRxns.removeAll(targetSamplesW1MutatedUpstreamRxn);
                    targetSamplesW0MutatedUpstreamRxns.removeAll(targetSamplesW2MutatedUpstreamRxns);
                    targetSamplesW0MutatedUpstreamRxns.removeAll(targetSamplesW3plusMutatedUpstreamRxns);
                    targetSamplesW0MutatedUpstreamRxns.removeAll(targetSamplesWTargetRxnMutations);

                    targetSamplesW1MutatedUpstreamRxn.removeAll(targetSamplesW2MutatedUpstreamRxns);
                    targetSamplesW1MutatedUpstreamRxn.removeAll(targetSamplesW3plusMutatedUpstreamRxns);
                    targetSamplesW1MutatedUpstreamRxn.removeAll(targetSamplesWTargetRxnMutations);

                    targetSamplesW2MutatedUpstreamRxns.removeAll(targetSamplesW3plusMutatedUpstreamRxns);
                    targetSamplesW2MutatedUpstreamRxns.removeAll(targetSamplesWTargetRxnMutations);

                    targetSamplesW3plusMutatedUpstreamRxns.removeAll(targetSamplesWTargetRxnMutations);

                    //we lose track of super relationships in upstream reaction pair context
                    targetIndirectMutations.removeAll(targetSuperIndirectMutations);
                    targetIndirectMutations.removeAll(targetSuperDirectMutations);

                    targetDirectMutations.removeAll(targetSuperDirectMutations);

                    targetSuperIndirectMutations.removeAll(targetSuperDirectMutations);

                    allTargetRxns.add(targetRxn);
                    allCooccurringUpstreamRxns.add(targetCooccurringUpstreamRxns);
                    allCooccurringUpstreamReactionFIs.add(targetCooccurringUpstreamReactionFIs);
                    allSamplesWTargetRxnMutations.add(targetSamplesWTargetRxnMutations);
                    allSamplesW0MutatedUpstreamRxns.add(targetSamplesW0MutatedUpstreamRxns.size());
                    allSamplesW1MutatedUpstreamRxn.add(targetSamplesW1MutatedUpstreamRxn);
                    allSamplesW2MutatedUpstreamRxns.add(targetSamplesW2MutatedUpstreamRxns);
                    allSamplesW3plusMutatedUpstreamRxns.add(targetSamplesW3plusMutatedUpstreamRxns);
                    allSuperIndirectMutatations.add(targetSuperIndirectMutations);
                    allIndirectMutations.add(targetIndirectMutations);
                    allSuperDirectMutations.add(targetSuperDirectMutations);
                    allDirectMutations.add(targetDirectMutations);
                    allPValues.add(
                            MathUtilities.boundDouble01(
                                    MathUtilities.combinePValuesWithFisherMethod(
                                            targetUpstreamReactionPValues)));
                }
            }

            iterationCounter++;
            if (iterationCounter % 500 == 0) {
                System.out.println(String.format("Processed %d of %d Target Reactions...",
                        iterationCounter, targetReactionCandidates.size()));
                System.out.flush();
            }
        }
        return new CooccurrenceResult(
                allTargetRxns,
                allCooccurringUpstreamRxns,
                allCooccurringUpstreamReactionFIs,
                allSamplesWTargetRxnMutations,
                allSamplesW0MutatedUpstreamRxns,
                allSamplesW1MutatedUpstreamRxn,
                allSamplesW2MutatedUpstreamRxns,
                allSamplesW3plusMutatedUpstreamRxns,
                allSuperIndirectMutatations,
                allIndirectMutations,
                allSuperDirectMutations,
                allDirectMutations,
                allPValues);
    }

    private void FindIntersectionUnionClusters(Set<Set<String>> fiIntersectingSetUnionClusters,
                                               Set<Set<String>> allFIs) {
        Iterator<Set<String>> allFIsItr = allFIs.iterator();
        while (allFIsItr.hasNext()) {
            Set<String> fiSet = allFIsItr.next();
            Set<String> union = new HashSet<>(fiSet);
            Iterator<Set<String>> fiClusterItr = fiIntersectingSetUnionClusters.iterator();
            while (fiClusterItr.hasNext()) {
                Set<String> fiCluster = fiClusterItr.next();
                Set<String> intersection = new HashSet<>(fiSet);
                intersection.retainAll(fiCluster);
                if (!intersection.isEmpty()) {
                    union.addAll(fiCluster);
                    fiClusterItr.remove(); //remove out-dated cluster
                }
            }
            fiIntersectingSetUnionClusters.add(union);
        }
    }

    private void CalculateClusterStats(Set<Set<String>> fiIntersectingSetUnionClusters) {
        Iterator<Set<String>> fiClusterItr = fiIntersectingSetUnionClusters.iterator();
        Double min = (double) fiClusterItr.next().size();
        Double max = min;
        Double sum = max;
        while (fiClusterItr.hasNext()) {
            Double sz = (double) fiClusterItr.next().size();
            min = sz < min ? sz : min;
            max = sz > max ? sz : max;
            sum += sz;
        }
        Double mean = sum / (double) fiIntersectingSetUnionClusters.size();

        System.out.println(String.format("Clusters min = %f, max = %f, mean = %f, sum = %f",
                min,
                max,
                mean,
                sum));
    }

    private void WriteClustersToFile(String outputDir,
                                     String outputFilePrefix,
                                     Set<Set<String>> fiIntersectingSetUnionClusters) {
        //perform pathway enrichment analysis among FI clusters
        FileUtility fileUtility0 = new FileUtility();
        String outFilePath0 = outputDir + outputFilePrefix + "fiClusters.csv";
        try {
            fileUtility0.setOutput(outFilePath0);
            Iterator<Set<String>> fiClusterItr0 = fiIntersectingSetUnionClusters.iterator();
            while (fiClusterItr0.hasNext()) {
                Set<String> fis0 = fiClusterItr0.next();
                fileUtility0.printLine(
                        fis0.size() > 1
                                ? org.gk.util.StringUtils.join(" ",
                                new ArrayList(fis0)).replace("\t", " ")
                                : Arrays.asList(fis0).get(0).toString()
                                .replace("[", "")
                                .replace("]", "")
                                .replace("\t", " "));
            }
            fileUtility0.close();
        } catch (IOException ioe) {
            logger.error(String.format("Couldn't use %s, %s: %s",
                    outFilePath0.toString(),
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    private void WriteReaction2FIsToFile(String outputDir,
                                         String outputFilePrefix,
                                         Map<Long, Set<String>> reaction2FiSet) {
        //write reaction2FIs to file
        FileUtility fileUtility1 = new FileUtility();
        String outFilePath1 = outputDir + outputFilePrefix + "reactions2FIs.csv";
        try {
            fileUtility1.setOutput(outFilePath1);
            fileUtility1.printLine("RxnId,Num FIs,Mapped FIs");
            Iterator<Map.Entry<Long, Set<String>>> reaction2FiItr = reaction2FiSet.entrySet().iterator();
            while (reaction2FiItr.hasNext()) {
                Map.Entry<Long, Set<String>> pair = reaction2FiItr.next();
                Long rxnId = pair.getKey();
                Set<String> mappedFis = pair.getValue();
                fileUtility1.printLine(String.format("%s,%d,%s",
                        rxnId.toString(),
                        mappedFis.size(),
                        mappedFis.size() > 1
                                ? org.gk.util.StringUtils.join(" ",
                                new ArrayList(mappedFis)).replace("\t", " ")
                                : Arrays.asList(mappedFis).get(0).toString()
                                .replace("[", "")
                                .replace("]", "")
                                .replace("\t", " ")));
            }
            fileUtility1.close();
        } catch (IOException ioe) {
            logger.error(String.format("Couldn't use %s, %s: %s",
                    outFilePath1.toString(),
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    private void WriteSamples2FIsToFile(String outputDir,
                                        String outputFilePrefix,
                                        Map<String, Set<String>> samples2FIs) {
        //write samples2FIs to file
        FileUtility fileUtility2 = new FileUtility();
        String outFilePath2 = outputDir + outputFilePrefix + "samples2FIs.csv";
        try {
            fileUtility2.setOutput(outFilePath2);
            fileUtility2.printLine("Sample Barcode,Num FIs,Mapped FIs");
            Iterator<Map.Entry<String, Set<String>>> sample2FiItr = samples2FIs.entrySet().iterator();
            while (sample2FiItr.hasNext()) {
                Map.Entry<String, Set<String>> pair = sample2FiItr.next();
                String sampleBarcode = pair.getKey();
                Set<String> mappedFis = pair.getValue();
                fileUtility2.printLine(String.format("%s,%d,%s",
                        sampleBarcode,
                        mappedFis.size(),
                        mappedFis.size() > 1
                                ? org.gk.util.StringUtils.join(" ",
                                new ArrayList(mappedFis)).replace("\t", " ")
                                : Arrays.asList(mappedFis).get(0).toString()
                                .replace("[", "")
                                .replace("]", "")
                                .replace("\t", " ")));
            }
            fileUtility2.close();
        } catch (IOException ioe) {
            logger.error(String.format("Couldn't use %s, %s: %s",
                    outFilePath2.toString(),
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    private void WriteFIs2SamplesToFile(String outputDir,
                                        String outputFilePrefix,
                                        Map<String, Set<String>> fis2Samples,
                                        Map<String, Set<String>> samples2FIs) {
        //write fis2Samples to file
        FileUtility fileUtility3 = new FileUtility();
        String outFilePath3 = outputDir + outputFilePrefix + "fis2Samples.csv";
        try {
            fileUtility3.setOutput(outFilePath3);
            fileUtility3.printLine("FI,Num Samples,FI Frequency,Mapped Samples");
            Iterator<Map.Entry<String, Set<String>>> fi2SampleItr = fis2Samples.entrySet().iterator();
            while (fi2SampleItr.hasNext()) {
                Map.Entry<String, Set<String>> pair = fi2SampleItr.next();
                String fi = pair.getKey();
                Set<String> mappedSamples = pair.getValue();
                fileUtility3.printLine(String.format("%s,%d,%f,%s",
                        fi,
                        mappedSamples.size(),
                        new Double(mappedSamples.size()) /
                                new Double(samples2FIs.keySet().size()),
                        mappedSamples.size() > 1
                                ? org.gk.util.StringUtils.join(" ",
                                new ArrayList(mappedSamples)).replace(
                                "\t", " ")
                                : Arrays.asList(mappedSamples).get(0).toString()
                                .replace("[", "")
                                .replace("]", "")
                                .replace("\t", " ")));
            }
            fileUtility3.close();
        } catch (IOException ioe) {
            logger.error(String.format("Couldn't use %s, %s: %s",
                    outFilePath3.toString(),
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    private void WriteTargetReactionSummariesToFile(String outputDir,
                                                    String outputFilePrefix,
                                                    Set<TargetReactionCandidate> targetReactionSummaries) {
        //write targetReactionSummaries to file
        FileUtility fileUtility = new FileUtility();
        String outFilePath = outputDir + outputFilePrefix + "dnUpReactions.csv";
        try {
            fileUtility.setOutput(outFilePath);
            fileUtility.printLine(TargetReactionCandidate.getHeaderLine());
            Iterator<TargetReactionCandidate> trsItr = targetReactionSummaries.iterator();
            while (trsItr.hasNext()) {
                TargetReactionCandidate targetReactionCandidate = trsItr.next();
                fileUtility.printLine(targetReactionCandidate.toString());
            }
            fileUtility.close();
        } catch (IOException ioe) {
            logger.error(String.format("Couldn't use %s, %s: %s",
                    outFilePath.toString(),
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }

        System.out.println(String.format("Output written to %s",
                outFilePath));
    }

    private void WriteRxn2SamplesToFile(String outputDir,
                                        String outputFilePrefix,
                                        ReactomeMechismoDataMap reactomeMechismoDataMap) {
        //write rxn2Samples to file
        FileUtility fileUtility = new FileUtility();
        String outFilePath = outputDir + outputFilePrefix + "rxn2Samples.csv";
        try {
            fileUtility.setOutput(outFilePath);
            fileUtility.printLine("Rxn,Num Samples,Rxn Frequency,Mapped Samples");
            for (Reaction reaction : reactomeMechismoDataMap.getSupportedReactions()) {
                Set<Patient> mappedSamples = reactomeMechismoDataMap.getPatients(reaction);
                fileUtility.printLine(String.format("%s,%d,%f,%s",
                        reaction.toString(),
                        mappedSamples.size(),
                        (double) mappedSamples.size() /
                                (double) reactomeMechismoDataMap.getSupportedReactions().size(),
                        mappedSamples.size() > 1
                                ? org.gk.util.StringUtils.join(" ",
                                new ArrayList<>(mappedSamples)).replace("\t", " ")
                                : Collections.singletonList(mappedSamples).get(0).toString()
                                .replace("[", "")
                                .replace("]", "")
                                .replace("\t", " ")));
            }
            fileUtility.close();
        } catch (IOException ioe) {
            logger.error(String.format("Couldn't use %s, %s: %s",
                    outFilePath,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }
}
