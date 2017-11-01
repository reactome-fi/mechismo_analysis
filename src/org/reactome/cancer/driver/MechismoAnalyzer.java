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
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.junit.Test;
import org.reactome.annotate.AnnotationHelper;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.cancer.MechismoOutputLoader;
import org.reactome.cancer.ReactomeReactionGraphLoader;
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
     * @param mechismoOutputFilePath
     * @param mechismoInteractionScorePattern
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
     * @param mechismoOutputFilePath
     * @param mechismoInteractionScorePattern
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

    private Set<Long[]> GenerateReactionSetPairs(Set<Long> reactionSet) {
        Set<Long[]> reactionSetPairs = new HashSet<>();
        List<Long> reactionList = new ArrayList<>(reactionSet);
        for (int i = 0; i < reactionList.size() - 1; i++) { //first to second-last
            for (int j = i + 1; j < reactionList.size(); j++) { //second to last
                if (j == i) {
                    continue;
                }
                Long[] pair = {new Long(reactionList.get(i)), new Long(reactionList.get(j))};
                reactionSetPairs.add(pair);
            }
        }
        return reactionSetPairs;
    }

    private DirectedGraph<Long, DefaultEdge> RewireReactionGraph(
            DefaultDirectedGraph<Long, DefaultEdge> reactionGraph) {
        return RewireReactionGraph(reactionGraph, 88L);
    }


    private DirectedGraph<Long, DefaultEdge> RewireReactionGraph(
            DefaultDirectedGraph<Long, DefaultEdge> reactionGraph,
            long randomSeed) {
        DefaultDirectedGraph<Long, DefaultEdge> rewiredReactionGraph = new DefaultDirectedGraph<>(DefaultEdge.class);
        Graphs.addGraph(rewiredReactionGraph, reactionGraph);

        ConnectivityInspector connectivityInspector = new ConnectivityInspector(rewiredReactionGraph);
        List<Set<Long>> connectedSets = connectivityInspector.connectedSets();

        int initialComponentCount = connectedSets.size();
        int initialEdgeCount = rewiredReactionGraph.edgeSet().size();
        int initialVtxCount = rewiredReactionGraph.vertexSet().size();

        Set<Long> largestComponentVertexSet = new HashSet<>();
        for (Set<Long> connectedSet : connectedSets) {
            if (connectedSet.size() > largestComponentVertexSet.size()) {
                largestComponentVertexSet = connectedSet;
            }
        }

        if (largestComponentVertexSet.size() >= initialVtxCount) {
            throw new IllegalStateException(String.format(
                    "largest component vtx size (%d) should be no less than initial (%d)",
                    largestComponentVertexSet.size(),
                    initialVtxCount));
        }

        Set<DefaultEdge> largestComponentEdges = new HashSet<>();
        for (Long vertex : largestComponentVertexSet) {
            largestComponentEdges.addAll(rewiredReactionGraph.incomingEdgesOf(vertex));
            largestComponentEdges.addAll(rewiredReactionGraph.outgoingEdgesOf(vertex));
        }

        if (largestComponentEdges.size() >= initialEdgeCount) {
            throw new IllegalStateException(String.format(
                    "largest component edge size (%d) should be no larger than initial (%d)",
                    largestComponentEdges.size(),
                    initialEdgeCount));
        }

        Random random = new Random(randomSeed);

        int E = largestComponentEdges.size();
        int Em = E - 1;
        int rand1, rand2;
        DefaultEdge randEdge1, randEdge2;
        ArrayList<DefaultEdge> largestComponentEdgesAry = new ArrayList<>(largestComponentEdges);
        Long sourceVtx1, sourceVtx2, targetVtx1, targetVtx2;

        //from https://en.wikipedia.org/wiki/Degree-preserving_randomization
        int rewires = (int) ((E / 2.0) * Math.log(1 / Math.pow(10, -7)));
        for (int i = 0; i < rewires; i++) {

            //remove edge 1
            rand1 = random.nextInt(E);
            randEdge1 = largestComponentEdgesAry.get(rand1);
            sourceVtx1 = rewiredReactionGraph.getEdgeSource(randEdge1);
            targetVtx1 = rewiredReactionGraph.getEdgeTarget(randEdge1);
            boolean success = rewiredReactionGraph.removeEdge(randEdge1);
            largestComponentEdgesAry.remove(rand1);

            int curEdgeCount = rewiredReactionGraph.edgeSet().size();
            if (curEdgeCount != (initialEdgeCount - 1)) {
                throw new IllegalStateException(String.format(
                        "current edge count (%d) should not differ from initial - 1 (%d)",
                        curEdgeCount,
                        initialEdgeCount - 1));
            }

            //remove edge2
            rand2 = random.nextInt(Em);
            randEdge2 = largestComponentEdgesAry.get(rand2);
            sourceVtx2 = rewiredReactionGraph.getEdgeSource(randEdge2);
            targetVtx2 = rewiredReactionGraph.getEdgeTarget(randEdge2);

            //otherwise we lose edges & degrees
            while (rewiredReactionGraph.getEdge(sourceVtx1, targetVtx2) != null ||
                    rewiredReactionGraph.getEdge(sourceVtx2, targetVtx1) != null) {
                rand2 = random.nextInt(Em);
                randEdge2 = largestComponentEdgesAry.get(rand2);
                sourceVtx2 = rewiredReactionGraph.getEdgeSource(randEdge2);
                targetVtx2 = rewiredReactionGraph.getEdgeTarget(randEdge2);
            }
            success = rewiredReactionGraph.removeEdge(randEdge2);
            largestComponentEdgesAry.remove(rand2);

            curEdgeCount = rewiredReactionGraph.edgeSet().size();
            if (curEdgeCount != (initialEdgeCount - 2)) {
                throw new IllegalStateException(String.format(
                        "current edge count (%d) should not differ from initial - 2 (%d)",
                        curEdgeCount,
                        initialEdgeCount - 2));
            }

            //rewire edge 1
            rewiredReactionGraph.addEdge(sourceVtx1, targetVtx2);
            largestComponentEdgesAry.add(
                    rewiredReactionGraph.getEdge(sourceVtx1, targetVtx2));

            //rewire edge 2
            rewiredReactionGraph.addEdge(sourceVtx2, targetVtx1);
            largestComponentEdgesAry.add(
                    rewiredReactionGraph.getEdge(sourceVtx2, targetVtx1));

            curEdgeCount = rewiredReactionGraph.edgeSet().size();
            if (curEdgeCount != initialEdgeCount) {
                throw new IllegalStateException(String.format(
                        "current edge count (%d) should not differ from initial (%d)",
                        curEdgeCount,
                        initialEdgeCount));
            }

            //if (i % 10000 == 0) {
            //    System.out.println(String.format("Rewired %d of %d Edges...",
            //            i, rewires));
            //}
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

        return rewiredReactionGraph;
    }

    private void MapReactionFIs(Map<String, Set<Long>> fi2ReactionSet,
                                Map<Long, Set<String>> reaction2FiSet,
                                Set<Long> rxns,
                                Set<String> fis,
                                MySQLAdaptor dba) throws Exception {
        Iterator<Long> rxnItr = rxns.iterator();
        int itCounter = 0;
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        Long rxnDbId;
        System.out.println("Processing Reactions...");
        System.out.flush();
        while (rxnItr.hasNext()) {
            rxnDbId = rxnItr.next();
            Set<String> rxnInteractions = new HashSet<>();
            reactomeAnalyzer.generateFIsForReactionsWithFeatures(dba,
                    new ArrayList<>(Arrays.asList(new Long[]{rxnDbId})),
                    null,
                    rxnInteractions);
            //reaction FIs now in interactions set
            Set<String> intersection = new HashSet<>(rxnInteractions);
            intersection.retainAll(fis);
            if (intersection.size() > 0) {
                //if Mechismo FIs in reaction interactions
                reaction2FiSet.put(rxnDbId, intersection);
                for (String interaction :
                        intersection) {
                    Set<Long> interactionRxns;
                    if (fi2ReactionSet.containsKey(interaction)) {
                        interactionRxns = fi2ReactionSet.get(interaction);
                    } else {
                        interactionRxns = new HashSet<>();
                    }
                    interactionRxns.add(rxnDbId);
                    fi2ReactionSet.put(interaction, interactionRxns);
                }
            }
            itCounter++;
            if (itCounter % 100 == 0) {
                //System.out.println(String.format("Processed %d of %d Reactions...",
                //        itCounter, rxns.size()));
                //System.out.flush();
            }
        }
    }

    private void SearchUpstreamReactions(int depth,
                                         Set<DefaultEdge> incomingEdges,
                                         DirectedGraph<Long, DefaultEdge> reactionGraph,
                                         Map<Long, Set<String>> reaction2FiSet,
                                         Long currentReactionId,
                                         Set<Set<String>> allFIs,
                                         Set<TargetReactionSummary> targetReactionSummaries,
                                         Long targetReactionId,
                                         Map<Long, String> longRxnDbIdToName) {
        Set<DefaultEdge> incomingEdgesCpy = new HashSet<>(incomingEdges);
        Set<Long> upstreamReactions = new HashSet<>();
        List<Set<Long>> upstreamReactionsByDepth = new ArrayList<>();
        for (int i = 1; i <= depth; i++) {
            Iterator<DefaultEdge> incomingEdgeItr = incomingEdgesCpy.iterator();
            while (incomingEdgeItr.hasNext()) {
                DefaultEdge incomingEdge = incomingEdgeItr.next();
                Long upstreamReactionId = reactionGraph.getEdgeSource(incomingEdge);
                if (upstreamReactions.contains(upstreamReactionId)) {
                    upstreamReactions.remove(upstreamReactionId);
                } else {
                    upstreamReactions.add(upstreamReactionId);
                }
            }
            if (i < depth) {
                upstreamReactionsByDepth.add(new HashSet<>(upstreamReactions));
                incomingEdgesCpy = new HashSet<>();
                Iterator<Long> upstreamReactionItr = upstreamReactionsByDepth.get(i - 1).iterator();
                while (upstreamReactionItr.hasNext()) {
                    Long upstreamReactionId = upstreamReactionItr.next();
                    incomingEdgesCpy.addAll(reactionGraph.incomingEdgesOf(upstreamReactionId));
                }
            }
        }
        upstreamReactionsByDepth.clear();

        //find intersection of upstream reactions and supported reactions
        Set<Long> supportedUpstreamReactions = new HashSet<>(upstreamReactions);
        supportedUpstreamReactions.retainAll(reaction2FiSet.keySet());
        /*if (!supportedUpstreamReactions.contains(currentReactionId)) {
            throw new IllegalStateException(
                    String.format("Initial rxn ID '%d' not in downstream upstream rxns",
                            currentReactionId)
            );
        }*/
        Set<String> supportedFIs = new HashSet<>();
        Iterator<Long> supportedUpstreamReactionItr = supportedUpstreamReactions.iterator();
        while (supportedUpstreamReactionItr.hasNext()) {
            Long supportedUpstreamReactionId = supportedUpstreamReactionItr.next();
            supportedFIs.addAll(reaction2FiSet.get(supportedUpstreamReactionId));
        }
        allFIs.add(supportedFIs);

        targetReactionSummaries.add(new TargetReactionSummary(
                targetReactionId,
                supportedUpstreamReactions.size(),
                upstreamReactions.size(),
                new Double(supportedUpstreamReactions.size()) /
                        new Double(upstreamReactions.size()),
                supportedUpstreamReactions,
                upstreamReactions,
                supportedFIs,
                longRxnDbIdToName));
    }

    private void SearchRxnNetwork(Set<TargetReactionSummary> targetReactionSummaries,
                                  Set<Set<String>> allFIs,
                                  Map<Long, Set<String>> reaction2FiSet,
                                  DirectedGraph<Long, DefaultEdge> reactionGraph,
                                  Map<Long, String> longRxnDbIdToName,
                                  int depth) {
        //foreach SupReaction in reaction2FiSet
        Iterator<Map.Entry<Long, Set<String>>> rxn2FiItr = reaction2FiSet.entrySet().iterator();
        while (rxn2FiItr.hasNext()) {
            Map.Entry<Long, Set<String>> pair = rxn2FiItr.next();
            Long currentReactionId = pair.getKey();
            //http://jgrapht.org/javadoc/org/jgrapht/graph/DefaultDirectedGraph.html
            //go one DnReaction downstream from SupReaction and
            Set<DefaultEdge> outgoingEdges = reactionGraph.outgoingEdgesOf(currentReactionId);
            Iterator<DefaultEdge> outgoingEdgeItr = outgoingEdges.iterator();
            while (outgoingEdgeItr.hasNext()) {
                DefaultEdge outgoingEdge = outgoingEdgeItr.next();
                Long targetReactionId = reactionGraph.getEdgeTarget(outgoingEdge);
                //one UpReaction from DnReaction
                SearchUpstreamReactions(depth,
                        reactionGraph.incomingEdgesOf(targetReactionId),
                        reactionGraph,
                        reaction2FiSet,
                        currentReactionId,
                        allFIs,
                        targetReactionSummaries,
                        targetReactionId,
                        longRxnDbIdToName);
            }
        }
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
        Double min = new Double(fiClusterItr.next().size());
        Double max = min;
        Double sum = max;
        while (fiClusterItr.hasNext()) {
            Double sz = new Double(fiClusterItr.next().size());
            min = sz < min ? sz : min;
            max = sz > max ? sz : max;
            sum += sz;
        }
        Double mean = sum / new Double(fiIntersectingSetUnionClusters.size());

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
                                                    Set<TargetReactionSummary> targetReactionSummaries) {
        //write targetReactionSummaries to file
        FileUtility fileUtility = new FileUtility();
        String outFilePath = outputDir + outputFilePrefix + "dnUpReactions.csv";
        try {
            fileUtility.setOutput(outFilePath);
            fileUtility.printLine(TargetReactionSummary.headerLine);
            Iterator<TargetReactionSummary> trsItr = targetReactionSummaries.iterator();
            while (trsItr.hasNext()) {
                TargetReactionSummary targetReactionSummary = trsItr.next();
                fileUtility.printLine(targetReactionSummary.toString());
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

    private void MapReactionsToSamples(Map<Long, Set<String>> reaction2FiSet,
                                       Map<Long, Set<String>> rxn2Samples,
                                       Map<String, Set<String>> fis2Samples) {
        Iterator<Map.Entry<Long, Set<String>>> rxn2FisItr = reaction2FiSet.entrySet().iterator();
        while (rxn2FisItr.hasNext()) {
            Map.Entry<Long, Set<String>> pair = rxn2FisItr.next();
            Long rxn = new Long(pair.getKey());
            Set<String> rxnFis = pair.getValue();
            Iterator<String> rxnFisItr = rxnFis.iterator();
            while (rxnFisItr.hasNext()) {
                String rxnFi = rxnFisItr.next();
                Set<String> rxnSamples;
                if (rxn2Samples.containsKey(rxn)) {
                    rxnSamples = new HashSet<>(rxn2Samples.get(rxn));
                } else {
                    rxnSamples = new HashSet<>();
                }
                rxnSamples.addAll(fis2Samples.get(rxnFi));
                rxn2Samples.put(rxn, rxnSamples);
            }
        }
    }

    private void WriteRxn2SamplesToFile(String outputDir,
                                        String outputFilePrefix,
                                        Map<Long, Set<String>> rxn2Samples,
                                        Map<String, Set<String>> samples2FIs) {
        //write rxn2Samples to file
        FileUtility fileUtility4 = new FileUtility();
        String outFilePath4 = outputDir + outputFilePrefix + "rxn2Samples.csv";
        try {
            fileUtility4.setOutput(outFilePath4);
            fileUtility4.printLine("Rxn,Num Samples,Rxn Frequency,Mapped Samples");
            Iterator<Map.Entry<Long, Set<String>>> rxn2SampleItr = rxn2Samples.entrySet().iterator();
            while (rxn2SampleItr.hasNext()) {
                Map.Entry<Long, Set<String>> pair = rxn2SampleItr.next();
                Long rxn = new Long(pair.getKey());
                Set<String> mappedSamples = pair.getValue();
                if (mappedSamples.size() != rxn2Samples.get(
                        new Long(rxn.toString())).size() ||
                        mappedSamples.size() < 1) {
                    throw new IllegalStateException("These should always match");
                }
                fileUtility4.printLine(String.format("%s,%d,%f,%s",
                        rxn.toString(),
                        mappedSamples.size(),
                        new Double(mappedSamples.size()) /
                                new Double(samples2FIs.keySet().size()),
                        mappedSamples.size() > 1
                                ? org.gk.util.StringUtils.join(" ",
                                new ArrayList(mappedSamples)).replace("\t", " ")
                                : Arrays.asList(mappedSamples).get(0).toString()
                                .replace("[", "")
                                .replace("]", "")
                                .replace("\t", " ")));
            }
            fileUtility4.close();
        } catch (IOException ioe) {
            logger.error(String.format("Couldn't use %s, %s: %s",
                    outFilePath4.toString(),
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    public void mapReactomeReactions(CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer,
                                     String mechismoOutputFilePath,
                                     String reactomeReactionNetworkFilePath,
                                     String outputDir,
                                     String outputFilePrefix,
                                     int depth,
                                     int numPermutations,
                                     double mechScoreLowerBoundInclusive,
                                     double eThresh,
                                     boolean ignoreDependentUpstreamReactions,
                                     boolean ignoreIndependentUpstreamReactions,
                                     boolean excludeMultipleImmediateUpstreamReactions) throws Exception {
        long startMethodTime = System.currentTimeMillis();

        //build jgrapht network from reaction file
        ReactomeReactionGraphLoader reactomeReactionGraphLoader =
                new ReactomeReactionGraphLoader(reactomeReactionNetworkFilePath);


        //extract FIs from Mechismo output
        MechismoOutputLoader mechismoOutputLoader =
                new MechismoOutputLoader(mechismoOutputFilePath,
                        mechScoreLowerBoundInclusive,
                        eThresh);

        //Set<FI>
        Set<String> fis = mechismoOutputLoader.ExtractMechismoFIs();

        System.out.println(String.format("Loaded %d Mechismo FIs",
                fis.size()));

        //Map<Sample Barcode,Set<FI>>
        Map<String, Set<String>> samples2FIs = mechismoOutputLoader.ExtractSamples2FIs();

        System.out.println(String.format("%s samples mapped to FIs",
                samples2FIs.keySet().size()));

        //Map<FI,Set<Sample Barcode>>
        Map<String, Set<String>> fis2Samples = mechismoOutputLoader.ExtractFIs2Samples();

        System.out.println(String.format("%s FIs mapped to samples",
                fis2Samples.keySet().size()));

        //find reactions containing Mechismo FIs
        //FI -> Set<Reaction Id>
        Map<String, Set<Long>> fi2ReactionSet = new HashMap<>();
        //Reaction Id -> Set<FIs>
        Map<Long, Set<String>> reaction2FiSet = new HashMap<>();
        MapReactionFIs(fi2ReactionSet,
                reaction2FiSet,
                reactomeReactionGraphLoader.getReactionSet(),
                fis,
                cancerDriverReactomeAnalyzer.getDBA());

        System.out.println(String.format("Mapped %d FIs to %d reactions",
                fi2ReactionSet.keySet().size(),
                reaction2FiSet.keySet().size()));

        Map<Long, String> longRxnDbIdToName =
                cancerDriverReactomeAnalyzer.loadReactionLongDBIDToName();

        Map<Long, Set<String>> rxn2Samples = new HashMap<>();
        MapReactionsToSamples(reaction2FiSet,
                rxn2Samples,
                fis2Samples);

        List<CooccurrenceResult> rewiredNetworkResults = new ArrayList<>();

        DefaultDirectedGraph<Long, DefaultEdge> reactionGraph =
                reactomeReactionGraphLoader.getReactionGraph();

        Map<String, Map<String, Set<String>>> samples2FIs2Muts = mechismoOutputLoader.ExtractSamples2FIs2Muts();
        double maxMemUsed = CalculateJavaMemFootprintGiB();
        for (int i = 0; i < numPermutations; i++) {
            long startLoopTime = System.currentTimeMillis();

            DirectedGraph<Long, DefaultEdge> rewiredReactionGraph =
                    RewireReactionGraph(reactionGraph, i);

            ConnectivityInspector<Long, DefaultEdge> connectivityInspector =
                    new ConnectivityInspector<>(rewiredReactionGraph);

            System.out.println(String.format("Created reaction graph with %d vertices, %d edges, %d components",
                    rewiredReactionGraph.vertexSet().size(),
                    rewiredReactionGraph.edgeSet().size(),
                    connectivityInspector.connectedSets().size()));

            //for(Set<Long> component : connectivityInspector.connectedSets()){
            //    System.out.println(component.size());
            //}

            CooccurrenceResult rewiredNetworkResult = CalculateCooccurrencePValues(
                    longRxnDbIdToName,
                    rewiredReactionGraph,
                    samples2FIs,
                    reaction2FiSet,
                    rxn2Samples,
                    samples2FIs2Muts,
                    depth,
                    false,
                    ignoreDependentUpstreamReactions,
                    ignoreIndependentUpstreamReactions,
                    excludeMultipleImmediateUpstreamReactions);

            rewiredNetworkResult.MagicallyShrinkMemoryFootprint();

            rewiredNetworkResults.add(rewiredNetworkResult);

            double curMemUsed = CalculateJavaMemFootprintGiB();
            maxMemUsed = curMemUsed > maxMemUsed ?
                    curMemUsed :
                    maxMemUsed;

            long endLoopTime = System.currentTimeMillis();
            System.out.println(String.format("Completed permutation iteration %d in %.2f minutes\n" +
                            "Total running time: %.2f minutes\n" +
                            "Max memory footprint: %.2f GiB",
                    i,
                    (endLoopTime - startLoopTime) / 60000.0,
                    (endLoopTime - startMethodTime) / 60000.0,
                    maxMemUsed));
        }

        CooccurrenceResult realResult = CalculateCooccurrencePValues(
                longRxnDbIdToName,
                reactionGraph,
                samples2FIs,
                reaction2FiSet,
                rxn2Samples,
                samples2FIs2Muts,
                depth,
                true,
                ignoreDependentUpstreamReactions,
                ignoreIndependentUpstreamReactions,
                excludeMultipleImmediateUpstreamReactions);

        realResult.CalculateNetworkPermutationFDR(rewiredNetworkResults);

        WriteRxn2SamplesToFile(outputDir,
                outputFilePrefix,
                rxn2Samples,
                samples2FIs);

        WriteCooccurrenceToFile(
                realResult,
                outputDir,
                outputFilePrefix,
                longRxnDbIdToName);

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

    private void WriteCooccurrenceToFile(CooccurrenceResult cooccurrenceResult,
                                         String outputDir,
                                         String outputFilePrefix,
                                         Map<Long, String> longRxnDbIdToName) {

        String outFilePath5 = outputDir + outputFilePrefix + "rxnCooccurrence.csv";
        FileUtility fileUtility = new FileUtility();

        try {
            fileUtility.setOutput(outFilePath5);
            fileUtility.printLine(
                    "Target Reaction," +
                            "Upstream Reactions," +
                            "A," +
                            "B," +
                            "C," +
                            "D," +
                            "Large D," +
                            "Large D muts," +
                            "B samples," +
                            "C samples," +
                            "D samples," +
                            "D sample FIs Rxn 1," +
                            "D sample FIs Rxn 2," +
                            "Co-occurrence P-value," +
                            "Co-occurrence Network Permutation FDR");

            for (int i = 0; i < cooccurrenceResult.targetReactionList.size(); i++) {
                WriteLineToFile(fileUtility,
                        cooccurrenceResult.targetReactionList,
                        i,
                        longRxnDbIdToName,
                        cooccurrenceResult.upstreamReactionList,
                        cooccurrenceResult.ABCDs,
                        cooccurrenceResult.BCDSamples,
                        cooccurrenceResult.DSampleFI1s,
                        cooccurrenceResult.DSampleFI2s,
                        cooccurrenceResult.pValues,
                        cooccurrenceResult.pValue2FDRMap);
            }
            fileUtility.close();
        } catch (IOException ioe) {
            logger.error(String.format("Couldn't use %s, %s: %s",
                    outFilePath5.toString(),
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    private void WriteLineToFile(FileUtility fileUtility5,
                                 List<String> targetReactionList,
                                 int i,
                                 Map<Long, String> longRxnDbIdToName,
                                 List<String> upstreamReactionList,
                                 List<int[]> ABCDs,
                                 List<Set<String>[]> BCDSamples,
                                 List<Set<String>> DSampleFI1s,
                                 List<Set<String>> DSampleFI2s,
                                 List<Double> pValues,
                                 Map<Double, Double> pValue2FDRMap) throws IOException {
        fileUtility5.printLine(String.format(
                "%s:%s,%s,%d,%d,%d,%d,%d,%d,%s,%s,%s,%s,%s,%.20e,%.20e",
                targetReactionList.get(i),
                longRxnDbIdToName.get(
                        Long.parseLong(
                                targetReactionList.get(i))).replace(
                        ",", "~"),
                upstreamReactionList.get(i),
                ABCDs.get(i)[0],
                ABCDs.get(i)[1],
                ABCDs.get(i)[2],
                ABCDs.get(i)[3],
                ABCDs.get(i)[4],
                ABCDs.get(i)[5],
                BCDSamples.get(i)[0].size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(BCDSamples.get(i)[0]))
                        : Arrays.asList(BCDSamples.get(i)[0]).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                BCDSamples.get(i)[1].size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(BCDSamples.get(i)[1]))
                        : Arrays.asList(BCDSamples.get(i)[1]).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                BCDSamples.get(i)[2].size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(BCDSamples.get(i)[2]))
                        : Arrays.asList(BCDSamples.get(i)[2]).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                DSampleFI1s.get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(DSampleFI1s.get(i)))
                        : Arrays.asList(DSampleFI1s.get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                DSampleFI2s.get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(DSampleFI2s.get(i)))
                        : Arrays.asList(DSampleFI2s.get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                pValues.get(i),
                pValue2FDRMap.get(pValues.get(i))));
    }

    private CooccurrenceResult CalculateCooccurrencePValues(
            Map<Long, String> longRxnDbIdToName,
            DirectedGraph<Long, DefaultEdge> reactionGraph,
            Map<String, Set<String>> samples2FIs,
            Map<Long, Set<String>> reaction2FiSet,
            Map<Long, Set<String>> rxn2Samples,
            Map<String, Map<String, Set<String>>> samples2fis2muts,
            int depth,
            boolean rememberEverything,
            boolean ignoreDependentUpstreamReactions,
            boolean ignoreIndependentUpstreamReactions,
            boolean excludeMultipleImmediateUpstreamReactions) {

        Set<TargetReactionSummary> targetReactionSummaries = new HashSet<>();
        Set<Set<String>> allFIs = new HashSet<>();
        SearchRxnNetwork(targetReactionSummaries,
                allFIs,
                reaction2FiSet,
                reactionGraph,
                longRxnDbIdToName,
                depth);

        System.out.println(String.format("Found %d supported Dn/Up reactions",
                targetReactionSummaries.size()));

        System.out.println(String.format("Found %d unique supporting FI sets",
                allFIs.size()));

        /*
        Set<Set<String>> fiIntersectingSetUnionClusters = new HashSet<>();
        FindIntersectionUnionClusters(fiIntersectingSetUnionClusters,
                allFIs);

        System.out.println(String.format("Aggregated %d unique FI sets into %d intersection-unions",
                allFIs.size(),
                fiIntersectingSetUnionClusters.size()));

        CalculateClusterStats(fiIntersectingSetUnionClusters);
        */

        //Map<Reaction,Set<Sample Barcode>>
        //report reaction frequencies

        System.gc();

        FisherExact fisherExact = new FisherExact(samples2FIs.keySet().size());
        List<String> targetReactionList = new ArrayList<>();
        List<String> upstreamReactionList = new ArrayList<>();
        List<Double> pValues = new ArrayList<>();
        List<int[]> abcds = new ArrayList<>();
        List<Set<String>[]> bcdSamples = new ArrayList<>();
        List<Set<String>> dSampleFI1s = new ArrayList<>();
        List<Set<String>> dSampleFI2s = new ArrayList<>();

        CalculateCooccurrence(targetReactionSummaries,
                longRxnDbIdToName,
                reactionGraph,
                rxn2Samples,
                samples2FIs,
                reaction2FiSet,
                fisherExact,
                targetReactionList,
                upstreamReactionList,
                pValues,
                abcds,
                bcdSamples,
                dSampleFI1s,
                dSampleFI2s,
                samples2fis2muts,
                rememberEverything,
                ignoreDependentUpstreamReactions,
                ignoreIndependentUpstreamReactions,
                excludeMultipleImmediateUpstreamReactions);

        //Map<Double, Double> inPValue2FDRMap = CalculateFDRs(pValues);

        return (new CooccurrenceResult(targetReactionList,
                upstreamReactionList,
                pValues,
                abcds,
                bcdSamples,
                dSampleFI1s,
                dSampleFI2s));
    }

    private void CalculateCooccurrence(Set<TargetReactionSummary> targetReactionSummaries,
                                       Map<Long, String> longRxnDbIdToName,
                                       DirectedGraph<Long, DefaultEdge> reactionGraph,
                                       Map<Long, Set<String>> rxn2Samples,
                                       Map<String, Set<String>> samples2FIs,
                                       Map<Long, Set<String>> reaction2FiSet,
                                       FisherExact fisherExact,
                                       List<String> targetReactionList,
                                       List<String> upstreamReactionList,
                                       List<Double> inPValues,
                                       List<int[]> inABCDs,
                                       List<Set<String>[]> inBCDSamples,
                                       List<Set<String>> inDSampleFI1s,
                                       List<Set<String>> inDSampleFI2s,
                                       Map<String, Map<String, Set<String>>> samples2fis2muts,
                                       boolean rememberEverything,
                                       boolean ignoreDependentUpstreamReactions,
                                       boolean ignoreIndependentUpstreamReactions,
                                       boolean excludeMultipleImmediateUpstreamReactions) {
        Iterator<TargetReactionSummary> trsItr = targetReactionSummaries.iterator();
        int trsCounter = 0;
        while (trsItr.hasNext()) {
            TargetReactionSummary trs = trsItr.next();
            Long targetReactionId = trs.rxnId;
            List<Long[]> upstreamReactionPairs =
                    new ArrayList<>(GenerateReactionSetPairs(trs.supRxns));
            List<String> upstreamReactionPairStrings = new ArrayList<>();
            List<Double> cooccurrencePvalues = new ArrayList<>();
            List<int[]> inABCD = new ArrayList<>();
            List<Set<String>[]> inBCDSample = new ArrayList<>();
            List<Set<String>> inDSampleFI1 = new ArrayList<>();
            List<Set<String>> inDSampleFI2 = new ArrayList<>();
            for (Long[] upstreamReactionPair : upstreamReactionPairs) {
                Long upstreamReaction1Id = upstreamReactionPair[0];
                Long upstreamReaction2Id = upstreamReactionPair[1];

                if(upstreamReaction1Id == 114548L ||
                        upstreamReaction2Id == 114548L){
                    int debug = 1;
                }

                if (ignoreDependentUpstreamReactions &&
                        !reactionGraph.getAllEdges(upstreamReaction1Id, upstreamReaction2Id).isEmpty() ||
                        !reactionGraph.getAllEdges(upstreamReaction2Id, upstreamReaction1Id).isEmpty()) {
                    continue;
                }

                if (ignoreIndependentUpstreamReactions &&
                        reactionGraph.getAllEdges(upstreamReaction1Id, upstreamReaction2Id).isEmpty() &&
                        reactionGraph.getAllEdges(upstreamReaction2Id, upstreamReaction1Id).isEmpty()) {
                    continue;
                }

                if (excludeMultipleImmediateUpstreamReactions &&
                        !reactionGraph.getAllEdges(upstreamReaction1Id, targetReactionId).isEmpty() &&
                        !reactionGraph.getAllEdges(upstreamReaction2Id, targetReactionId).isEmpty()) {
                    continue;
                }


                upstreamReactionPairStrings.add(
                        ConvertReactionPairToString(upstreamReaction1Id,
                                upstreamReaction2Id
                                , longRxnDbIdToName));

                //one row
                Set<String> A = new HashSet<>(samples2FIs.keySet());
                Set<String> B = new HashSet<>(samples2FIs.keySet());
                Set<String> C = new HashSet<>(samples2FIs.keySet());
                Set<String> D = new HashSet<>(samples2FIs.keySet());
                Set<String> upstreamRxn1Samples = new HashSet<>(rxn2Samples.get(upstreamReaction1Id));
                Set<String> upstreamRxn2Samples = new HashSet<>(rxn2Samples.get(upstreamReaction2Id));

                D.retainAll(upstreamRxn1Samples);
                D.retainAll(upstreamRxn2Samples);

                C.removeAll(D);
                C.retainAll(upstreamRxn2Samples);

                B.removeAll(D);
                B.retainAll(upstreamRxn1Samples);

                A.removeAll(B);
                A.removeAll(C);
                A.removeAll(D);

                if ((A.size() + B.size() + C.size() + D.size()) != samples2FIs.keySet().size()) {
                    int debug = 1;
                }

                Set<String> largeD = new HashSet<>();
                Set<String> largeDMuts = new HashSet<>();

                //TODO:filter d by mutation
                Iterator<String> dSampleItr = D.iterator();
                while (dSampleItr.hasNext()) {
                    String sample = dSampleItr.next();
                    Set<String> rxn1muts = new HashSet<>();
                    Set<String> sampleFIsSupportingRxn1 = findSampleFIsSupportingRxn(
                            sample,
                            samples2FIs,
                            reaction2FiSet,
                            upstreamReaction1Id);
                    for (String fi : sampleFIsSupportingRxn1) {
                        rxn1muts.addAll(
                                samples2fis2muts.get(sample).get(fi)
                        );
                    }
                    Set<String> rxn2muts = new HashSet<>();
                    Set<String> sampleFIsSupportingRxn2 = findSampleFIsSupportingRxn(
                            sample,
                            samples2FIs,
                            reaction2FiSet,
                            upstreamReaction2Id);
                    for (String fi : sampleFIsSupportingRxn2) {
                        rxn2muts.addAll(
                                samples2fis2muts.get(sample).get(fi)
                        );
                    }

                    Set<String> sampleRxnMutIntersection = new HashSet<>(rxn1muts);
                    sampleRxnMutIntersection.retainAll(rxn2muts);

                    if (!sampleRxnMutIntersection.isEmpty()) {
                        largeD.add(sample);
                        largeDMuts.addAll(sampleRxnMutIntersection);
                        dSampleItr.remove();
                    }
                }

                int[] abcd = {A.size(), B.size(), C.size(), D.size(), largeD.size(), largeDMuts.size()};
                inABCD.add(abcd);

                Set[] bcdSample = new Set[]{B, C, D};
                inBCDSample.add(bcdSample);

                inDSampleFI1.add(findSampleFIsSupportingRxn(D,
                        samples2FIs,
                        reaction2FiSet,
                        upstreamReaction1Id));

                inDSampleFI2.add(findSampleFIsSupportingRxn(D,
                        samples2FIs,
                        reaction2FiSet,
                        upstreamReaction2Id));

                cooccurrencePvalues.add(
                        fisherExact.getTwoTailedP(
                                A.size(),
                                B.size(),
                                C.size(),
                                D.size()));

            }
            if (rememberEverything) {
                fillAccumulators(cooccurrencePvalues,
                        targetReactionList,
                        targetReactionId,
                        upstreamReactionList,
                        upstreamReactionPairStrings,
                        inPValues,
                        inABCDs,
                        inABCD,
                        inBCDSamples,
                        inBCDSample,
                        inDSampleFI1s,
                        inDSampleFI1,
                        inDSampleFI2s,
                        inDSampleFI2);
            } else {
                fillPvalueAccumulator(cooccurrencePvalues,
                        inPValues);
            }
            trsCounter++;
            if (trsCounter % 500 == 0) {
                System.out.println(String.format("Processed %d of %d Target Reactions...",
                        trsCounter, targetReactionSummaries.size()));
                System.out.flush();
            }
        }
    }

    private Set<String> findSampleFIsSupportingRxn(Set<String> samples,
                                                   Map<String, Set<String>> samples2FIs,
                                                   Map<Long, Set<String>> reaction2FiSet,
                                                   Long upstreamReactionId) {
        Set<String> sampleFIsSupportingRxn = new HashSet<>();
        for (String sample : samples) {
            Set<String> sampleFIs = new HashSet<>(samples2FIs.get(sample));
            Set<String> supportingFIs = new HashSet<>(reaction2FiSet.get(upstreamReactionId));
            sampleFIs.retainAll(supportingFIs);
            sampleFIsSupportingRxn.addAll(sampleFIs);
        }
        return sampleFIsSupportingRxn;
    }

    private Set<String> findSampleFIsSupportingRxn(String sample,
                                                   Map<String, Set<String>> samples2FIs,
                                                   Map<Long, Set<String>> reaction2FiSet,
                                                   Long upstreamReactionId) {
        Set<String> sampleSet = new HashSet<>();
        sampleSet.add(sample);
        return findSampleFIsSupportingRxn(sampleSet,
                samples2FIs,
                reaction2FiSet,
                upstreamReactionId);
    }

    /*
    private void CalculateEXCooccurrence(Set<TargetReactionSummary> targetReactionSummaries,
                                         Map<Long, String> longRxnDbIdToName,
                                         DefaultDirectedGraph<Long, DefaultEdge> reactionGraph,
                                         Map<String, Set<String>> samples2FIs,
                                         Map<Long, Set<String>> reaction2FiSet,
                                         FisherExact fisherExact,
                                         Map<String, Set<String>> fis2Samples,
                                         List<String> targetReactionList,
                                         List<String> upstreamReactionList,
                                         List<Double> exPValues,
                                         List<int[]> exABCDs,
                                         List<Set<String>[]> exBCDSamples,
                                         List<Set<String>> exDSampleFI1s,
                                         List<Set<String>> exDSampleFI2s,
                                         Map<String, Map<String, Set<String>>> samples2fis2muts) {
        Iterator<TargetReactionSummary> trsItr = targetReactionSummaries.iterator();
        int trsCounter = 0;
        while (trsItr.hasNext()) {
            TargetReactionSummary trs = trsItr.next();
            Long targetReactionId = trs.rxnId;
            List<Long[]> upstreamReactionPairs =
                    new ArrayList<>(GenerateReactionSetPairs(trs.supRxns));
            List<String> upstreamReactionPairStrings =
                    ConvertReactionPairToString(upstreamReactionPairs, longRxnDbIdToName);
            List<Double> cooccurrencePvalues = new ArrayList<>();
            List<int[]> exABCD = new ArrayList<>();
            List<Set<String>[]> exBCDSample = new ArrayList<>();
            List<Set<String>> exDSampleFI1 = new ArrayList<>();
            List<Set<String>> exDSampleFI2 = new ArrayList<>();
            for (Long[] upstreamReactionPair : upstreamReactionPairs) {
                Long upstreamReaction1Id = upstreamReactionPair[0];
                Long upstreamReaction2Id = upstreamReactionPair[1];



                //ignore interdependent pairs
                if (reactionGraph.getAllEdges(upstreamReaction1Id, upstreamReaction2Id).isEmpty() &&
                        reactionGraph.getAllEdges(upstreamReaction2Id, upstreamReaction1Id).isEmpty() &&
                        !upstreamReaction1Id.equals(upstreamReaction2Id)) {

                    //one row for EX
                    Set<String> dnUp1FisExclusive =
                            new HashSet<>(reaction2FiSet.get(upstreamReaction1Id));
                    dnUp1FisExclusive.removeAll(reaction2FiSet.get(upstreamReaction2Id));

                    Set<String> dnUp2FisExclusive =
                            new HashSet<>(reaction2FiSet.get(upstreamReaction2Id));
                    dnUp2FisExclusive.removeAll(reaction2FiSet.get(upstreamReaction1Id));

                    Set<String> fi1sGenes = new HashSet<>();
                    for (String fi : dnUp1FisExclusive) {
                        String[] genes = fi.split("\t");
                        fi1sGenes.add(genes[0]);
                        fi1sGenes.add(genes[1]);
                    }

                    Set<String> fi2sGenes = new HashSet<>();
                    for (String fi : dnUp2FisExclusive) {
                        String[] genes = fi.split("\t");
                        fi2sGenes.add(genes[0]);
                        fi2sGenes.add(genes[1]);
                    }

                    Set<String> fisGeneIntersection = new HashSet<>(fi1sGenes);
                    fisGeneIntersection.retainAll(fi2sGenes);

                    Set<String> dnUp1SamplesExclusive = new HashSet<>();
                    for (String dnUp1Fi : dnUp1FisExclusive) {
                        String[] genes = dnUp1Fi.split("\t");
                        if (!fisGeneIntersection.contains(genes[0]) &&
                                !fisGeneIntersection.contains(genes[1])) {
                            dnUp1SamplesExclusive.addAll(fis2Samples.get(dnUp1Fi));
                        }
                    }
                    Set<String> dnUp2SamplesExclusive = new HashSet<>();
                    for (String dnUp2Fi : dnUp2FisExclusive) {
                        String[] genes = dnUp2Fi.split("\t");
                        if (!fisGeneIntersection.contains(genes[0]) &&
                                !fisGeneIntersection.contains(genes[1])) {
                            dnUp2SamplesExclusive.addAll(fis2Samples.get(dnUp2Fi));
                        }
                    }

                    Set<String> acE = new HashSet<>(samples2FIs.keySet());
                    acE.removeAll(dnUp1SamplesExclusive);

                    Set<String> bdE = new HashSet<>(dnUp1SamplesExclusive);

                    Set<String> aE = new HashSet<>(acE);
                    aE.removeAll(dnUp2SamplesExclusive);

                    Set<String> cE = new HashSet<>(acE);
                    cE.retainAll(dnUp2SamplesExclusive);

                    Set<String> bE = new HashSet<>(bdE);
                    bE.removeAll(dnUp2SamplesExclusive);

                    Set<String> dE = new HashSet<>(bdE);
                    dE.retainAll(dnUp2SamplesExclusive);

                    //TODO: filter dE by mutation

                    int[] abcdE = {aE.size(), bE.size(), cE.size(), dE.size()};
                    exABCD.add(abcdE);

                    Set[] bcdSampleE = new Set[]{bE, cE, dE};
                    exBCDSample.add(bcdSampleE);

                    exDSampleFI1.add(findSampleFIsSupportingRxn(dE,
                            samples2FIs,
                            reaction2FiSet,
                            upstreamReaction1Id));

                    exDSampleFI2.add(findSampleFIsSupportingRxn(dE,
                            samples2FIs,
                            reaction2FiSet,
                            upstreamReaction2Id));

                    cooccurrencePvalues.add(
                            fisherExact.getTwoTailedP(
                                    aE.size(),
                                    bE.size(),
                                    cE.size(),
                                    dE.size()));
                } else {
                    //System.out.println(String.format("Ignoring interdependent upstream reaction pair %s-%s",
                    //        dnUp1,dnUp2));
                }
            }
            fillAccumulators(cooccurrencePvalues,
                    targetReactionList,
                    targetReactionId,
                    upstreamReactionList,
                    upstreamReactionPairStrings,
                    exPValues,
                    exABCDs,
                    exABCD,
                    exBCDSamples,
                    exBCDSample,
                    exDSampleFI1s,
                    exDSampleFI1,
                    exDSampleFI2s,
                    exDSampleFI2);
            trsCounter++;
            if (trsCounter % 100 == 0) {
                System.out.println(String.format("Processed %d of %d Target Reactions...",
                        trsCounter, targetReactionSummaries.size()));
                System.out.flush();
            }
        }
    }
    */

    private void fillAccumulators(List<Double> cooccurrencePvalues,
                                  List<String> targetReactionList,
                                  Long targetReactionId,
                                  List<String> upstreamReactionList,
                                  List<String> upstreamReactionPairStrings,
                                  List<Double> pValues,
                                  List<int[]> ABCDs,
                                  List<int[]> ABCD,
                                  List<Set<String>[]> BCDSamples,
                                  List<Set<String>[]> BCDSample,
                                  List<Set<String>> DSampleRxn1FIs,
                                  List<Set<String>> DSampleRxn1FI,
                                  List<Set<String>> DSampleRxn2FIs,
                                  List<Set<String>> DSampleRxn2FI) {

        if (cooccurrencePvalues.size() != upstreamReactionPairStrings.size() ||
                cooccurrencePvalues.size() != ABCD.size() ||
                cooccurrencePvalues.size() != BCDSample.size() ||
                cooccurrencePvalues.size() != DSampleRxn1FI.size() ||
                cooccurrencePvalues.size() != DSampleRxn2FI.size()) {
            throw new IllegalStateException(String.format("These should all have size == %d",
                    cooccurrencePvalues.size()));
        }

        for (int i = 0; i < cooccurrencePvalues.size(); i++) {
            targetReactionList.add(targetReactionId.toString());
            upstreamReactionList.add(upstreamReactionPairStrings.get(i));
            pValues.add(cooccurrencePvalues.get(i));
            ABCDs.add(ABCD.get(i));
            BCDSamples.add(BCDSample.get(i));
            DSampleRxn1FIs.add(DSampleRxn1FI.get(i));
            DSampleRxn2FIs.add(DSampleRxn2FI.get(i));
        }
    }

    private void fillPvalueAccumulator(List<Double> cooccurrencePvalues, List<Double> pValues) {
        for (Double pValue : cooccurrencePvalues) {
            pValues.add(pValue);
        }
    }

    private Map<Double, Double> CalculateFDRs(List<Double> pValues) {
        //The FDR calculation has the side effect of sorting the passed list...
        List<Double> pValuesSorted = new ArrayList<>(pValues);

        List<Double> FDRs =
                MathUtilities.calculateFDRWithBenjaminiHochberg(pValuesSorted);

        Map<Double, Double> pValue2FDRMap = new HashMap<>();
        for (int i = 0; i < pValues.size(); i++) {
            pValue2FDRMap.put(pValuesSorted.get(i), FDRs.get(i));
        }
        return pValue2FDRMap;
    }

    private Double sumSet(Set<Double> doubleSet) {
        Double sum = 0.0;
        for (Double dbl : doubleSet) {
            sum += dbl;
        }
        return sum;
    }

    private String ConvertReactionPairToString(Long reaction1Id,
                                                Long reaction2Id,
                                               Map<Long, String> longRxnDbIdToName) {
        return String.format("%d:%s|%d:%s",
                    reaction1Id,
                    longRxnDbIdToName.get(
                            reaction1Id).replace(",", "~"),
                    reaction2Id,
                    longRxnDbIdToName.get(
                            reaction2Id).replace(",", "~"));
    }

    private class CooccurrenceResult {
        List<String> targetReactionList;
        List<String> upstreamReactionList;
        List<Double> pValues;
        List<int[]> ABCDs;
        List<Set<String>[]> BCDSamples;
        List<Set<String>> DSampleFI1s;
        List<Set<String>> DSampleFI2s;
        Map<Double, Double> pValue2FDRMap;

        private CooccurrenceResult(List<String> targetReactionList,
                                   List<String> upstreamReactionList,
                                   List<Double> pValues,
                                   List<int[]> inABCDs,
                                   List<Set<String>[]> inBCDSamples,
                                   List<Set<String>> inDSampleFI1s,
                                   List<Set<String>> inDSampleFI2s) {
            this.targetReactionList = targetReactionList;
            this.upstreamReactionList = upstreamReactionList;
            this.pValues = pValues;
            this.ABCDs = inABCDs;
            this.BCDSamples = inBCDSamples;
            this.DSampleFI1s = inDSampleFI1s;
            this.DSampleFI2s = inDSampleFI2s;
        }

        //TODO: make this faster by sorting all these lists together
        private void CalculateNetworkPermutationFDR(List<CooccurrenceResult> rewiredNetworkResults) {


            List<Double> rewiredNetworkPvalues = new ArrayList<>();

            //loop over rewirings
            for (CooccurrenceResult rewiredNetworkResult : rewiredNetworkResults) {
                rewiredNetworkPvalues.addAll(rewiredNetworkResult.pValues);
            }

            this.pValue2FDRMap = new HashMap<>();

            //loop over pValues
            for (Double pValue : this.pValues) {
                this.pValue2FDRMap.put(pValue,
                        calculateFDR(
                                NumValsBelowThresh(rewiredNetworkPvalues, pValue),
                                NumValsBelowThresh(this.pValues, pValue),
                                rewiredNetworkPvalues.size(),
                                this.pValues.size()));
            }
        }

        private double calculateFDR(int randomLessThan,
                                    int realLessThan,
                                    int randomTotal,
                                    int realTotal) {
            double fdr = ((double) randomLessThan /
                    (double) randomTotal) /
                    ((double) realLessThan /
                            (double) realTotal);
            if (fdr > 1.0)
                fdr = 1.0;
            return fdr;
        }

        private int NumValsBelowThresh(List<Double> vals, Double threshold) {
            int count = 0;
            for (Double val : vals) {
                count = val < threshold ?
                        count + 1 :
                        count;
            }
            return count;
        }

        private void MagicallyShrinkMemoryFootprint() {
            //leave only pValues
            this.targetReactionList = null;
            this.upstreamReactionList = null;
            this.ABCDs = null;
            this.BCDSamples = null;
            this.DSampleFI1s = null;
            this.DSampleFI2s = null;
            this.pValue2FDRMap = null;
            System.gc();
        }
    }

    private class TargetReactionSummary {
        private static final String headerLine =
                "RxnId," +
                        "Num Sup Dn/Up Rxns," +
                        "Num Sup Dn/Up Combos," +
                        "Num All Dn/Up Rxns," +
                        "Support Ratio," +
                        "Sup Dn/Up Rxns," +
                        "Sup Dn/Up Combos," +
                        "All Dn/Up Rxns," +
                        "Supporting FIs";
        private Long rxnId;
        private Integer numSupUpRxn;
        private Integer numUpRxn;
        private Double supUpRatio;
        private Set<Long> supRxns;
        private Set<Long> upRxns;
        private Set<String> supFIs;
        private Map<Long, String> longRxnDbIdToName;

        private TargetReactionSummary(Long rxnId,
                                      Integer numSupUpRxn,
                                      Integer numUpRxn,
                                      Double supUpRatio,
                                      Set<Long> supRxns,
                                      Set<Long> upRxns,
                                      Set<String> supFIs,
                                      Map<Long, String> longRxnDbIdToName) {
            this.rxnId = rxnId;
            this.numSupUpRxn = numSupUpRxn;
            this.numUpRxn = numUpRxn;
            this.supUpRatio = supUpRatio;
            this.supRxns = supRxns;
            this.upRxns = upRxns;
            this.supFIs = supFIs;
            this.longRxnDbIdToName = longRxnDbIdToName;
        }

        @Override
        public boolean equals(Object o) {
            if (o == this) return true;
            if (!(o instanceof TargetReactionSummary)) {
                return false;
            }
            TargetReactionSummary targetReactionSummary = (TargetReactionSummary) o;
            return Objects.equals(this.rxnId, targetReactionSummary.rxnId) &&
                    Objects.equals(this.numSupUpRxn, targetReactionSummary.numSupUpRxn) &&
                    Objects.equals(this.numUpRxn, targetReactionSummary.numUpRxn) &&
                    Objects.equals(this.supUpRatio, targetReactionSummary.supUpRatio) &&
                    Objects.equals(this.supRxns, targetReactionSummary.supRxns) &&
                    Objects.equals(this.upRxns, targetReactionSummary.upRxns) &&
                    Objects.equals(this.supFIs, targetReactionSummary.supFIs);
        }

        @Override
        public int hashCode() {
            return Objects.hash(this.rxnId,
                    this.numSupUpRxn,
                    this.numUpRxn,
                    this.supUpRatio,
                    this.supRxns,
                    this.supFIs);
        }

        /*
        @Override
        public String toString() {
            List<String> supCombos = ConvertReactionPairToString(
                    new ArrayList<>(GenerateReactionSetPairs(this.supRxns)),
                    this.longRxnDbIdToName);
            return String.format(
                    "%d," + //rxnId
                            "%d," + //numSupUpRxn
                            "%d," + //numSupUpCombos
                            "%d," + //numUpRxn
                            "%f," + //supUpRatio
                            "%s," + //supRxns
                            "%s," + //supCombos
                            "%s," + //upRxns
                            "%s",   //supFIs
                    this.rxnId,
                    this.numSupUpRxn,
                    supCombos.size(),
                    this.numUpRxn,
                    this.supUpRatio,
                    this.supRxns.size() > 1
                            ? org.gk.util.StringUtils.join("~",
                            new ArrayList<>(this.supRxns))
                            : Arrays.asList(this.supRxns).get(0),
                    supCombos.size() > 1
                            ? org.gk.util.StringUtils.join("~",
                            new ArrayList<>(supCombos))
                            : Arrays.asList(supCombos).get(0),
                    this.upRxns.size() > 1
                            ? org.gk.util.StringUtils.join("~",
                            new ArrayList<>(this.upRxns))
                            : Arrays.asList(this.upRxns).get(0),
                    this.supFIs.size() > 1
                            ? org.gk.util.StringUtils.join("~",
                            new ArrayList<>(this.supFIs))
                            : Arrays.asList(this.supFIs).get(0));
        }
        */
    }
}
