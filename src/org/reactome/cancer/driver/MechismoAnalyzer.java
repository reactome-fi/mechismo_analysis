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

import javax.swing.text.html.HTMLDocument;
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
    private String dirName = "datasets/Mechismo/";
    private String outputFileName = dirName + "COSMICv74_somatic_noSNPs_GWS_mechismo_output.tsv";
    private String pciContactFile = dirName + "human_pci_contact_hits.tsv";
    private final static Logger logger = Logger.getLogger(MechismoOutputLoader.class);
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

    public void mapReactomeReactions(CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer,
                                     String mechismoOutputFilePath,
                                     String reactomeReactionNetworkFilePath,
                                     String outputDir) throws Exception {
        //build jgrapht network from reaction file
        ReactomeReactionGraphLoader reactomeReactionGraphLoader = new ReactomeReactionGraphLoader(reactomeReactionNetworkFilePath);
        DefaultDirectedGraph<Long, DefaultEdge> reactionGraph = reactomeReactionGraphLoader.getReactionGraph();

        System.out.println(String.format("Created reaction graph with %d vertices and %d edges",
                reactionGraph.vertexSet().size(),
                reactionGraph.edgeSet().size()));

        //extract FIs from Mechismo output
        //Set<FI>
        Set<String> fis = new MechismoOutputLoader(mechismoOutputFilePath).extractMechismoFIs();

        System.out.println(String.format("Loaded %d Mechismo FIs",
                fis.size()));

        //find reactions containing Mechismo FIs
        //FI -> Set<Reaction Id>
        Map<String, Set<Long>> fi2ReactionSet = new HashMap<>();
        //Reaction Id -> Set<FIs>
        Map<Long, Set<String>> reaction2FiSet = new HashMap<>();
        Set<Long> rxns = reactomeReactionGraphLoader.getReactionSet();
        Iterator<Long> rxnItr = rxns.iterator();
        int itCounter = 0;
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        Long rxnDbId;
        System.out.println("Processing Reactions...");
        System.out.flush();
        while (rxnItr.hasNext()) {
            rxnDbId = rxnItr.next();
            Set<String> rxnInteractions = new HashSet<>();
            reactomeAnalyzer.generateFIsForReactionsWithFeatures(cancerDriverReactomeAnalyzer.getDBA(),
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
                System.out.println(String.format("Processed %d of %d Reactions...",
                        itCounter, rxns.size()));
                System.out.flush();
            }
        }

        System.out.println(String.format("Mapped %d FIs to %d reactions",
                fi2ReactionSet.keySet().size(),
                reaction2FiSet.keySet().size()));

        //query jgrapht network to search 1 reaction downstream/upstream
        Set<TargetReactionSummary> targetReactionSummaries = new HashSet<>();
        Set<Set<String>> allFIs = new HashSet<>();
        //foreach SupReaction in reaction2FiSet
        Iterator rxn2FiItr = reaction2FiSet.entrySet().iterator();
        while (rxn2FiItr.hasNext()) {
            Map.Entry pair = (Map.Entry) rxn2FiItr.next();
            Long rxnId = (Long) pair.getKey();
            //http://jgrapht.org/javadoc/org/jgrapht/graph/DefaultDirectedGraph.html
            //go one DnReaction downstream from SupReaction and
            Set<DefaultEdge> dnEdges = reactionGraph.outgoingEdgesOf(rxnId);
            Iterator<DefaultEdge> dnEdgeItr = dnEdges.iterator();
            while (dnEdgeItr.hasNext()) {
                DefaultEdge dnEdge = dnEdgeItr.next();
                Long dnVtx = reactionGraph.getEdgeTarget(dnEdge);
                //one UpReaction from DnReaction
                Set<DefaultEdge> upEdges = reactionGraph.incomingEdgesOf(dnVtx);
                Set<Long> dnUpRxns = new HashSet<>();
                Iterator<DefaultEdge> upEdgeItr = upEdges.iterator();
                while (upEdgeItr.hasNext()) {
                    DefaultEdge upEdge = upEdgeItr.next();
                    Long upVtx = reactionGraph.getEdgeSource(upEdge);
                    dnUpRxns.add(upVtx);
                }
                //find intersection of (dnUpReaction - SupReaction) and reaction2FiSet keys
                Set<Long> upSupIntsct = new HashSet(dnUpRxns);
                upSupIntsct.retainAll(reaction2FiSet.keySet());
                if (!upSupIntsct.contains(rxnId)) {
                    throw new IllegalStateException(
                            String.format("Initial rxn ID '%d' not in downstream upstream rxns",
                                    rxnId)
                    );
                }
                Set<String> supFIs = new HashSet<>();
                Iterator<Long> upSupIntsctItr = upSupIntsct.iterator();
                while (upSupIntsctItr.hasNext()) {
                    Long upSupRxnId = upSupIntsctItr.next();
                    supFIs.addAll(reaction2FiSet.get(upSupRxnId));
                }
                allFIs.add(supFIs);

                targetReactionSummaries.add(new TargetReactionSummary(
                        dnVtx,
                        upSupIntsct.size(),
                        dnUpRxns.size(),
                        new Double(upSupIntsct.size()) / new Double(dnUpRxns.size()),
                        upSupIntsct,
                        dnUpRxns,
                        supFIs));
            }
        }

        System.out.println(String.format("Found %d supported Dn/Up reactions",
                targetReactionSummaries.size()));

        System.out.println(String.format("Found %d unique supporting FI sets",
                allFIs.size()));

        Set<Set<String>> fiIntersectingSetUnionClusters = new HashSet<>();
        Iterator<Set<String>> allFIsItr = allFIs.iterator();
        while(allFIsItr.hasNext()){
            Set<String> fiSet = allFIsItr.next();
            Set<String> union = new HashSet<>(fiSet);
            Iterator<Set<String>> fiClusterItr = fiIntersectingSetUnionClusters.iterator();
            while(fiClusterItr.hasNext()){
                Set<String> fiCluster = fiClusterItr.next();
                Set<String> intersection = new HashSet<>(fiSet);
                intersection.retainAll(fiCluster);
                if (!intersection.isEmpty()){
                    union.addAll(fiCluster);
                    fiClusterItr.remove(); //remove out-dated cluster
                }
            }
            fiIntersectingSetUnionClusters.add(union);
        }

        System.out.println(String.format("Aggregated %d unique FI sets into %d intersection-unions",
                allFIs.size(),
                fiIntersectingSetUnionClusters.size()));

        Iterator<Set<String>> fiClusterItr = fiIntersectingSetUnionClusters.iterator();
        Double min = new Double(fiClusterItr.next().size());
        Double max = min;
        Double sum = max;
        while(fiClusterItr.hasNext())
        {
            Double sz = new Double(fiClusterItr.next().size());
            min = sz < min ? sz : min;
            max = sz > max ? sz : max;
            sum += sz;
            System.out.println(sz);
        }
        Double mean = sum / new Double(fiIntersectingSetUnionClusters.size());

        System.out.println(String.format("Clusters min = %f, max = %f, mean = %f, sum = %f",
                min,
                max,
                mean,
                sum));

        //perform pathway enrichment analysis among FI clusters
        FileUtility fileUtility0 = new FileUtility();
        String outFilePath0 = outputDir + "fiClusters.csv";
        try{
            fileUtility0.setOutput(outFilePath0);
            Iterator<Set<String>> fiClusterItr0 =fiIntersectingSetUnionClusters.iterator();
            while(fiClusterItr0.hasNext()) {
                Set<String> fis0 = fiClusterItr0.next();
                fileUtility0.printLine(
                        fis0.size() > 1
                            ? org.gk.util.StringUtils.join(" ",
                                new ArrayList(fis0)).replace("\t"," ")
                            : Arrays.asList(fis0).get(0).toString()
                                .replace("[","")
                                .replace("]","")
                .replace("\t"," "));
            }
            fileUtility0.close();
        }catch(IOException ioe){
            logger.error(String.format("Couldn't use %s, %s: %s",
                    outFilePath0.toString(),
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }

        //write targetReactionSummaries to file
        FileUtility fileUtility = new FileUtility();
        String outFilePath = outputDir + "dnUpReactions.csv";
        try {
            fileUtility.setOutput(outFilePath);
            fileUtility.printLine(TargetReactionSummary.headerLine);
            Iterator<TargetReactionSummary> trsItr = targetReactionSummaries.iterator();
            while (trsItr.hasNext()) {
                TargetReactionSummary targetReactionSummary = trsItr.next();
                fileUtility.printLine(targetReactionSummary.toString());
            }
            fileUtility.close();
        }catch(IOException ioe){
            logger.error(String.format("Couldn't use %s, %s: %s",
                    outFilePath.toString(),
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }

        System.out.println(String.format("Output written to %s",
                outFilePath));
    }

    private class TargetReactionSummary {
        private Long rxnId;
        private Integer numSupUpRxn;
        private Integer numUpRxn;
        private Double supUpRatio;
        private Set<Long> supRxns;
        private Set<Long> upRxns;
        private Set<String> supFIs;
        private static final String headerLine =
                "RxnId," +
                        "Num Sup Dn/Up Rxns," +
                        "Num All Dn/Up Rxns," +
                        "Support Ratio," +
                        "Sup Dn/Up Rxns," +
                        "All Dn/Up Rxns," +
                        "Supporting FIs";

        private TargetReactionSummary(Long rxnId,
                                      Integer numSupUpRxn,
                                      Integer numUpRxn,
                                      Double supUpRatio,
                                      Set<Long> supRxns,
                                      Set<Long> upRxns,
                                      Set<String> supFIs) {
            this.rxnId = rxnId;
            this.numSupUpRxn = numSupUpRxn;
            this.numUpRxn = numUpRxn;
            this.supUpRatio = supUpRatio;
            this.supRxns = supRxns;
            this.upRxns = upRxns;
            this.supFIs = supFIs;
        }

        @Override
        public boolean equals(Object o){
            if(o == this) return true;
            if(!(o instanceof TargetReactionSummary)){
                return false;
            }
            TargetReactionSummary targetReactionSummary = (TargetReactionSummary)o;
            return Objects.equals(this.rxnId,targetReactionSummary.rxnId) &&
                    Objects.equals(this.numSupUpRxn,targetReactionSummary.numSupUpRxn) &&
                    Objects.equals(this.numUpRxn,targetReactionSummary.numUpRxn) &&
                    Objects.equals(this.supUpRatio,targetReactionSummary.supUpRatio) &&
                    Objects.equals(this.supRxns,targetReactionSummary.supRxns) &&
                    Objects.equals(this.upRxns,targetReactionSummary.upRxns) &&
                    Objects.equals(this.supFIs,targetReactionSummary.supFIs);
        }

        @Override
        public int hashCode(){
            return Objects.hash(this.rxnId,
                    this.numSupUpRxn,
                    this.numUpRxn,
                    this.supUpRatio,
                    this.supRxns,
                    this.supFIs);
        }

        @Override
        public String toString() {
            return String.format(
                    "%d," + //rxnId
                            "%d," + //numSupUpRxn
                            "%d," + //numUpRxn
                            "%f," + //supUpRatio
                            "%s," + //supRxns
                            "%s," + //upRxns
                            "%s",   //supFIs
                    this.rxnId,
                    this.numSupUpRxn,
                    this.numUpRxn,
                    this.supUpRatio,
                    this.supRxns.size() > 1
                            ? org.gk.util.StringUtils.join("~",
                                new ArrayList(this.supRxns))
                            : Arrays.asList(this.supRxns).get(0),
                    this.upRxns.size() > 1
                            ? org.gk.util.StringUtils.join("~",
                                new ArrayList(this.upRxns))
                            : Arrays.asList(this.upRxns).get(0),
                    this.supFIs.size() > 1
                            ? org.gk.util.StringUtils.join("~",
                            new ArrayList(this.supFIs))
                            : Arrays.asList(this.supFIs).get(0));
        }
    }
}
