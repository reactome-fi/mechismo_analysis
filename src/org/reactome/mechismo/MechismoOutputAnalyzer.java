package org.reactome.mechismo;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.gk.model.GKInstance;
import org.junit.Test;
import org.reactome.cancer.driver.CancerDriverAnalyzer;
import org.reactome.cancer.driver.CancerDriverReactomeAnalyzer;
import org.reactome.r3.Interactome3dAnalyzer;
import org.reactome.r3.ReactionMapGenerator;
import org.reactome.r3.ReactomeAnalyzer;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.Plotter;

public class MechismoOutputAnalyzer {
    private String dirName = "datasets/Mechismo/";
    private String outputFileName = dirName + "COSMICv74_somatic_noSNPs_GWS_mechismo_output.tsv";
    private String pciContactFile = dirName + "human_pci_contact_hits.tsv";
    private FileUtility fu = new FileUtility();
    
    public MechismoOutputAnalyzer() {
    }
    
    /**
     * Look through all PPIs subject to mechismo analysis. Some of them
     * may don't have any mutations and not reported in Francesco's results.
     * @throws Exception
     */
    @Test
    public void checkFIsSubjectToAnalysis() throws Exception {
        try (Stream<String> stream = Files.lines(Paths.get(dirName, "TCGA", "073018", "m2b_contact_hits.tsv"))) {
            Set<String> confidenceLevels = new HashSet<>();
            Set<String> ppis = stream.skip(1)
                              .map(line -> line.split("\t"))
                              .filter(tokens -> tokens.length > 5)
                              .filter(tokens -> tokens[1].equals("PPI")) // PPI only
                              .filter(tokens -> !tokens[3].equals(tokens[5])) // Not self interaction
//                              .filter(tokens -> tokens[14].equals("high")) // Filter to high only
                              .map(tokens -> {
                                  confidenceLevels.add(tokens[14]);
                                  String protein1 = tokens[3];
                                  int index = protein1.indexOf("-");
                                  if (index > 0)
                                      protein1 = protein1.substring(0, index); // Don't want to count isoforms now
                                  String protein2 = tokens[5];
                                  index = protein2.indexOf("-");
                                  if (index > 0)
                                      protein2 = protein2.substring(0, index);
                                  return InteractionUtilities.generateFIFromGene(protein1, protein2);
                              })
                              .collect(Collectors.toSet());
            System.out.println("Confidence levels: " + confidenceLevels);
            System.out.println("Total ppis: " + ppis.size());
//            int count = 0;
//            for (String ppi : ppis) {
//                System.out.println(ppi);
//                count ++;
//                if (count == 100)
//                    break;
//            }
            
            FINetworkAnalzyer fiAnalyzer = new FINetworkAnalzyer();
            Set<String> fis = fiAnalyzer.loadReactionFIsInProteins();
            System.out.println("Total Reactome FIs: " + fis.size());
            
            Set<String> sharedFIs = InteractionUtilities.getShared(fis, ppis);
            System.out.println("Total shared: " + sharedFIs.size());
            
            // Output shared FIs
            String outFileName = "results/ProteinFIsInReactions_073118_Mechismo.txt";
            fiAnalyzer.filterReactionFIs(sharedFIs, outFileName);
        }
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
//              "AKT1"
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
//                  System.out.println();
//        }
//        fu.close();
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

}
