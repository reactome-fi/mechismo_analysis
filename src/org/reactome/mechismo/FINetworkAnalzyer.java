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
import java.util.StringJoiner;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.junit.Test;
import org.reactome.cancer.driver.MechismoBFSHelper;
import org.reactome.fi.util.FileUtility;
import org.reactome.mechismo.model.AnalysisResult;
import org.reactome.mechismo.model.CancerType;
import org.reactome.mechismo.model.Interaction;
import org.reactome.mechismo.model.Sample;
import org.reactome.mechismows.MechismowsReader;
import org.reactome.r3.graph.BreadthFirstSearch;
import org.reactome.r3.graph.GraphAnalyzer;
import org.reactome.r3.graph.NetworkBuilderForGeneSet;
import org.reactome.r3.util.InteractionUtilities;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

/**
 * This class is used to analyze the mechismo results in the perspective of the 
 * Reactome FI network, composed FIs extracted from Reactome only.
 * @author wug
 *
 */
public class FINetworkAnalzyer {
    // As of August, 2018, this cutoff is used for pancancer data only.
    private final Double FDR_CUTOFF = 0.05d;
    private final int MUTATION_COUNT_CUTOFF = 5;
    private final FileUtility fu = new FileUtility();

    public FINetworkAnalzyer() {
    }
    
    private Set<String> loadGenesInSignificantReactomeFIs() throws IOException {
        FileUtility fu = new FileUtility();
        String file = "results/ReactomeSignificantFIs_082118.txt";
        Set<String> genes = new HashSet<>();
        fu.setInput(file);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            genes.add(tokens[2]);
            genes.add(tokens[3]);
        }
        fu.close();
        return genes;
    }
    
    /**
     * Load genes from Francesco's summary file.
     * @return
     * @throws IOException
     */
    private Set<String> loadGenesInSignificantReactomeFIsInSummary() throws IOException {
        String file = "results/tcga_mechismo_stat_cancer_wise_interfaces_sig_021921.tsv";
        Set<String> genes = new HashSet<>();
        fu.setInput(file);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            genes.add(tokens[1]);
            genes.add(tokens[2]);
        }
        fu.close();
        return genes;
    }
    
    /**
     * Check cancers and significant FIs assigned to cancers.
     * @throws IOException
     */
    @Test
    public void checkCancerToSigFIs() throws IOException {
        Map<String, Set<String>> cancerToFIs = loadCancerToSigFIs();
        List<String> cancers = new ArrayList<>(cancerToFIs.keySet());
        cancers.sort((c1, c2) -> cancerToFIs.get(c1).size() - cancerToFIs.get(c2).size());
        System.out.println("Cancer\tNumber\tSignificantFIs");
        cancers.forEach(key -> {
            Set<String> set = cancerToFIs.get(key);
            List<String> list = new ArrayList<>(set);
            list.sort(Comparator.naturalOrder());
            System.out.println(key + "\t" + set.size() + "\t" + String.join(", ", list));
        });
        
        String outFileName = "results/SigFIs_based_Cancer_Dist_100818.txt";
        
        generatePairWiseNetworkDistance(cancerToFIs, outFileName);
    }

    protected void generatePairWiseNetworkDistance(Map<String, Set<String>> sampleToFIs,
                                                   String outFileName) throws IOException {
        // Calculate pair-wise distances between cancers
        // Use a helper object
        String networkFile = "results/FINodesNetwork_First_Component_010719.txt";
        Set<String> network = fu.loadInteractions(networkFile);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> idToPartners = bfs.generateIdToPartnersMap(network);
        // Need a little bit filter
        sampleToFIs.forEach((cancer, fis) -> fis.retainAll(idToPartners.keySet()));
        // Check if there is any cancer having no FIs after the above filtering
        System.out.println("\nThe following cancers don't have FIs in the largest component:");
        for (Iterator<String> it = sampleToFIs.keySet().iterator(); it.hasNext();) {
            String sample = it.next();
            Set<String> fis = sampleToFIs.get(sample);
            if (fis.size() == 0) {
                System.out.println(sample);
                it.remove();
            }
        }
        
        MechismoBFSHelper helper = new MechismoBFSHelper();
        System.out.println("\nCalculating pair-wise distances among cancers:");
        Map<String, Integer> pairToDist = helper.calculateShortestPath(sampleToFIs,
                                                                       idToPartners,
                                                                       bfs);
        List<String> samples = sampleToFIs.keySet().stream().sorted().collect(Collectors.toList());
        StringBuilder builder = new StringBuilder();
        builder.append("Sample");
        samples.forEach(cancer -> builder.append("\t").append(cancer));
        System.out.println(builder.toString());
        fu.setOutput(outFileName);
        fu.printLine(builder.toString());
        builder.setLength(0);
        for (int i = 0; i < samples.size(); i++) {
            String cancer = samples.get(i);
            builder.append(cancer);
            Set<String> set1 = sampleToFIs.get(cancer);
            for (int j = 0; j < samples.size(); j++) {
                String cancer2 = samples.get(j);
                builder.append("\t");
                if (i == j)
                    builder.append(0.0d);
                else {
                    Set<String> set2 = sampleToFIs.get(cancer2);
                    double dist = helper.calculateMinShortestPath(set1,
                            set2,
                            pairToDist);
                    builder.append(dist);
                }
            }
            System.out.println(builder.toString());
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    @Test
    public void calculateSampleDistancesOnFINodeNetwork() throws IOException {

        String cancer = "GBM";
        cancer = "PAAD";
//        cancer = "SKCM";
                cancer = "LGG";
                cancer = "UCEC";
//        cancer = "COADREAD";
                cancer = "THCA";
                cancer = "STAD";
                cancer = "LUAD";
                cancer = "HNSC";

        // Generate distances based on significant FIs for individual cancers
        String srcFileName = "datasets/Mechismo/FrancescoResults/080918/tcga_mechismo_stat_cancer_wise_undirected.tsv";
        double fdrCutoff = 0.05d;
        boolean hasHeader = true;

        Map<String, Set<String>> sampleToFIs = loadSampleToFIs(srcFileName,
                                                               cancer,
                                                               fdrCutoff,
                                                               hasHeader);

        // A quick quality check
        sampleToFIs.forEach((sample, set) -> System.out.println(sample + "\t" + set.size()));

        //        filterSamplesWithoutReactions(sampleToReactions);

        String date = "010719";

        String fileName = "results/MechismoSamplePairWiseFINodeNetworkDist_" + cancer + "_" + date + ".txt";

        generatePairWiseNetworkDistance(sampleToFIs, fileName);
    }
    
    
    @Test
    public void testLoadSamplesToFIs() throws IOException {
        String fileName = "datasets/Mechismo/FrancescoResults/080918/tcga_mechismo_stat_cancer_wise_undirected.tsv";
        String cancer = "COADREAD";
        double cutoff = 0.05d;
        Map<String, Set<String>> sampleToFIs = loadSampleToFIs(fileName, cancer, cutoff, true);
        sampleToFIs.forEach((sample, fis) -> System.out.println(sample + "\t" + fis.size() + "\t" + fis));
    }
    
    /**
     * Use this method to load sample to hit FIs.
     *
     * @param cancer
     * @param fdrCutoff
     * @return
     * @throws IOException
     */
    private Map<String, Set<String>> loadSampleToFIs(String srcFileName,
                                                     String cancer,
                                                     double fdrCutoff,
                                                     boolean hasHeader) throws IOException {
        Map<String, Set<String>> sampleToFIs = new HashMap<>();
        try (Stream<String> stream = Files.lines(Paths.get(srcFileName))) {
            stream.skip(hasHeader ? 1 : 0)
                    .map(line -> line.split("\t"))
                    .filter(tokens -> tokens[0].equals(cancer)) // Filter to cancer type
                    .filter(tokens -> !tokens[13].equals("-")) // Need to have reactions
                    .filter(tokens -> Double.parseDouble(tokens[11]) < fdrCutoff) // Less than the predefined fdr cutoff
                    .forEach(tokens -> {
                        // Get FI
                        String fi = InteractionUtilities.generateFIFromGene(tokens[1], tokens[2]);
                        String fi1 = fi.replace('\t', ' ');
                        // Get samples
                        int index1 = tokens[15].indexOf("\'");
                        int index2 = tokens[15].lastIndexOf("\'");
                        String tmp = tokens[15].substring(index1, index2 + 1);
                        String[] samples = tmp.split(", ");
                        Arrays.asList(samples).forEach(sample -> {
                            // Need to remove single quotes
                            sample = sample.substring(1, sample.length() - 1);
                            sampleToFIs.compute(sample, (key, set) -> {
                                if (set == null)
                                    set = new HashSet<>();
                                set.add(fi1);
                                return set;
                            });
                        });
                    });
        }
        return sampleToFIs;
    }
    
    public Map<String, Set<String>> loadCancerToSigFIs() throws IOException {
        FileUtility fu = new FileUtility();
        String file = "results/ReactomeSignificantFIs_082118.txt";
        Map<String, Set<String>> cancerToFIs = new HashMap<>();
        fu.setInput(file);
        String line = fu.readLine();
        String[] header = line.split("\t");
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String fi = InteractionUtilities.generateFIFromGene(tokens[2], tokens[3]);
            String fi1 = fi.replaceAll("\t", " ");
            for (int i = 9; i < tokens.length; i++) {
                String token = tokens[i];
                if (token.length() == 0)
                    continue;
                Double fdr = new Double(token);
                if (fdr >= FDR_CUTOFF)
                    continue;
                String cancer = header[i];
                cancerToFIs.compute(cancer, (key, set) -> {
                    if (set == null)
                        set = new HashSet<>();
                    set.add(fi1);
                    return set;
                });
            }
        }
        fu.close();
        return cancerToFIs;
    }
    
    /**
     * Perform some check on the network constructed by using FIs as nodes and edges if two FIs having
     * shared genes. See method convertGeneNetworkToFINetwork().
     * @throws Exception
     */
    @Test
    public void checkFINodeNetwork() throws Exception {
        String fileName = "results/FINodesNetwork_100818.txt";
        Set<String> edges = new HashSet<>();
        Set<String> nodes = new HashSet<>();
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            line = line.replaceAll("\t", " ");
            String[] tokens = line.split(",");
            edges.add(tokens[0] + "\t" + tokens[1]);
            nodes.add(tokens[0]);
            nodes.add(tokens[1]);
        }
        fu.close();
        System.out.println("Total edges: " + edges.size());
        System.out.println("Total nodes: " + nodes.size());
        
        GraphAnalyzer graphAnalyzer = new GraphAnalyzer();
        List<Set<String>> components = graphAnalyzer.calculateGraphComponents(edges);
        System.out.println("Total components: " + components.size());
        components.forEach(comp -> System.out.println(comp.size()));
        // Output the first component
        String outFileName = "results/FINodesNetwork_First_Component_100818.txt";
        // All FIs are sorted
        outFileName = "results/FINodesNetwork_First_Component_010719.txt";
        fu.setOutput(outFileName);
        Set<String> firstCompNodes = components.get(0);
        for (String edge : edges) {
            String[] fis = edge.split("\t");
            if (firstCompNodes.contains(fis[0])) {// We need edges
                String fi1 = sortFI(fis[0]);
                String fi2 = sortFI(fis[1]);
                String sortedEdge = InteractionUtilities.generateFIFromGene(fi1, fi2);
                fu.printLine(sortedEdge);
            }
        }
        fu.close();
    }
    
    private String sortFI(String fi) {
        String[] tokens = fi.split(" ");
        return Stream.of(tokens).sorted().collect(Collectors.joining(" "));
    }
    
    /**
     * This method is used to generate a connected FI subnetwork contains all significant FIs.
     * @throws IOException
     */
    @Test
    public void generateFINetworkForSignificantFIs() throws IOException {
        Set<String> significantGenes = loadGenesInSignificantReactomeFIsInSummary();
        System.out.println("Total significant genes: " + significantGenes.size());
        Set<String> fis = loadReactomeFIs();
        System.out.println("Total FIs: " + fis.size());
        NetworkBuilderForGeneSet builder = new NetworkBuilderForGeneSet();
        builder.setAllFIs(fis);
        // Work on a sub-network
        // Since the original FI network may have multiple components, the generated
        // sub-network may have multiple components.
        fis = builder.constructFINetworkForGeneSet(significantGenes);
        System.out.println("Size of sub-network: " + fis.size());
        // Output the FI network
        String outputFileName = "results/SignificantFIsNetwork_022721.txt";
        fu.saveInteractions(fis, outputFileName);
        // The following genes are not in the network
        Set<String> genesInFIs = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Total genes in the constructed network: " + genesInFIs.size());
        significantGenes.removeAll(genesInFIs);
        System.out.println("Significant genes not in the network: " + significantGenes.size());
        significantGenes.stream().sorted().forEach(System.out::println);
    }
    
    @Test
    public void selectFIsInReactionForCancer() throws IOException {
        String srcFileName = "results/tcga_mechismo_stat_cancer_wise_reactions_sig_CGC_021921.tsv";
        String cancer = "STAD";
        fu.setInput(srcFileName);
        String line = fu.readLine();
        String output = "results/" + cancer + "_FIsInSignificantReactions_030121.txt";
        Set<String> fis = new HashSet<>();
        Map<String, Set<String>> fiToRxts = new HashMap<>();
        Set<String> genes = new HashSet<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (!tokens[0].equals(cancer))
                continue;
            String[] tokens1 = tokens[15].split(",");
            for (String token : tokens1) {
                String[] tokens2 = token.split(":|\\|");
                String fi = InteractionUtilities.generateFIFromGene(tokens2[0], tokens2[1]);
                genes.add(tokens2[0]);
                genes.add(tokens2[1]);
                fiToRxts.compute(fi, (key, set) -> {
                    if (set == null)
                        set = new HashSet<>();
                    set.add(tokens[1]);
                    return set;
                });
                fi = fi.replace("\t", "\tFI\t");
                fis.add(fi);
            }
        }
        fu.close();
        fu.saveInteractions(fis, output);
        System.out.println("Total selected FIs: " + fis.size());
        // Build a connected network
        Set<String> allFIs = loadReactomeFIs();
        System.out.println("Total FIs: " + allFIs.size());
        NetworkBuilderForGeneSet builder = new NetworkBuilderForGeneSet();
        builder.setAllFIs(allFIs);
        // Work on a sub-network
        // Since the original FI network may have multiple components, the generated
        // sub-network may have multiple components.
        fis = builder.constructFINetworkForGeneSet(genes);
        System.out.println("Size of sub-network: " + fis.size());
        // The returned FIs are not sorted in individual FIs
        fis = fis.stream()
                 .map(fi -> fi.split("\t"))
                 .map(tokens -> InteractionUtilities.generateFIFromGene(tokens[0], tokens[1]))
                 .collect(Collectors.toSet());
        System.out.println("Re-formated FIs: " + fis.size());
        // Output the FI network
        String outputFileName = "results/" + cancer + "_FIsInSignificantReactions_Connected_030121.txt";
        fu.setOutput(outputFileName);
        fu.printLine("Protein1\tType\tProtein2\tSigReactionNumber\tSigReactions");
        for (String fi : fis) {
            String[] tokens = fi.split("\t");
            Set<String> rxts = fiToRxts.get(fi);
            fu.printLine(tokens[0] + "\tFI\t" + tokens[1] + "\t" + 
                         (rxts == null ? "0" : rxts.size()) + "\t" + 
                         (rxts == null ? "" : String.join(",", rxts)));
        }
        fu.close();
    }
    
    @Test
    public void selectSignficantReactionFIsInCancer() throws IOException {
        String srcFileName = "results/tcga_mechismo_stat_cancer_wise_interfaces_sig_021921.tsv";
        String cancer = "STAD";
        fu.setInput(srcFileName);
        String output = "results/" + cancer + "_ReactomeSigFIs_021921.txt";
        fu.setOutput(output);
        String line = fu.readLine();
        Set<String> genes = new HashSet<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (!tokens[0].equals(cancer) ||
                tokens[13].equals("-")) // Not from Reaction
                continue; 
            fu.printLine(tokens[1] + "\t" + tokens[2]);
            String fi = InteractionUtilities.generateFIFromGene(tokens[1], tokens[2]);
            fi = fi.replace("\t", " (FI) ");
            System.out.println(fi + "\t1");
            if (!tokens[1].startsWith("[)"))
                genes.add(tokens[1]);
            if (!tokens[2].startsWith("["))
                genes.add(tokens[2]);
        }
        fu.close();
        System.out.println("Selected genes: " + genes.size());
        if (true)
            return;
        // Build a connected network
        Set<String> fis = loadReactomeFIs();
        System.out.println("Total FIs: " + fis.size());
        NetworkBuilderForGeneSet builder = new NetworkBuilderForGeneSet();
        builder.setAllFIs(fis);
        // Work on a sub-network
        // Since the original FI network may have multiple components, the generated
        // sub-network may have multiple components.
        fis = builder.constructFINetworkForGeneSet(genes);
        System.out.println("Size of sub-network: " + fis.size());
        // Output the FI network
        String outputFileName = "results/" + cancer + "_ReactomeSigFIs_Connected_021921.txt";
        fu.saveInteractions(fis, outputFileName);
        
        // For selection
        genes.stream().sorted().forEach(System.out::println);
    }
    
    /**
     * Convert from a gene-based network (aka genes are nodes) to FI-based network (aka FIs 
     * are nodes)
     * Note: A huge network file will be generated by using all FIs in the Reactome FI network. Instead,
     * focus constructing a linked network for genes involved in the significant FIs collected
     * in file, ReactomeSignificantFIs_082118.txt, and then try to link them together. 
     * @throws IOException
     */
    @Test
    public void convertGeneNetworkToFINetwork() throws IOException {
//        Set<String> significantGenes = loadGenesInSignificantReactomeFIs();
        Set<String> significantGenes = loadGenesInSignificantReactomeFIsInSummary();
        System.out.println("Total significant genes: " + significantGenes.size());
        Set<String> fis = loadReactomeFIs();
        System.out.println("Total FIs: " + fis.size());
        NetworkBuilderForGeneSet builder = new NetworkBuilderForGeneSet();
        builder.setAllFIs(fis);
        // Work on a sub-network
        // Since the original FI network may have multiple components, the generated
        // sub-network may have multiple components.
        fis = builder.constructFINetworkForGeneSet(significantGenes);
        System.out.println("Size of sub-network: " + fis.size());
        // Cache everything for quick performance
        Map<String, String[]> fiToGenes = new HashMap<>();
        fis.forEach(fi -> {
            String[] tokens = fi.split("\t");
            fiToGenes.put(fi, tokens);
        });
        System.out.println("Finished spliting!");
        List<String> fiList = new ArrayList<>(fis);
        Collections.sort(fiList);
//        String outFileName = "results/FINodesNetwork_100818.txt";
        String outFileName = "results/FINodesNetwork_022721.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(outFileName);
        for (int i = 0; i < fiList.size() - 1; i++) {
            String fi1 = fiList.get(i);
            String[] genes1 = fiToGenes.get(fi1);
            for (int j = i + 1; j < fiList.size(); j++) {
                String fi2 = fiList.get(j);
                String[] genes2 = fiToGenes.get(fi2);
                if (hasSharedGenes(genes1, genes2))
                    fu.printLine(fi1 + "," + fi2);
            }
        }
        fu.close();
    }
    
    private boolean hasSharedGenes(String[] genes1, String[] genes2) {
        for (String gene1 : genes1) {
            for (String gene2 : genes2) {
                if (gene1.equals(gene2))
                    return true;
            }
        }
        return false;
    }

    private Set<String> loadReactomeFIs() throws IOException {
//        String fileName = "results/ProteinFIsInReactions_032017.txt";
//        String fileName = "results/ProteinFIsInReactions_073118.txt";
        // This is the file used to generate tcga_mechismo_stat_cancer_wise_interfaces_sig_021921.tsv by
        // Francesco (N.B. by GW on Feb 27, 2021)
        String fileName = "results/ProteinFIsInReactionsWithPPIEvidence_032017.txt";
        try (Stream<String> lines = Files.lines(Paths.get(fileName))) {
            Set<String> fis = lines.map(line -> {
                String[] tokens = line.split("\t");
                if (tokens[2].equals("null") || tokens[3].equals("null"))
                    return null;
//                if (tokens[2].startsWith("UB") || tokens[3].equals("UB")) // Too many FIs
//                    return null;
                String fi = InteractionUtilities.generateFIFromGene(tokens[2], tokens[3]);
                return fi;
            })
                    .filter(fi -> fi != null)
                    .collect(Collectors.toSet());
            return fis;
        }
    }
    
    public void filterReactionFIs(Set<String> targetFIs, String outFileName) throws IOException {
        String inFileName = "results/ProteinFIsInReactions_073118.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(inFileName);
        fu.setOutput(outFileName);
        String line = fu.readLine();
        fu.printLine(line);
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String fi = InteractionUtilities.generateFIFromGene(tokens[0], tokens[1]);
            if (targetFIs.contains(fi))
                fu.printLine(line);
        }
        fu.close();
    }
    
    public Set<String> loadReactionFIsInProteins() throws IOException {
        String fileName = "results/ProteinFIsInReactions_073118.txt";
        try (Stream<String> lines = Files.lines(Paths.get(fileName))) {
            Set<String> fis = lines.skip(1)
                    .map(line -> line.split("\t"))
                    .map(tokens -> InteractionUtilities.generateFIFromGene(tokens[0], tokens[1]))
                    .collect(Collectors.toSet());
            return fis;
        }
    }
    
    @Test
    public void outputSelectedFIs() throws IOException {
        String fileName = "results/FI_Cancer_FDR_Filtered_080918.txt";
        fileName = "results/FI_Cancer_FDR_Filtered_082118.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String header = fu.readLine();
        System.out.println("Name\tEDGE_TYPE\t" + header);
        fu.close();
        try (Stream<String> stream = Files.lines(Paths.get(fileName))) {
            stream.skip(1)
                  .map(line -> line.split("\t")) // Escape the header line
                  .filter(tokens -> tokens[2].equals("1")) // Make sure only Reactome FIs are used
                  .filter(tokens -> tokens.length > 6) // At least one value should be there
                  .filter(tokens -> (Integer.parseInt(tokens[3]) > 0 || (tokens[6].length() > 0 && Double.parseDouble(tokens[6]) < 0.05d))) // Make sure only significant in at least on cancer or pancancer is selected
                  .map(tokens -> {
                      StringJoiner joiner = new StringJoiner("\t");
                      joiner.add(tokens[0] + " (FI) " + tokens[1]);
                      joiner.add("FI");
                      for (String token : tokens)
                          joiner.add(token);
                      return joiner.toString();
                   })
                  .forEach(System.out::println);
        }
    }
    
    /**
     * Generate a tab delimited file having FDR values for individual cancer types.
     * @throws Exception
     */
    @Test
    public void sortInteractions() throws Exception {
        Set<String> reactomeFIs = loadReactomeFIs();
        // Output all FIs in mechismo output file
        String fileName = "results/FI_Cancer_FDR_050918.txt";
        // Output FIs having analysis results attached (no protein/chemicals)
        fileName = "results/FI_Cancer_FDR_Filtered_050918.txt";
        fileName = "results/FI_Cancer_FDR_Filtered_051018.txt";
        fileName = "results/FI_Cancer_FDR_Filtered_073118_1.txt";
        fileName = "results/FI_Cancer_FDR_Filtered_080918.txt";
        fileName = "results/FI_Cancer_FDR_Filtered_082118.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        MechismowsReader reader = new MechismowsReader();
        AnnotationConfigApplicationContext context = reader.getContext();
        List<CancerType> cancerTypes = reader.loadCancerTypes(context);
        System.out.println("Total cancer types: " + cancerTypes.size());
        List<String> sortedCancerNames = sortCancerTypes(cancerTypes);
        List<Interaction> interactions = reader.loadInteractions(context);
        System.out.println("Total interactions: " + interactions.size());
        double smallestFDR = getSmallestNonZeroFDR(interactions);
        double largestPScore = -Math.log10(smallestFDR / 10.0d); // Treat as the largest p-score
        // Perform filtering
        interactions = interactions.stream()
                                   .filter(i -> i.getAnalysisResults() != null) // Make sure analysis results there
                                   .filter(i -> !i.getName().contains("[CHEM:")) // Make sure they are protein FIs
                                   .collect(Collectors.toList());
        fu.printLine(generateHeader(sortedCancerNames));
        interactions.forEach(interaction -> {
            String text = outputResultsAsString(interaction, sortedCancerNames, reactomeFIs, largestPScore);
            if (text == null || text.length() == 0)
                return;
            try {
                fu.printLine(text); 
            }
            catch(IOException e) {
                e.printStackTrace();
            }
        });
        context.close();
        fu.close();
    }
    
    private Double getSmallestNonZeroFDR(List<Interaction> interactions) {
        Double smallest = Double.MAX_VALUE;
        for (Interaction interaction : interactions) {
            Set<AnalysisResult> results = interaction.getAnalysisResults();
            for (AnalysisResult result : results) {
                if (result.getFdr() > 0.0d && result.getFdr() < smallest)
                    smallest = result.getFdr();
            }
        }
        return smallest;
    }
    
    private List<String> sortCancerTypes(List<CancerType> types) {
        List<String> rtn = types.stream()
                                .map(type -> type.getAbbreviation())
                                .sorted((t1, t2) -> {
                                    if (t1.equals("PANCAN"))
                                        return -1;
                                    if (t2.equals("PANCAN"))
                                        return 1;
                                    return t1.compareTo(t2);
                                })
                                .collect(Collectors.toList());
        return rtn;
    }
    
    private String generateHeader(List<String> cancerTypes) {
        StringBuilder builder = new StringBuilder();
        builder.append("Partner1\tPartner2\tInReactome\tSignificantCancers\tRankScore\tPanCancer_PScore");
        cancerTypes.forEach(type -> builder.append("\t").append(type));
        return builder.toString();
    }
    
    /**
     * A rank score is calculated as sum of negative log fdr value. If fdr == 0.0,
     * the largest score will be used. Pancancer value is not calculated here.
     * @param results
     * @return
     */
    private double calculateRankScore(Set<AnalysisResult> results, 
                                      double largestScore, 
                                      Map<String, Integer> cancerToCounts) {
        double score = 0.0d;
        for (AnalysisResult result : results) {
            // Need to escape pancancer
            if (result.getCancerType().getAbbreviation().equals("PANCAN"))
                continue;
            Integer counts = cancerToCounts.get(result.getCancerType().getAbbreviation());
//            if (counts < MUTATION_COUNT_CUTOFF)
//                continue;
            Double fdr = result.getFdr();
            if (fdr.equals(0.0d))
                score += largestScore;
            else
                score -= Math.log(fdr);
        }
        return score;
    }
    
    private int calculateSignificantCancers(Set<AnalysisResult> results,
                                            Map<String, Integer> cancerToCounts) {
        int counter = 0;
        for (AnalysisResult result : results) {
            // Need to escape pancancer
            if (result.getCancerType().getAbbreviation().equals("PANCAN"))
                continue;
            Integer counts = cancerToCounts.get(result.getCancerType().getAbbreviation());
//            if (counts < MUTATION_COUNT_CUTOFF)
//                continue;
            if (result.getFdr() <= 0.05d)
                counter ++;
        }
        return counter;
    }
    
    private double calculatePancanPScore(Map<String, AnalysisResult> typeToResult, 
            double largestPScore,
            Map<String, Integer> cancerToCounts) {
        AnalysisResult result = typeToResult.get("PANCAN");
        if (result == null)
            return 0.0d;
        Integer counts = cancerToCounts.get("PANCAN");
        if (counts < MUTATION_COUNT_CUTOFF)
            return 0.0d;
        Double fdr = result.getFdr();
        if (fdr.equals(0.0d))
            return largestPScore;
        return -Math.log10(fdr);
    }
    
    /**
     * To make the following method work, annotation for Interaction and Mutation classes in mechismows should be modified 
     * as following (use EAGER, instead of LAZY, which is the default).
     * @OneToMany(cascade = CascadeType.ALL, fetch = FetchType.EAGER) // Have to use ALL. Otherwise, mutations cannot be saved automatically for some reason.
     * @JsonIgnoreProperties({"hibernateLazyInitializer", "handler"})
     *      * private Set<Mutation> mutations;
     *     @ManyToMany(cascade= CascadeType.ALL, fetch = FetchType.EAGER)
     * @JsonIgnoreProperties({"hibernateLazyInitializer", "handler"})
     *  private Set<Sample> samples;
     * 
     * @param interaction
     * @return
     */
    private Map<String, Integer> extractCancerToCounts(Interaction interaction) {
        Map<String, Integer> cancerToCounts = new HashMap<>();
        Set<Sample> allSamples = new HashSet<>();
        interaction.getMutations().forEach(mutation -> {
            Set<Sample> samples = mutation.getSamples();
            samples.forEach(sample -> {
                cancerToCounts.compute(sample.getCancerType().getAbbreviation(), (key, value) -> {
                    if (value == null)
                        return 1;
                    return ++ value;
                });
            });
            allSamples.addAll(samples);
        });
        cancerToCounts.put("PANCAN", allSamples.size());
        return cancerToCounts;
    }
    
    private String outputResultsAsString(Interaction interaction,
                                         List<String> cancerTypes,
                                         Set<String> reactomeFIs,
                                         double largestPScore) {
        Set<AnalysisResult> results = interaction.getAnalysisResults();
        if (results == null || results.size() == 0)
            return null;
        Map<String, Integer> cancerToCounts = extractCancerToCounts(interaction);
        // Make sure there is something having enough mutation count
//        boolean isGood = cancerToCounts.values().stream().filter(c -> c >= MUTATION_COUNT_CUTOFF).findAny().isPresent();
//        if (!isGood)
//            return null;
        StringBuilder builder = new StringBuilder();
        builder.append(interaction.getName());
        builder.append("\t").append(reactomeFIs.contains(interaction.getName()) ? 1 : 0);
        int significantCancers = calculateSignificantCancers(results, cancerToCounts);
        builder.append("\t").append(significantCancers);
        double rankScore = calculateRankScore(results, largestPScore, cancerToCounts);
        builder.append("\t").append(rankScore);
        // CancerType objects may be loaded from different hibernate session,
        // therefore abbreviations are used as keys.
        Map<String, AnalysisResult> typeToResult = new HashMap<>();
        results.forEach(result -> typeToResult.put(result.getCancerType().getAbbreviation(),
                                                   result));
        double pancanPScore = calculatePancanPScore(typeToResult, largestPScore, cancerToCounts);
        builder.append("\t").append(pancanPScore);
        cancerTypes.forEach(type -> {
            AnalysisResult result = typeToResult.get(type);
            builder.append("\t");
            if (result == null)
                builder.append("");
            else
                builder.append(result.getFdr());
        });
        return builder.toString();
    }
    
}
