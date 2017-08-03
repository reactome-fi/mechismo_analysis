/*
 * Created on Sep 6, 2016
 *
 */
package org.reactome.cancer.driver;

import org.apache.log4j.Logger;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.gk.model.GKInstance;
import org.gk.persistence.MySQLAdaptor;
import org.gk.util.StringUtils;
import org.junit.Test;
import org.reactome.cancer.MAFFileLoader;
import org.reactome.px.util.InteractionUtilities;
import org.reactome.r3.*;
import org.reactome.r3.Interactome3dAnalyzer.PDBUniProtMatch;
import org.reactome.r3.ProteinSequenceHandler.Sequence;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.MutationObservation;

import java.io.*;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * This class is used to handle protein 3D structures based on PDB using biojava APIs.
 *
 * @author gwu
 */
public class Interactome3dDriverAnalyzer {
    private FileUtility fu = new FileUtility();
    // Cahced some values for quick performance
    private Map<String, Sequence> accToSeq;
    private Map<String, Integer> geneToLength;
    // a flag
    private boolean needDetails = false;
    // A temp to hold the min-pvalue from checking interface
    private double minPValue;
    private final static Logger logger = Logger.getLogger(Interactome3dDriverAnalyzer.class);

    /**
     * Default constructor.
     */
    public Interactome3dDriverAnalyzer() {
    }

    /**
     * Load a map from reaction to minimum p-value as -log10.
     *
     * @return
     * @throws IOException
     */
    public Map<String, Double> loadReactionTo3dScore() throws IOException {
        String fileName = "results/checkAllHumanReactions_summary_052417.txt";
        Map<String, Double> reactionToScore = Files.lines(Paths.get(fileName))
                .skip(1) // Skip the first header line
                .map(line -> line.split("\t"))
                .collect(Collectors.toMap(tokens -> tokens[1],
                        tokens -> -Math.log10(Double.parseDouble(tokens[2])),
                        (v1, v2) -> Math.max(v1, v2))); // In case of two reactions having the same name
        System.out.println("Total reactions in reaction to 3d score: " + reactionToScore.size());
        return reactionToScore;
    }

    /**
     * This method is used to generate a network view for structures.
     *
     * @throws Exception
     */
    @Test
    public void generateFINetworkForPDBs() throws Exception {
        String dirName = "datasets/interactome3d/2016_06/reactions_fdr_05/interactions";
        Map<String, File> fiToFile = new Interactome3dAnalyzer().loadPPIToPDBFile(dirName,
                false);
        System.out.println("Total FIs downloaded from interactome3d: " + fiToFile.size());
//        for (String fi : fiToFile.keySet())
//            System.out.println(fi);
        // Convert FIs into genes
        Map<String, String> proteinToGene = loadProteinToGene();
        System.out.println("Total proteins: " + proteinToGene.keySet().size());
        Set<String> fisInGenes = new HashSet<String>(); // Use set to merge FIs having same genes
        for (String fi : fiToFile.keySet()) {
            String[] tokens = fi.split("\t");
            String gene1 = proteinToGene.get(tokens[0]);
            String gene2 = proteinToGene.get(tokens[1]);
            String fiInGene = InteractionUtilities.generateFIFromGene(gene1, gene2);
            fisInGenes.add(fiInGene);
        }
        System.out.println("Total FIs in genes: " + fisInGenes.size());
        for (String fi : fisInGenes)
            System.out.println(fi);
    }

    private Map<String, String> loadProteinToGene() throws Exception {
        CancerDriverReactomeAnalyzer reactomeAnalyzer = new CancerDriverReactomeAnalyzer();
        Map<Long, Set<String>> rxtIdsToFIsWithFeatures = reactomeAnalyzer.loadReactionIdToFIsWithFeatures();
        System.out.println("Total reaction ids: " + rxtIdsToFIsWithFeatures.size());
        Map<String, String> proteinToGene = new HashMap<String, String>();
        for (Set<String> fisWithFeature : rxtIdsToFIsWithFeatures.values()) {
            for (String fiWithFeature : fisWithFeature) {
                String[] tokens = fiWithFeature.split("\t");
                proteinToGene.put(tokens[0], tokens[2]);
                proteinToGene.put(tokens[1], tokens[3]);
            }
        }
        return proteinToGene;
    }

    /**
     * Extract a list of analysis results for a set of FIs
     *
     * @throws Exception
     */
    @Test
    public void checkReactionFIsAnalysisResults() throws Exception {
        String dirName = "results/DriverGenes/Drivers_0816/";
        //String fileName = dirName + "MutationAnalysisResultsForReactionsFDR_05_PPI_092916.txt";
        String fileName = dirName + "MutationAnalysisResultsForReactionsFDR_05_PPI_101016.txt";
        fileName = dirName + "MutationAnalysisResultsForReactionsFDR_05_101016.txt";

        fu.setInput(fileName);
        String line = null;
        Set<String> pdbFiles = new HashSet<String>();
        System.out.println("Edge\tSource\tTarget\tSource_Threshold\tTarget_Threshold");
        while ((line = fu.readLine()) != null) {
            if (!line.endsWith(".pdb"))
                continue;
            Map<String, Double> geneToPValue = parseGeneToPValue(fu);
            if (geneToPValue.size() < 2)
                continue;
            List<String> geneList = new ArrayList<String>(geneToPValue.keySet());
            String gene1 = geneList.get(0);
            String gene2 = geneList.get(1);
            if (gene1.compareTo(gene2) > 0) {
                gene1 = geneList.get(1);
                gene2 = geneList.get(0);
            }
            System.out.println(gene1 + " (interacts with) " + gene2 + "\t" +
                    geneToPValue.get(gene1) + "\t" +
                    geneToPValue.get(gene2) + "\t" +
                    (geneToPValue.get(gene1) <= 0.01 ? 1 : 0) + "\t" +
                    (geneToPValue.get(gene2) <= 0.01 ? 1 : 0));
        }
        System.out.println("Total analyed PDB files: " + pdbFiles.size());
    }

    private Map<String, Double> parseGeneToPValue(FileUtility fu) throws IOException {
        Map<String, Double> geneToPValue = new HashMap<String, Double>();
        String line = fu.readLine(); // Escape the header
        line = fu.readLine(); // Value in the first chain
        if (_parseGeneToPValue(line, geneToPValue)) {
            line = fu.readLine();
            line = fu.readLine();
        } else {
            line = fu.readLine(); // Should pick up next line
        }
        if (line.length() == 0)
            return geneToPValue;
        _parseGeneToPValue(line, geneToPValue);
        return geneToPValue;
    }

    private boolean _parseGeneToPValue(String line,
                                       Map<String, Double> geneToValue) {
        if (!line.contains("\t")) {
            String[] tokens = line.split(" ");
            geneToValue.put(tokens[0], 1.0d);
            return false;
        }
        String[] tokens = line.split("\t");
        String gene = tokens[1];
        Double pvalue = new Double(tokens[tokens.length - 2]);
        Double aaPValue = new Double(tokens[tokens.length - 1]);
        if (aaPValue < pvalue)
            pvalue = aaPValue;
        geneToValue.put(gene, pvalue);
        return true;
    }

    @Test
    public void extractPValuesFromAnalysisResults() throws Exception {
        String dirName = "results/DriverGenes/Drivers_0816/";
        String fileName = dirName + "MutationAnalysisResultsForReactionsFDR_05_PPI_120516.txt";

        String titleLine = "Gene\tLength\tContacts\tRatio\tMutations\tMutationsInContacts";
        boolean isInUse = false;

        fu.setInput(fileName);
        String line = null;
        final Map<String, Double> geneToPValue = new HashMap<String, Double>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith(titleLine)) {
                isInUse = true;
            } else if (isInUse) {
                if (line.trim().length() == 0) {
                    isInUse = false;
                    continue;
                }
                String[] tokens = line.split("\t");
                Double pvalue = new Double(tokens[tokens.length - 2]);
                Double prePValue = geneToPValue.get(tokens[0]);
                if (prePValue == null || prePValue > pvalue)
                    geneToPValue.put(tokens[0], pvalue);
                pvalue = new Double(tokens[tokens.length - 1]);
                prePValue = geneToPValue.get(tokens[0]);
                if (prePValue > pvalue)
                    geneToPValue.put(tokens[0], pvalue);
            }
        }
        List<String> geneList = new ArrayList<String>(geneToPValue.keySet());
        Collections.sort(geneList, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Double pvalue1 = geneToPValue.get(gene1);
                Double pvalue2 = geneToPValue.get(gene2);
                return pvalue1.compareTo(pvalue2);
            }
        });

        Set<String> driverGenes = new CancerDriverAnalyzer().getDriverGenes(null);
        System.out.println("Total cancer driver genes: " + driverGenes.size());
        System.out.println();

        System.out.println("Gene\tP-Value\tCancerGenes");
        for (String gene : geneList)
            System.out.println(gene + "\t" + geneToPValue.get(gene) + "\t" +
                    (driverGenes.contains(gene) ? "1" : "0"));
    }

    /**
     * A method to analyze the results generated from method checkInteractionInAllReactions.
     *
     * @throws Exception
     */
    @Test
    public void checkAllReactionsAnalysisResults() throws Exception {
        String dirName = "results/DriverGenes/Drivers_0816/";
        //String fileName = dirName + "MutationAnalysisResultsForReactionsFDR_05_PPI_092916.txt";
        String fileName = dirName + "MutationAnalysisResultsForReactionsFDR_05_PPI_101016.txt";
        fileName = dirName + "MutationAnalysisResultsForReactionsFDR_05_101016.txt";
        fileName = dirName + "MutationAnalysisResultsForReactionsFDR_05_PPI_120516.txt";

        boolean isInSummary = false;
        fu.setInput(fileName);
        String line = null;
        // Reactions to be escaped
        Set<String> escapedReactions = new HashSet<String>();
        escapedReactions.add("5654402");
        escapedReactions.add("5654413");
        // Get all new genes
        Set<String> totalNewGenes = new HashSet<String>();
        Set<String> totalSigGenes = new HashSet<String>();
        Set<String> totalCheckedDriverGenes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("Summary")) {
                isInSummary = true;
                // Escape the title line
                line = fu.readLine();
            } else if (isInSummary) {
                if (line.length() == 0)
                    continue;
                String[] tokens = line.split("\t");
                if (escapedReactions.contains(tokens[0]))
                    continue;
                collectGenes(tokens[5], totalNewGenes);
                collectGenes(tokens[4], totalCheckedDriverGenes);
                collectGenes(tokens[3], totalSigGenes);
            }
        }
        fu.close();
        System.out.println("Total significant genes: " + totalSigGenes.size());
        System.out.println("Total checked known cancer genes: " + totalCheckedDriverGenes.size());
        System.out.println("Total new cancer genes: " + totalNewGenes.size());
        List<String> geneList = new ArrayList<String>(totalNewGenes);
        Collections.sort(geneList);
        for (String gene : geneList)
            System.out.print(gene + ", ");
        System.out.println();
        // Get known cancer genes that cannot be mapped
        List<String> missedDriverGenes = new ArrayList<String>(totalCheckedDriverGenes);
        missedDriverGenes.removeAll(totalSigGenes);
        Collections.sort(missedDriverGenes);
        System.out.println("Total missed driver genes: " + missedDriverGenes.size());
        for (String gene : missedDriverGenes)
            System.out.print(gene + ", ");
        System.out.println();

        // Generate output for Cytoscape
        Set<String> allGenes = new HashSet<String>();
        allGenes.addAll(totalSigGenes);
        allGenes.addAll(totalCheckedDriverGenes);
        allGenes.addAll(totalNewGenes);
        System.out.println("\n\nGene\tSignificant\tDriver\tScore");
        for (String gene : allGenes) {
            int driverScore = totalCheckedDriverGenes.contains(gene) ? 1 : 0;
            int sigScore = totalSigGenes.contains(gene) ? 2 : 0;
            int totalScore = driverScore + sigScore;
            System.out.println(gene + "\t" +
                    totalSigGenes.contains(gene) + "\t" +
                    totalCheckedDriverGenes.contains(gene) + "\t" +
                    totalScore);
        }
    }

    private void collectGenes(String txt,
                              Set<String> geneSet) {
        txt = txt.substring(1, txt.length() - 1); // Remove [ and ].
        if (txt.length() == 0)
            return;
        String[] genes = txt.split(", ");
        for (String gene : genes)
            geneSet.add(gene);
    }

    @Test
    public void checkInteractionInAllReactions() throws Exception {
        // A flag to exclude FIs having no physical interaction evidences
        boolean usePhysicalInteractionsOnly = false;
        CancerDriverReactomeAnalyzer reactomeAnalyzer = new CancerDriverReactomeAnalyzer();
        Map<Long, Set<String>> rxtIdsToFIsWithFeatures = reactomeAnalyzer.loadReactionIdToFIsWithFeatures();
        System.out.println("Total reaction ids: " + rxtIdsToFIsWithFeatures.size());
        // Want to get the total FIs
        Set<String> totalFIs = new HashSet<String>();
        for (Set<String> fis : rxtIdsToFIsWithFeatures.values())
            totalFIs.addAll(fis);
        System.out.println("Total FIs: " + totalFIs.size());

        String dirName = "datasets/interactome3d/2016_06/reactions_fdr_05/interactions";
        Map<String, File> fiToFile = new Interactome3dAnalyzer().loadPPIToPDBFile(dirName,
                false);
        System.out.println("Total FIs downloaded from interactome3d: " + fiToFile.size());

        Set<String> driverGenes = new CancerDriverAnalyzer().getDriverGenes(null);
        CosmicAnalyzer cosmicAnalyzer = new CosmicAnalyzer();
        Set<String> matchedGenes = getMatchedGenes(rxtIdsToFIsWithFeatures, fiToFile);
        Map<String, List<MutationObservation>> geneToCosmicEntries = cosmicAnalyzer.loadMutations(matchedGenes);

        MySQLAdaptor dba = reactomeAnalyzer.getDBA();
        Set<String> matchedFIs = new HashSet<String>();
        Set<String> totalCheckedFIs = new HashSet<String>();
        Set<String> genes = new HashSet<String>();
        Map<String, String> proteinToGene = new HashMap<String, String>();
        Set<File> pdbFiles = new HashSet<File>();
        StringBuilder summary = new StringBuilder();
        summary.append("DB_ID\tName\tTotalCheckedGenes\tSignificantGenes\tDriverGenes\tNewGenes\n");
        for (Long rxtId : rxtIdsToFIsWithFeatures.keySet()) {

//            // For debug
//            needDetails = true;
////            if (rxtId != 5675431L) // PP2A dephosphorylates RAF1
//            if (rxtId != 2671868L) // Phosphorylated LEPR Binds STAT3  
//                continue;

            genes.clear();
            proteinToGene.clear();
            pdbFiles.clear();
            GKInstance rxt = dba.fetchInstance(rxtId);
            if (rxt == null)
                throw new IllegalStateException(rxtId + " cannot be found!");
            Set<String> fisWithFeatures = rxtIdsToFIsWithFeatures.get(rxtId);
            int matched = 0;
            for (String fiWithFeatures : fisWithFeatures) {
                String[] tokens = fiWithFeatures.split("\t");
                // Check if PPI should be checked
                if (usePhysicalInteractionsOnly) {
                    boolean ppi = new Boolean(tokens[4]);
                    if (!ppi)
                        continue;
                }
                String fi = InteractionUtilities.generateFIFromGene(tokens[0], tokens[1]);
                totalCheckedFIs.add(fi);
                File structureFile = fiToFile.get(fi);
                if (structureFile == null)
                    continue;
                matched++;
                matchedFIs.add(fi);
                genes.add(tokens[2]);
                genes.add(tokens[3]);
                proteinToGene.put(tokens[0], tokens[2]);
                proteinToGene.put(tokens[1], tokens[3]);
                pdbFiles.add(structureFile);
            }
            double ratio = (double) matched / fisWithFeatures.size();
            System.out.println("\n#######################################");
            System.out.println(rxtId + "\t" +
                    rxt.getDisplayName() + "\t" +
                    fisWithFeatures.size() + "\t" +
                    matched + "\t" +
                    String.format("%.4f", ratio));
            if (genes.size() == 0) {
                System.out.println();
                continue;
            }
            // Print out cancer genes
            Set<String> cancerGenes = InteractionUtilities.getShared(driverGenes, genes);
            System.out.println("Driver genes: " + cancerGenes);
            Set<String> sigGenes = checkInterfaces(proteinToGene.keySet(),
                    proteinToGene,
                    geneToCosmicEntries,
                    pdbFiles);
            System.out.println();
            System.out.println("Significant genes: " + sigGenes);
            System.out.println("Driver genes: " + cancerGenes);
            Set<String> newSigGenes = new HashSet<String>(sigGenes);
            newSigGenes.removeAll(cancerGenes);
            System.out.println("Significant genes but not driver: " + newSigGenes);
            summary.append(rxtId).append("\t").append(rxt.getDisplayName()).append("\t");
            summary.append(genes.size()).append("\t");
            summary.append(sigGenes).append("\t").append(cancerGenes).append("\t").append(newSigGenes);
            summary.append("\n");
        }
        System.out.println("Total checked FIs: " + totalCheckedFIs.size());
        System.out.println("Total FIs having structures: " + matchedFIs.size());
        System.out.println("\nSummary:");
        System.out.println(summary.toString());
    }

    private Set<String> getMatchedGenes(Map<Long, Set<String>> rxtIdsToFIsWithFeatures,
                                        Map<String, File> fiToFile) {
        Set<String> genes = new HashSet<String>();
        for (Long rxtId : rxtIdsToFIsWithFeatures.keySet()) {
            Set<String> fisWithFeatures = rxtIdsToFIsWithFeatures.get(rxtId);
            for (String fiWithFeatures : fisWithFeatures) {
                String[] tokens = fiWithFeatures.split("\t");
                String fi = InteractionUtilities.generateFIFromGene(tokens[0], tokens[1]);
                File structureFile = fiToFile.get(fi);
                if (structureFile == null)
                    continue;
                genes.add(tokens[2]);
                genes.add(tokens[3]);
            }
        }
        return genes;
    }

    /**
     * This method implements a gene-central analysis by first loading mutation profiles annotated
     * in COSMIC or other places, find Reactome reactions having annotated genes involved, expand
     * reactions into a set of FIs, and check enrichment of mutations for FIs based on interactome3d
     * structures.
     *
     * @throws Exception
     */
    @Test
    public void checkAllHumanGenes() throws Exception {
        CancerDriverAnalyzer cancerDriverAnalyzer = new CancerDriverAnalyzer();
        Set<String> knownCancerGenes = cancerDriverAnalyzer.getDriverGenes(null);
        logger.info("Total known cancer driver genes: " + knownCancerGenes.size());
    }

    /**
     * This method overloads {@link #checkAllHumanReactions(CancerDriverReactomeAnalyzer)
     * checkAllHumanReactions} method.
     * for testing
     *
     * @throws Exception
     */
    @Test
    public void checkAllHumanReactions() throws Exception {
        checkAllHumanReactions(new CancerDriverReactomeAnalyzer());
    }

    /**
     * This method takes a reaction directly from the Reactome database, expands it into a set
     * of interactions in genes, loads structures using interactome3d data, and then performs mutation
     * interface analysis based on the cosmic annotation. The actual running is performed for all
     * human reactions in the database.
     *
     * @throws Exception
     */
    public void checkAllHumanReactions(CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer) throws Exception {
        // To control output
        boolean needReactionOutput = false;

        // Load all the non-disease reactions from Reactome
        // dbID + string description + other stuff
        List<GKInstance> reactions = cancerDriverReactomeAnalyzer.loadHumanReactions();
        System.out.println("Total reactions: " + reactions.size());

        // Map all interactions to a pdb file path
        // "<uniprot ID>\t<uniprot ID>" -> "path/to/complex.pdb"
        Interactome3dAnalyzer interactomeAnalyser = new Interactome3dAnalyzer();
        String interactomeDirName = "datasets/interactome3d/2016_06/prebuilt/representative/";
        Map<String, File> fiToPDB = interactomeAnalyser.loadPPIToPDBFile(interactomeDirName,
                false);

        // Store interactions in a set
        // "<uniprot ID>\t<uniprot ID>"
        ReactomeAnalyzer reactomeDataAnalyzer = new ReactomeAnalyzer();
        boolean hasStructure = false;
        int total = 0;
        Set<String> totalFIs = new HashSet<>();
        if (needReactionOutput)
            System.out.println("DB_ID\tName\tTotal_FIs\tHasStrucute");
        for (GKInstance reaction : reactions) {
            Set<String> fis = reactomeDataAnalyzer.generateTentativePPIsForReaction(reaction,
                    false);
            if (fis.size() == 0)
                continue;
            // Check if there is a structure available from interactome3d
            hasStructure = false;
            for (String fi : fis) {
                if (fiToPDB.containsKey(fi)) {
                    hasStructure = true;
                    total++;
                    break;
                }
            }
            if (needReactionOutput)
                System.out.println(reaction.getDBID() + "\t" +
                        reaction.getDisplayName() + "\t" +
                        fis.size() + "\t" +
                        hasStructure);
            totalFIs.addAll(fis);
        }
        System.out.println("Total reactions having FI structure: " + total);
        System.out.println("Total FIs: " + totalFIs.size());

        // Map uniprot ID's to gene symbols (COSMIC uses gene symbols)
        // "<uniprot ID>" -> "<gene symbol>"
        UniProtAnalyzer uniprotAnalyzer = new UniProtAnalyzer();
        Map<String, String> uniprotIdToGene = uniprotAnalyzer.getUniProtAccessionToGeneName();

        // Extract and store unique proteins from interactions in a set
        // "<uniprot ID>"
        Set<String> totalProteins = InteractionUtilities.grepIDsFromInteractions(totalFIs);
        System.out.println("Total proteins: " + totalProteins.size());

        // Lookup and store gene symbols for proteins in a set
        // "<gene symbol>"
        Set<String> totalGenes = new HashSet<>();
        for (String protein : totalProteins) {
            String gene = uniprotIdToGene.get(protein);
            if (gene == null) {
                //TODO: figure out which proteins don't map to genes and why
//                System.err.println("Cannot find gene for " + protein);
                continue;
            } else
                totalGenes.add(gene);
        }
        System.out.println("Total genes: " + totalGenes.size());

        // Map gene symbols to COSMIC mutation profiles
        // "<gene symbol>" -> [<MutationObservation>]
        CosmicAnalyzer cosmicAnalyzer = new CosmicAnalyzer();
        Map<String, List<MutationObservation>> geneToCosmicEntries = cosmicAnalyzer.loadMutations(totalGenes);
        System.out.println("Total genes in cosmic: " + geneToCosmicEntries.size());

        // Check interactions one by one
        Map<GKInstance, Double> reactionToMinPValue = new HashMap<>();
        Map<GKInstance, Integer> reactionToFINumber = new HashMap<>();
        Map<GKInstance, Set<String>> reactionToSigGenes = new HashMap<>();
        for (GKInstance reaction : reactions) {
            // Store interactions in a set
            // "<uniprot ID>\t<uniprot ID>"
            Set<String> fis = reactomeDataAnalyzer.generateTentativePPIsForReaction(reaction,
                    false);
            if (fis.size() == 0)
                continue;

            // Store pdb files for interactions in a set
            // <File>
            Set<File> pdbFiles = new HashSet<>();
            for (String fi : fis) {
                File file = fiToPDB.get(fi);
                if (file == null)
                    continue;
                pdbFiles.add(file);
            }
            if (pdbFiles.size() == 0)
                continue;

            // Extract and store unique proteins from interactions in a set
            // "<uniprot ID>"            
            System.out.println("\nChecking reaction " + reaction);
            Set<String> accessions = InteractionUtilities.grepIDsFromInteractions(fis);

            // Calculate significant genes?
            // this should be more clearly documented/described
            Set<String> sigGenes = checkInterfaces(accessions,
                    uniprotIdToGene,
                    geneToCosmicEntries,
                    pdbFiles);
            System.out.println(reaction.getDBID() + "\t" +
                    reaction.getDisplayName() + "\t" +
                    sigGenes.size() + "\t" +
                    sigGenes);

            // The whole reaction is being used as a key?
            // We should use the reaction ID for this..
            reactionToMinPValue.put(reaction, minPValue);
            reactionToFINumber.put(reaction, fis.size());
            reactionToSigGenes.put(reaction, sigGenes);
        }

        System.out.println("\nReactionToMinPValue:");
        System.out.println("DB_ID\tName\tMin_PValue\tFINumber\tSigGenes");

        for (GKInstance reaction : reactionToMinPValue.keySet())
            System.out.println(reaction.getDBID() + "\t" +
                    reaction.getDisplayName() + "\t" +
                    reactionToMinPValue.get(reaction) + "\t" +
                    reactionToFINumber.get(reaction) + "\t" +
                    StringUtils.join(", ", new ArrayList<String>(reactionToSigGenes.get(reaction))));
    }

    /**
     * This method overloads {@link #findInteractionsWithMutatedInterfaces(CancerDriverReactomeAnalyzer, String, String, String, String)
     * findInteractionsWithMutatedInterfaces} method for actual running.
     *
     * @throws Exception
     */
    @Test
    public void findInteractionsWithMutatedInterfaces() throws Exception {
        findInteractionsWithMutatedInterfaces(new CancerDriverReactomeAnalyzer(),
                null,
                null,
                null,
                null);
    }

    public void calculateMechismoInteractionCorrelationWithInterfaceMutationEnrichment(CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer,
                                                                                       String mechismoOutputFilePath,
                                                                                       String mechismoInteractionScorePattern,
                                                                                       String cancerTypeFilePath,
                                                                                       String cancerTypeInteractionPattern,
                                                                                       String mechismoOutFilePath,
                                                                                       String mafDirectoryPath,
                                                                                       String mafFileNamePattern,
                                                                                       String pdbDirectoryPath,
                                                                                       String reactomeOutFilePath,
                                                                                       String intersectionOutFilePath) throws Exception {

        //TODO: manual check -- write map to .csv and compare with interface enrichment (below)

        Map<String, InteractionMutationProfile> interfaceMutationRatios = findInteractionsWithMutatedInterfaces(cancerDriverReactomeAnalyzer,
                mafDirectoryPath,
                mafFileNamePattern,
                pdbDirectoryPath,
                reactomeOutFilePath);

        Map<String, Double> mechismoInteractionsInReactomeFIs = findMechismoInteractionsInReactome(cancerDriverReactomeAnalyzer,
                mechismoOutputFilePath,
                mechismoInteractionScorePattern,
                cancerTypeFilePath,
                cancerTypeInteractionPattern,
                mechismoOutFilePath);

        FileOutputStream fop = null;
        try {
            File file = new File(intersectionOutFilePath);
            file.createNewFile();
            fop = new FileOutputStream(intersectionOutFilePath);

            fop.write(("PPI," +
                    "Gene 1 Name, Gene1 Interface Length,Gene1 Length,Gene1 Interface Mutation Count,Gene1 Mutation Count," +
                    "Gene 2 Name, Gene2 Interface Length,Gene2 Length,Gene2 Interface Mutation Count,Gene2 Mutation Count," +
                    "P Value,FDR," +
                    "Mechismo Score\n").getBytes());

            for (String interactionWithInterface : interfaceMutationRatios.keySet()) {
                String[] interactionArray = interactionWithInterface.split("\t");
                String q1 = String.format("%s\t%s", interactionArray[0], interactionArray[1]);
                String q2 = String.format("%s\t%s", interactionArray[1], interactionArray[0]);
                double mechismoScore;
                if (mechismoInteractionsInReactomeFIs.containsKey(q1)) {
                    mechismoScore = mechismoInteractionsInReactomeFIs.get(q1);
                } else if (mechismoInteractionsInReactomeFIs.containsKey(q2)) {
                    mechismoScore = mechismoInteractionsInReactomeFIs.get(q2);
                } else {
                    logger.warn(String.format("Mechismo Interactions in ReactomeFIs missing: %s",
                            interactionWithInterface));
                    continue;
                }

                fop.write(String.format("%s,%s,%s,%f,%f,%f\n",
                        interactionWithInterface,
                        interfaceMutationRatios.get(interactionWithInterface).gene1MutationProfile.toString(),
                        interfaceMutationRatios.get(interactionWithInterface).gene2MutationProfile.toString(),
                        interfaceMutationRatios.get(interactionWithInterface).getPValue(),
                        interfaceMutationRatios.get(interactionWithInterface).getFdr(),
                        mechismoScore).getBytes());
            }
        } catch (IOException ioe) {
            logger.error(String.format("%s: %s",
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        } finally {
            if (fop != null) {
                fop.flush();
                fop.close();
            }
        }

        //TODO: correlate interface enrichment (below) with mechismo score
        //TODO: Heatmaps for interactions (Reactome vs Mechismo) per patient
        //TODO: Heatmaps for reactions per patient
        //TODO: Heatmaps for pathways per patient
        //TODO: Rerun analysis for all cancer types (along with interface enrichment comparison)
    }

    public void compareKnownDrivers(CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer,
    String pdbDirectoryPath,
    String outputDir,
    String firehoseCancerTypesDir,
    String knownDriverPath) throws Exception {
        //------------------
        //load known drivers
        //------------------
        FileInputStream fis = new FileInputStream(knownDriverPath);
        try{
            FileInputStream fileInputStream = new FileInputStream(knownDriverPath);

        }catch(IOException ioe){

        }finally {

        }
        //---------------------------------
        //report reactions of known drivers
        //---------------------------------

        //Reactome FI's
        // Load all the non-disease reactions from Reactome
        // dbID + string description + other stuff
        List<GKInstance> reactions = cancerDriverReactomeAnalyzer.loadHumanReactions();
        System.out.println(String.format("Total reactions: %d", reactions.size()));

        // Store interactions in a set
        // "<uniprot ID>\t<uniprot ID>"
        Set<String> totalFIs = extractInteractionsFromReactions(reactions);
        System.out.println(String.format("Total FIs: %d", totalFIs.size()));

        // Map uniprot ID's to gene symbols (COSMIC uses gene symbols)
        // "<uniprot ID>" -> "<gene symbol>"
        Map<String, String> uniprotIdToGene = new UniProtAnalyzer().getUniProtAccessionToGeneName();

        // Extract and store unique proteins from interactions in a set
        // "<uniprot ID>"
        Set<String> totalProteins = InteractionUtilities.grepIDsFromInteractions(totalFIs);
        System.out.println(String.format("Total proteins: %d", totalProteins.size()));

        // Keep only uniprot ID's from interactions in uniprotIdToGene map
        uniprotIdToGene.keySet().retainAll(totalProteins);

        //Keep only Reactome FI's containing >= 1 gene mutated in MAF
        // Lookup and store gene symbols for proteins in a set
        // "<gene symbol>"
        Set<String> totalGenes = keepGenesWithMappedProteins(totalProteins, uniprotIdToGene);
        System.out.println(String.format("Total genes: %d", totalGenes.size()));

        //Map Interactions to sets of Reactions
        Map<String, Set<GKInstance>> interactionReactionSetMap = generateInteractionReactionSetMap(uniprotIdToGene, reactions);

        ///////some function like this.... probs have to write another one/////////////
        //findSampleReactionsWithMutatedInteractionInterfaces(
        //        cancerTypeSampleInterfaceMutationRatios,
        //        interactionReactionSetMap,
        //        reactionOutFilePath);

        //-----------------------------------------
        //report mutation hotspots in known drivers
        //-----------------------------------------

        //-----------------------------------------
        //report known driver PPI's with structures
        //-----------------------------------------

        //------------------------------------------------------------------
        //report known driver PPI's interface mutation ratios & significance
        //------------------------------------------------------------------

        //------------------------------------------------
        //report reactions affected by interface mutations
        //------------------------------------------------

        //---------------------------------------------
        //report pathways containing affected reactions
        //---------------------------------------------

        //------------------------------------------------
        //report pan-cancer & cancer type specific results
        //------------------------------------------------
    }

    public void prepareHeatmapData(CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer,
                                   String pdbDirectoryPath,
                                   String outputDir,
                                   String firehoseCancerTypesDir,
                                   String mechismoCancerTypesDir) throws Exception {

        //Reactome FI's
        // Load all the non-disease reactions from Reactome
        // dbID + string description + other stuff
        List<GKInstance> reactions = cancerDriverReactomeAnalyzer.loadHumanReactions();
        System.out.println(String.format("Total reactions: %d", reactions.size()));

        // Store interactions in a set
        // "<uniprot ID>\t<uniprot ID>"
        Set<String> totalFIs = extractInteractionsFromReactions(reactions);
        System.out.println(String.format("Total FIs: %d", totalFIs.size()));

        // Map uniprot ID's to gene symbols (COSMIC uses gene symbols)
        // "<uniprot ID>" -> "<gene symbol>"
        Map<String, String> uniprotIdToGene = new UniProtAnalyzer().getUniProtAccessionToGeneName();

        // Extract and store unique proteins from interactions in a set
        // "<uniprot ID>"
        Set<String> totalProteins = InteractionUtilities.grepIDsFromInteractions(totalFIs);
        System.out.println(String.format("Total proteins: %d", totalProteins.size()));

        // Keep only uniprot ID's from interactions in uniprotIdToGene map
        uniprotIdToGene.keySet().retainAll(totalProteins);

        //Keep only Reactome FI's containing >= 1 gene mutated in MAF
        // Lookup and store gene symbols for proteins in a set
        // "<gene symbol>"
        Set<String> totalGenes = keepGenesWithMappedProteins(totalProteins, uniprotIdToGene);
        System.out.println(String.format("Total genes: %d", totalGenes.size()));

        // Map all pdb files matching uniprot ID's to gene names;
        // "<gene symbol>\t<gene symbol>" -> "path/to/complex.pdb"
        Map<String, File> interaction2PDBMap = getInteraction2PDBMap(pdbDirectoryPath, uniprotIdToGene);
        System.out.println(String.format("Total PDB's: %d", interaction2PDBMap.size()));

        //Extract Interfaces
        // "<gene symbol\tgene symbol>" -> HashMap<"<gene symbol>", Set<interface coordinate>>
        Map<String, Map<String, ProteinChainSummary>> interactionsWithInterfaces = extractInterfaces(interaction2PDBMap, uniprotIdToGene);

        //Map Interactions to sets of Reactions
        Map<String, Set<GKInstance>> interactionReactionSetMap = generateInteractionReactionSetMap(uniprotIdToGene, reactions);

        //For each cancer type in firehose
        File firehoseCancerTypes = new File(firehoseCancerTypesDir);
        String[] firehoseCancerTypeDirs =
                firehoseCancerTypes.list((current, name) -> new File(current, name).isDirectory());
        Map<String, Map<String, InteractionMutationProfile>> cancerTypeSampleInterfaceMutationRatios;
        Map<String, InteractionMutationProfile> cancerTypeTotalInterfaceMutationRatios;
        // The MAF filenames look like:
        // TCGA-AY-4070-01.hg19.oncotator.hugo_entrez_remapped.maf.txt
        // TCGA-AY-4071-01.hg19.oncotator.hugo_entrez_remapped.maf.txt
        String mafFileNamePattern = "^.+\\.maf\\.txt$";
        double fdrThresh = 0.05;
        for (String firehoseCancerTypeDir : firehoseCancerTypeDirs) {

            /**
             * CANCER TYPE RESOLUTION
             */
            String mafDirPath = String.format("%s/%s",
                    firehoseCancerTypesDir,
                    firehoseCancerTypeDir);

            //MAF Mutations
            // "<gene symbol>" -> HashMap<"<AA Coord> -> <sample support>>
            Map<String, Map<Integer, Set<MutationObservation>>> allSamplesMutatedGeneMap = summarizeAllMafs(mafFileNamePattern, mafDirPath);

            String cancerTypeOutFilePath = String.format("%s/firehose_%s_total_interactions_with_mutated_interfaces.csv",
                    outputDir,
                    firehoseCancerTypeDir);
            cancerTypeTotalInterfaceMutationRatios = processInteractionsWithMutatedInterfaces(interactionsWithInterfaces,
                    allSamplesMutatedGeneMap,
                    cancerTypeOutFilePath);

            Set<String> significantInteractions = new HashSet<>();
            Map<String, Map<String, Double>> totalInteractionSignificance = new HashMap<>();
            for (String interaction : cancerTypeTotalInterfaceMutationRatios.keySet()) {
                if (cancerTypeTotalInterfaceMutationRatios.get(interaction).getFdr() < fdrThresh) {
                    significantInteractions.add(interaction);
                    Map<String, Double> interactionSignificance = new HashMap<>();
                    interactionSignificance.put("PValue",cancerTypeTotalInterfaceMutationRatios.get(interaction).getPValue());
                    interactionSignificance.put("FDR",cancerTypeTotalInterfaceMutationRatios.get(interaction).getFdr());
                    totalInteractionSignificance.put(interaction,interactionSignificance);
                }
            }

            //Keep only significant interactions for individual sample-level analysis
            Map<String, Map<String, ProteinChainSummary>> significantInteractionsWithInterfaces = new HashMap<>(interactionsWithInterfaces);
            significantInteractionsWithInterfaces.keySet().retainAll(significantInteractions);

            /**
             * SAMPLE RESOLUTION
             */

            if (significantInteractionsWithInterfaces.size() > 0) {
                //heatmap of cancer type with each sample using p-value as value

                //ppi
                String ppiOutFilePath = String.format("%s/firehose_%s_individual_interactions_with_mutated_interfaces.csv",
                        outputDir,
                        firehoseCancerTypeDir);
                cancerTypeSampleInterfaceMutationRatios =
                        findSamplesInteractionsWithMutatedInterfaces(totalGenes,
                                significantInteractionsWithInterfaces,
                                totalInteractionSignificance,
                                mafDirPath,
                                mafFileNamePattern,
                                ppiOutFilePath);
                //reaction
                String reactionOutFilePath = String.format("%s/firehose_%s_individual_reactions_with_mutated_interfaces.csv",
                        outputDir,
                        firehoseCancerTypeDir);
                findSampleReactionsWithMutatedInteractionInterfaces(
                        cancerTypeSampleInterfaceMutationRatios,
                        interactionReactionSetMap,
                        reactionOutFilePath);

                //pathway
                //List<GeneSetAnnotation> annotateGenesWithReactomePathways(Collection<String> genes...);
            } else {
                logger.warn(String.format("No significant interactions found in %s.", firehoseCancerTypeDir));
            }
        }
        //heatmap of firehose data with cancer type using p-value as value? Maybe do this with R using output above

        //For each cancer type in mechismo
        //heatmap of cancer type with each sample using p-value as value
        //heatmap of cancer type with each sample using mechismo score as value

        //heatmap of mechismo data with cancer type using p-value as value
        //heatmap of mechismo data with cancer type using mechismo score as value


    }

    public Map<String, Double> findMechismoInteractionsInReactome(CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer,
                                                                  String mechismoOutputFilePath,
                                                                  String mechismoInteractionScorePattern,
                                                                  String cancerTypeFilePath,
                                                                  String cancerTypeInteractionPattern,
                                                                  String outFilePath) throws Exception {
        //Reactome FI's
        // Load all the non-disease reactions from Reactome
        // dbID + string description + other stuff
        List<GKInstance> reactions = cancerDriverReactomeAnalyzer.loadHumanReactions();
        System.out.println(String.format("Total reactions: %d", reactions.size()));

        // Store Reactome interactions in a set
        // "<uniprot ID>\t<uniprot ID>"
        Set<String> reactomeFIs = extractInteractionsFromReactions(reactions);
        System.out.println(String.format("Reactome FIs: %d", reactomeFIs.size()));

        // Store Mechismo interactions in a map
        // "<uniprot ID>\t<uniprot ID> => <mechismo score>"
        Map<String, Double> mechismoInteractions = mechismoInteractionMap(mechismoOutputFilePath, mechismoInteractionScorePattern);
        System.out.println(String.format("Mechismo PPIs: %d", mechismoInteractions.size()));

        // Report Coverage
        Map<String, Double> mechismoInteractionsInReactomeFIs = findMechismoReactomeIntersection(reactomeFIs, mechismoInteractions);
        System.out.println(String.format("%d of %d (%.3f %%) PPIs in Mechismo covered by Reactome FIs",
                mechismoInteractionsInReactomeFIs.keySet().size(),
                mechismoInteractions.keySet().size(),
                100.0 * ((double) mechismoInteractionsInReactomeFIs.keySet().size() / (double) mechismoInteractions.keySet().size())));

        // Map uniprot ID's to gene symbols (COSMIC uses gene symbols)
        // "<uniprot ID>" -> "<gene symbol>"
        Map<String, String> uniprotIdToGene = new UniProtAnalyzer().getUniProtAccessionToGeneName();
        FileOutputStream fop = null;
        try {
            File file = new File(outFilePath);
            file.createNewFile();
            fop = new FileOutputStream(outFilePath);

            fop.write(("Interaction,Gene 1, Gene 2,Mechismo Score\n").getBytes());
            for (String interaction : mechismoInteractionsInReactomeFIs.keySet()) {
                String[] interactionArray = interaction.split("\t");
                String gene1 = uniprotIdToGene.get(interactionArray[0]);
                String gene2 = uniprotIdToGene.get(interactionArray[1]);
                fop.write((String.format("%s\t%s,%s,%s,%f\n",
                        gene1,
                        gene2,
                        gene1,
                        gene2,
                        mechismoInteractionsInReactomeFIs.get(interaction)).getBytes()));
            }
        } catch (IOException ioe) {
            logger.error(String.format("%s: %s",
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        } finally {
            if (fop != null) {
                fop.flush();
                fop.close();
            }
        }
        return mechismoInteractionsInReactomeFIs;
    }

    private Map<String, Double> findMechismoReactomeIntersection(Set<String> reactomeFIs, Map<String, Double> mechismoInteractions) {
        Set<String> sortedFIs = new HashSet<>();
        for (String fi : reactomeFIs) {
            String[] fiArray = fi.split("\t");
            Arrays.sort(fiArray);
            sortedFIs.add(String.format("%s\t%s", fiArray[0], fiArray[1]));
        }

        Map<String, Double> mechismoInteractionsInReactomeFIs = new HashMap<>();
        for (String interaction : mechismoInteractions.keySet()) {
            String[] interactionArray = interaction.split("\t");
            Arrays.sort(interactionArray);
            String sortedInteraction = String.format("%s\t%s", interactionArray[0], interactionArray[1]);
            if (sortedFIs.contains(sortedInteraction)) {
                mechismoInteractionsInReactomeFIs.put(sortedInteraction, mechismoInteractions.get(interaction));
            }
        }
        return mechismoInteractionsInReactomeFIs;
    }

    private Map<String, Double> mechismoInteractionMap(String mechismoOutputFilePath,
                                                       String mechismoInteractionScorePattern) throws IOException {
        Pattern mechIntScorePattern = Pattern.compile(mechismoInteractionScorePattern);
        String prot1 = "prot1";
        String prot2 = "prot2";
        String mech = "mech";

        Map<String, Double> mechismoInteractions = new HashMap<>();

        FileUtility fu = new FileUtility();
        fu.setInput(mechismoOutputFilePath);
        String line;
//        while ((line = fu.readLine()) != null) {
//            Matcher matcher = mechIntScorePattern.matcher(line);
//            if(matcher.find()) {
//                mechismoInteractions.put(String.format("%s\t%s",matcher.group(prot1),matcher.group(prot2)), Double.parseDouble(matcher.group(mech)));
//            }
//        }
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 20)
                continue;
            if (tokens[19].trim().length() == 0)
                continue;
            mechismoInteractions.put(tokens[1].trim() + "\t" + tokens[19].trim(),
                    new Double(tokens[17]));
        }
        fu.close();
        return mechismoInteractions;
    }

    public Map<String, InteractionMutationProfile> findInteractionsWithMutatedInterfaces(CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer,
                                                                                         String mafDirectoryPath,
                                                                                         String mafFileNamePattern,
                                                                                         String pdbDirectoryPath,
                                                                                         String outFilePath) throws Exception {
        //Reactome FI's
        // Load all the non-disease reactions from Reactome
        // dbID + string description + other stuff
        List<GKInstance> reactions = cancerDriverReactomeAnalyzer.loadHumanReactions();
        System.out.println(String.format("Total reactions: %d", reactions.size()));

        // Store interactions in a set
        // "<uniprot ID>\t<uniprot ID>"
        Set<String> totalFIs = extractInteractionsFromReactions(reactions);
        System.out.println(String.format("Total FIs: %d", totalFIs.size()));

        // Map uniprot ID's to gene symbols (COSMIC uses gene symbols)
        // "<uniprot ID>" -> "<gene symbol>"
        Map<String, String> uniprotIdToGene = new UniProtAnalyzer().getUniProtAccessionToGeneName();

        // Extract and store unique proteins from interactions in a set
        // "<uniprot ID>"
        Set<String> totalProteins = InteractionUtilities.grepIDsFromInteractions(totalFIs);
        System.out.println(String.format("Total proteins: %d", totalProteins.size()));

        // Keep only uniprot ID's from interactions in uniprotIdToGene map
        uniprotIdToGene.keySet().retainAll(totalProteins);

        //Keep only Reactome FI's containing >= 1 gene mutated in MAF
        // Lookup and store gene symbols for proteins in a set
        // "<gene symbol>"
        Set<String> totalGenes = keepGenesWithMappedProteins(totalProteins, uniprotIdToGene);
        System.out.println(String.format("Total genes: %d", totalGenes.size()));

        //MAF Mutations
        // "<gene symbol>" -> HashMap<"<AA Coord> -> <sample support>>
        Map<String, Map<Integer, Set<MutationObservation>>> allSamplesGeneMap = summarizeAllMafs(mafFileNamePattern, mafDirectoryPath);

        //TODO: ensure this works... caps? some genes have multiple names (Akt/PKB) etc.
        allSamplesGeneMap.keySet().retainAll(totalGenes);

        // Map all pdb files matching uniprot ID's to gene names;
        // "<gene symbol>\t<gene symbol>" -> "path/to/complex.pdb"
        Map<String, File> interaction2PDBMap = getInteraction2PDBMap(pdbDirectoryPath, uniprotIdToGene);
        System.out.println(String.format("Total PDB's: %d", interaction2PDBMap.size()));

        //Extract Interfaces
        // "<gene symbol\tgene symbol>" -> HashMap<"<gene symbol>", Set<interface coordinate>>
        Map<String, Map<String, ProteinChainSummary>> interactionsWithInterfaces = extractInterfaces(interaction2PDBMap, uniprotIdToGene);

        return (processInteractionsWithMutatedInterfaces(interactionsWithInterfaces,
                allSamplesGeneMap,
                outFilePath));
    }

    /**
     * This method...
     *
     * @throws Exception
     */
    public Map<String, InteractionMutationProfile> processInteractionsWithMutatedInterfaces(Map<String, Map<String, ProteinChainSummary>> interactionsWithInterfaces,
                                                                                            Map<String, Map<Integer, Set<MutationObservation>>> allSamplesMutatedGeneMap,
                                                                                            String outFilePath) throws Exception {
        //Reactome FI's Interface Mutation Ratio
        // "<gene symbol\tgene symbol>" -> HashMap<<"gene symbol">,<Mutation Ratio>>
        Map<String, InteractionMutationProfile> interfaceMutationRatios = new HashMap<>();

        for(String interactionWithInterface: interactionsWithInterfaces.keySet()) {
            String gene1Name = interactionWithInterface.split("\t")[0];
            String gene2Name = interactionWithInterface.split("\t")[1];

            GeneMutationProfile gene1MutationProfile = getGeneMutationProfile(allSamplesMutatedGeneMap, interactionsWithInterfaces, interactionWithInterface, gene1Name);
            GeneMutationProfile gene2MutationProfile = getGeneMutationProfile(allSamplesMutatedGeneMap, interactionsWithInterfaces, interactionWithInterface, gene2Name);

            if (gene1MutationProfile.isValid() && gene2MutationProfile.isValid()) {

                int interactionInterfaceMutationCount = gene1MutationProfile.interfaceMutationCount
                        + gene2MutationProfile.interfaceMutationCount;

                if (interactionInterfaceMutationCount > 0) {

                    double pValue = 1.0;
                    //report lowest p value
                    pValue = gene1MutationProfile.calculatePValue() < gene2MutationProfile.calculatePValue()
                            ? gene1MutationProfile.calculatePValue()
                            : gene2MutationProfile.calculatePValue();

                    interfaceMutationRatios.put(interactionWithInterface, new InteractionMutationProfile(
                            gene1MutationProfile,
                            gene2MutationProfile,
                            pValue,
                            1.0,
                            interactionInterfaceMutationCount > 0
                    ));

                } else {
                    //skip interactions without mutations
                }
            } else {
                logger.warn(String.format("Invalid gene mutation profile pair: %s - %s",
                        gene1MutationProfile.toString(),
                        gene2MutationProfile.toString()));
            }
        }
        List<Double> pValues = new ArrayList<>();
        for (String interaction : interfaceMutationRatios.keySet()) {
            pValues.add(interfaceMutationRatios.get(interaction).getPValue());
        }
        List<Double> fdrs = MathUtilities.calculateFDRWithBenjaminiHochberg(pValues);
        Map<Double, Double> pValueFDRMap = new HashMap<>();
        for (int i = 0; i < pValues.size(); i++) {
            pValueFDRMap.put(pValues.get(i), fdrs.get(i));
        }
        for (String interaction : interfaceMutationRatios.keySet()) {
            double pValue = interfaceMutationRatios.get(interaction).getPValue();
            interfaceMutationRatios.get(interaction).setFdr(pValueFDRMap.get(pValue));
        }

        FileOutputStream fop = null;
        try {
            File file = new File(outFilePath);
            file.createNewFile();
            fop = new FileOutputStream(outFilePath);
            fop.write(("Interaction," +
                    "Gene1 Name,Gene1 Interface Length,Gene1 Length,Gene1 Interface Mutation Count,Gene1 Mutation Count," +
                    "Gene2 Name,Gene2 Interface Length,Gene2 Length,Gene2 Interface Mutation Count,Gene2 Mutation Count," +
                    "P Value,FDR\n").getBytes());

            for (String interaction: interfaceMutationRatios.keySet()) {
                InteractionMutationProfile interactionMutationProfile = interfaceMutationRatios.get(interaction);
                        fop.write(String.format("%s,%s,%s,%f,%f\n",
                                interaction,
                                interactionMutationProfile.getGene1MutationProfile().toString(),
                                interactionMutationProfile.getGene2MutationProfile().toString(),
                                interactionMutationProfile.getPValue(),
                                interactionMutationProfile.getFdr()).getBytes());
            }
            logger.info(String.format("Output written to %s", outFilePath));
        } catch (IOException ioe) {
            logger.error(String.format("%s: %s",
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        } finally {
            if (fop != null) {
                fop.flush();
                fop.close();
            }
        }

        return interfaceMutationRatios;
    }

    /**
     * This method...
     *
     * @throws Exception
     */
    public Map<String, Map<String, InteractionMutationProfile>> findSamplesInteractionsWithMutatedInterfaces(Set<String> totalGenes,
                                                                                                             Map<String, Map<String, ProteinChainSummary>> interactionsWithInterfaces,
                                                                                                             Map<String, Map<String, Double>> totalInteractionSignificance,
                                                                                                             String mafDirectoryPath,
                                                                                                             String mafFileNamePattern,
                                                                                                             String outFilePath) throws Exception {
        //MAF Mutations
        // "<gene symbol>" -> HashMap<"<AA Coord> -> <sample support>>
        Map<String, Map<String, Map<Integer, Set<MutationObservation>>>> individualSamplesGeneMap = summarizeIndividualMafs(mafFileNamePattern, mafDirectoryPath);

        //TODO: ensure this works... caps? some genes have multiple names (Akt/PKB) etc.
        for (String individualSample : individualSamplesGeneMap.keySet()) {
            individualSamplesGeneMap.get(individualSample).keySet().retainAll(totalGenes);
        }


        //Reactome FI's Interface Mutation Ratio
        // "<gene symbol\tgene symbol>" -> HashMap<<"gene symbol">,<Mutation Ratio>>
        Map<String, Map<String, InteractionMutationProfile>> sampleInterfaceMutationRatios = new HashMap<>();


        FileOutputStream fop = null;
        try {
            File file = new File(outFilePath);
            file.createNewFile();
            fop = new FileOutputStream(outFilePath);

            fop.write((
                    "Sample,Interaction," +
                            "Gene1 Name,Gene1 Interface Length,Gene1 Length,Gene1 Interface Mutation Count,Gene1 Mutation Count," +
                            "Gene2 Name,Gene2 Interface Length,Gene2 Length,Gene2 Interface Mutation Count,Gene2 Mutation Count," +
                            "PValue,FDR,Interface Mutation Present\n").getBytes());

            for (String individualSample : individualSamplesGeneMap.keySet()) {
                sampleInterfaceMutationRatios.put(individualSample, new HashMap<>());
                for (String interactionWithInterface : interactionsWithInterfaces.keySet()) {

                    double pValue = totalInteractionSignificance.get(interactionWithInterface).get("PValue");
                    double fdr = totalInteractionSignificance.get(interactionWithInterface).get("FDR");

                    String gene1Name = interactionWithInterface.split("\t")[0];
                    String gene2Name = interactionWithInterface.split("\t")[1];

                    GeneMutationProfile gene1MutationProfile = getGeneMutationProfile(individualSamplesGeneMap.get(individualSample), interactionsWithInterfaces, interactionWithInterface, gene1Name);
                    GeneMutationProfile gene2MutationProfile = getGeneMutationProfile(individualSamplesGeneMap.get(individualSample), interactionsWithInterfaces, interactionWithInterface, gene2Name);

                    if (gene1MutationProfile.isValid() && gene2MutationProfile.isValid()) {

                        int interactionInterfaceMutationCount = gene1MutationProfile.interfaceMutationCount
                                + gene2MutationProfile.interfaceMutationCount;

                        //Binary:
                        //1 = interaction has interface mutations
                        //0 = interaction does not have interface mutations
                        boolean hasInterfaceMuts = interactionInterfaceMutationCount > 0 ? true : false;

                        sampleInterfaceMutationRatios.get(individualSample).put(interactionWithInterface, new InteractionMutationProfile(
                                gene1MutationProfile,
                                gene2MutationProfile,
                                pValue,
                                fdr,
                                hasInterfaceMuts
                        ));

                        fop.write(String.format("%s,%s,%s,%s,%f,%f,%b\n",
                                individualSample,
                                interactionWithInterface,
                                gene1MutationProfile.toString(),
                                gene2MutationProfile.toString(),
                                pValue,
                                fdr,
                                hasInterfaceMuts).getBytes());
                    } else {
                        logger.warn(String.format("Invalid gene mutation profile pair: %s - %s",
                                gene1MutationProfile.toString(),
                                gene2MutationProfile.toString()));
                    }
                }
            }
            logger.info(String.format("Output written to %s", outFilePath));
        } catch (IOException ioe) {
            logger.error(String.format("%s: %s",
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        } finally {
            if (fop != null) {
                fop.flush();
                fop.close();
            }
        }

        return sampleInterfaceMutationRatios;
    }

    /***
     *
     * @param cancerTypeSampleInterfaceMutationRatios "sample id" -> Map<"gene\tgene" -> InteractionMutationProfile>
     * @param interactionReactionSetMap "gene\tgene" -> Set<reaction>
     * @param outFilePath
     */
    private void findSampleReactionsWithMutatedInteractionInterfaces(Map<String, Map<String, InteractionMutationProfile>> cancerTypeSampleInterfaceMutationRatios,
                                                                     Map<String, Set<GKInstance>> interactionReactionSetMap,
                                                                     String outFilePath) throws Exception {
        //map of reactions to interactions found to have
        //significantly mutated interface in at least one sample
        Map<GKInstance, Set<String>> reactionToInteractionMap = new HashMap<>();
        for (String sample : cancerTypeSampleInterfaceMutationRatios.keySet()) {
            for (String interaction : cancerTypeSampleInterfaceMutationRatios.get(sample).keySet()) {
                Set<GKInstance> reactions = interactionReactionSetMap.get(interaction);
                if (reactions != null && !reactions.isEmpty()) {
                    Iterator<GKInstance> gkInstanceIterator = reactions.iterator();
                    while (gkInstanceIterator.hasNext()) {
                        GKInstance reaction = gkInstanceIterator.next();
                        Set<String> interactions;
                        if (reactionToInteractionMap.keySet().contains(reaction)) {
                            interactions = reactionToInteractionMap.get(reaction);
                        } else {
                            interactions = new HashSet<>();
                        }
                        interactions.add(interaction);
                        reactionToInteractionMap.put(reaction, interactions);
                    }
                } else {
                    //The interaction isn't found in our reaction set
                    //logger.warn(String.format("No reactions mapped to %s...", interaction));
                }
            }
        }
        System.out.println(String.format("Reactions with interactions: %d", reactionToInteractionMap.keySet().size()));

        //write reactions and maximum interaction interface mutation significance for samples
        FileOutputStream fop = null;
        try {
            File file = new File(outFilePath);
            file.createNewFile();
            fop = new FileOutputStream(outFilePath);

            fop.write((
                    "Sample~Reaction ID~Reaction Description~AffectedInteractions~InterfaceMutationCount\n").getBytes());
            for (String sample : cancerTypeSampleInterfaceMutationRatios.keySet()) {
                for (GKInstance reaction : reactionToInteractionMap.keySet()) {
                    String affectedInteractions = "";
                    int interactionInterfaceMutationCount = 0;

                    Iterator<String> interactionIterator = reactionToInteractionMap.get(reaction).iterator();
                    while (interactionIterator.hasNext()) {
                        String interaction = interactionIterator.next();
                        if (cancerTypeSampleInterfaceMutationRatios.get(sample).containsKey(interaction)) {
                            if (cancerTypeSampleInterfaceMutationRatios.get(sample).get(interaction).getHasInterfaceMutations())
                                ;
                            if (interactionInterfaceMutationCount == 0) {
                                affectedInteractions = interaction;
                            } else {
                                affectedInteractions = String.format("%s,%s", affectedInteractions, interaction);
                            }
                            interactionInterfaceMutationCount += 1;
                        }
                    }
                    fop.write(String.format("%s~%s~%s~%s~%d\n",
                            sample,
                            reaction.getDBID(),
                            reaction.getDisplayName(),
                            affectedInteractions,
                            interactionInterfaceMutationCount).getBytes());
                }
            }
        } catch (IOException ioe) {
            logger.error(String.format("%s: %s",
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        } finally {
            if (fop != null) {
                fop.flush();
                fop.close();
            }
        }
    }

    private class ProteinChainSummary {
        private Set<Integer> interfaceCoordinates;
        private int chainLength;

        ProteinChainSummary(Set<Integer> interfaceCoordinates,
                            int chainLength) {
            this.interfaceCoordinates = interfaceCoordinates;
            this.chainLength = chainLength;
            if (this.chainLength <= 0) {
                logger.error("Invalid chain length detected...");
            }
        }

        public Set<Integer> getInterfaceCoordinates() {
            return this.interfaceCoordinates;
        }

        public int getChainLength() {
            return this.chainLength;
        }
    }

    private class InteractionMutationProfile {
        private GeneMutationProfile gene1MutationProfile;
        private GeneMutationProfile gene2MutationProfile;
        private double pValue;
        private boolean hasIntefaceMutations;
        private double fdr;

        InteractionMutationProfile(GeneMutationProfile gene1MutationProfile,
                                   GeneMutationProfile gene2MutationProfile,
                                   double pValue,
                                   double fdr,
                                   boolean hasInterfaceMutations) {
            this.gene1MutationProfile = gene1MutationProfile;
            this.gene2MutationProfile = gene2MutationProfile;
            this.pValue = pValue;
            this.fdr = fdr;
            this.hasIntefaceMutations = hasInterfaceMutations;
        }

        public GeneMutationProfile getGene1MutationProfile() {
            return this.gene1MutationProfile;
        }

        public GeneMutationProfile getGene2MutationProfile() {
            return this.gene2MutationProfile;
        }

        public double getPValue() {
            return this.pValue;
        }

        public boolean getHasInterfaceMutations() {
            return this.hasIntefaceMutations;
        }

        public double getFdr() {
            return this.fdr;
        }

        public void setFdr(double fdr) {
            this.fdr = fdr;
        }
    }

    private class GeneMutationProfile {
        private String geneName;
        private int interfaceLength;
        private int geneLength;
        private int interfaceMutationCount;
        private int geneMutationCount;

        GeneMutationProfile(String geneName,
                            int interfaceLength,
                            int geneLength,
                            int interfaceMutationCount,
                            int geneMutationCount) {
            this.geneName = geneName;
            this.interfaceLength = interfaceLength;
            this.geneLength = geneLength;
            this.interfaceMutationCount = interfaceMutationCount;
            this.geneMutationCount = geneMutationCount;
        }

        public double calculatePValue() {
            double ratio = (double) this.interfaceLength / (double) this.geneLength;
            if (ratio > 0.0 &&
                    this.geneMutationCount > 0 &&
                    this.interfaceMutationCount > 0) {
                try {
                    return MathUtilities.calculateBinomialPValue(
                            ratio,
                            this.geneMutationCount,
                            this.interfaceMutationCount);
                } catch (IllegalArgumentException iae) {
                    logger.error(String.format("Could not compute p-value for mutation profile %s: " +
                                    "gene length = %d," +
                                    "interface length = %d," +
                                    "gene mutation count = %d," +
                                    "interface mutation count = %d",
                            this.geneName,
                            this.geneLength,
                            this.interfaceLength,
                            this.geneMutationCount,
                            this.interfaceMutationCount,
                            iae.getMessage(),
                            Arrays.toString(iae.getStackTrace())));
                }
            }
            return 1.0d;
        }

        public String getGeneName() {
            return this.geneName;
        }

        public int getInterfaceLength() {
            return this.interfaceLength;
        }

        public int getGeneLength() {
            return this.geneLength;
        }

        public int getInterfaceMutationCount() {
            return this.interfaceMutationCount;
        }

        public int getGeneMutationCount() {
            return this.geneMutationCount;
        }

        public boolean isValid() {
            return this.geneName != null
                    && !this.geneName.equals("")
                    && this.geneLength > 0
                    && this.interfaceLength <= this.geneLength;
        }

        public String toString() {
            return String.format("%s,%d,%d,%d,%d",
                    this.geneName,
                    this.interfaceLength,
                    this.geneLength,
                    this.interfaceMutationCount,
                    this.geneMutationCount);
        }
    }

    private GeneMutationProfile getGeneMutationProfile(Map<String, Map<Integer, Set<MutationObservation>>> allSamplesGeneMap,
                                                       Map<String, Map<String, ProteinChainSummary>> interactionInterfaces,
                                                       String interactionInterface,
                                                       String geneName) {
        int interfaceLength = -1;
        int geneLength = 0;
        int interfaceMutationCount = 0;
        int geneMutationCount = 0;
        if (interactionInterfaces == null) {
            throw new IllegalStateException();
        }
        Map<String, ProteinChainSummary> mspcs = interactionInterfaces.get(interactionInterface);
        if (mspcs == null) {
            throw new IllegalStateException();
        }
        if (mspcs.size() > 0) {
            ProteinChainSummary pcs = mspcs.get(geneName);
            if (pcs == null) {
                throw new IllegalStateException();
            }
            Set<Integer> si = pcs.getInterfaceCoordinates();
            if (si == null) {
                throw new IllegalStateException();
            }
            interfaceLength = si.size();
            geneLength = pcs.getChainLength();
            if (allSamplesGeneMap.containsKey(geneName)) {
                Map<Integer, Set<MutationObservation>> mismo = allSamplesGeneMap.get(geneName);
                for (Integer aaCoord : mismo.keySet()) {
                    if (aaCoord == null) {
                        throw new IllegalStateException();
                    }
                    Set<MutationObservation> smo = mismo.get(aaCoord);
                    if (smo == null) {
                        throw new IllegalStateException();
                    }
                    geneMutationCount += smo.size();
                    Set<Integer> interfaceCoordinates = pcs.getInterfaceCoordinates();
                    if (interfaceCoordinates == null) {
                        throw new IllegalStateException();
                    }
                    if (interfaceCoordinates.contains(aaCoord)) {
                        interfaceMutationCount += smo.size();
                    }
                }
            } else {
                // not all genes that participate in interactions are found in the MAF's
                //logger.info(String.format("No mutations observed in %s",
                //        geneName));
            }
        } else {
            logger.info(String.format("Interaction '%s' gene '%s' has no interface coordinates", interactionInterface, geneName));
        }

        GeneMutationProfile geneMutationProfile = new GeneMutationProfile(geneName,
                interfaceLength,
                geneLength,
                interfaceMutationCount,
                geneMutationCount);

        if (!geneMutationProfile.isValid()) {
            logger.warn(String.format("Invalid gene mutation profile created: %s", geneMutationProfile.toString()));
        }

        return geneMutationProfile;
    }

    private Map<String, Map<String, ProteinChainSummary>> extractInterfaces(Map<String, File> interaction2PDBMap,
                                                                            Map<String, String> accessionToGene) throws IOException, StructureException {
        Interactome3dAnalyzer interactome3dAnalyzer = new Interactome3dAnalyzer();
        Map<String, Map<String, ProteinChainSummary>> interactionInterfaces = new HashMap<>();
        Map<Chain, List<Integer>> chainToContactCoordinates;
        Map<Chain, PDBUniProtMatch> chainToMatch;
        Map<String, ProteinChainSummary> interactionInterface;
        Structure structure;
        File pdbFile;

        Map<String, Integer> ProteinLengthMap = new UniProtAnalyzer().loadGeneToProteinLength();

        //sort of hate doing stuff like this...
        //TODO: make sure to clean code with refactoring recommendations later
        if (accToSeq == null) {
            ProteinSequenceHandler sequenceHandler = new ProteinSequenceHandler();
            accToSeq = sequenceHandler.loadSwissProtSequences();
        }

        String[] proteinIDs, proteinIDs2;
        for (String interaction : interaction2PDBMap.keySet()) {
            //if (interaction.equals("ZNF417\tZNF587")) {
            //    int startDebugging = 1;
            //}
            pdbFile = interaction2PDBMap.get(interaction);
            structure = StructureIO.getStructure(pdbFile.getAbsolutePath());
            chainToContactCoordinates = interactome3dAnalyzer.extractContacts(structure);
            proteinIDs = pdbFile.getName().split("-");
            proteinIDs2 = new String[]{proteinIDs[0], proteinIDs[1]};
            assert !proteinIDs.equals(proteinIDs2);
            if (!proteinIDs[2].equals("EXP")) {
                chainToMatch = interactome3dAnalyzer.mapCoordinatesToUniProtInPDB(
                        structure,
                        proteinIDs2,
                        accToSeq,
                        accessionToGene,
                        false);
            } else {
                chainToMatch = interactome3dAnalyzer.getMatchForExpStructure(
                        structure,
                        proteinIDs2,
                        accessionToGene);
            }
            //this changes chainToContactCoordinates sometimes?
            //we should return a value rather than rely on side effects
            assert chainToMatch.size() == 2;
            remapCoordinates(chainToMatch, chainToContactCoordinates);

            //TODO: kindof sloppy.. are we sure chains always map correctly to
            //genes in other methods? if so, this works
            //see Interactome3dAnalyzer.mapCoordinatesToUniProtInPDB() and
            //Interactome3dAnalyzer.getMatchForExpStructure()
            //assert chainToContactCoordinates.size() == 2;
            if(chainToContactCoordinates.size() != 2){
                logger.error(String.format(
                "Incorrect Number (%d) of Chains Detected (should be 2)",
                        chainToContactCoordinates.size()));
            }
            interactionInterface = new HashMap<>();
            try {
                if(chainToMatch == null){
                    throw new IllegalStateException();
                }
                for (Chain chainID : chainToMatch.keySet()) {
                    if(chainID == null){
                        throw new IllegalStateException();
                    }
                    String uniprotID = chainToMatch.get(chainID).getUniprot();
                    if(uniprotID == null){
                        throw new IllegalStateException();
                    }
                    String geneID = accessionToGene.get(uniprotID);
                    if(geneID == null){
                        throw new IllegalStateException();
                    }
                    Integer lengthFromPdb = chainToMatch.get(chainID).getChainResSequence().length();
                    if(lengthFromPdb == null){
                        throw new IllegalStateException();
                    }
                    Integer lengthFromMap = ProteinLengthMap.get(geneID);
                    if(lengthFromMap != null) {

                        //map geneID to set of contact coordinates for this interaction
                        if (interactionInterface == null) {
                            throw new IllegalStateException();
                        }
                        if (chainToMatch.get(chainID) == null) {
                            throw new IllegalStateException();
                        }
                        if (chainToMatch.get(chainID).getGene() == null) {
                            throw new IllegalStateException();
                        }
                        if (chainToContactCoordinates.get(chainID) == null) {
                            throw new IllegalStateException();
                        }
                        interactionInterface.put(chainToMatch.get(chainID).getGene(),
                                new ProteinChainSummary(
                                        new HashSet<>(chainToContactCoordinates.get(chainID)),
                                        lengthFromMap));
                    }else{
                        logger.warn(String.format(
                                "Detected protein without known length. " +
                                        "Uniprot ID: %s, " +
                                        "Gene ID: %s",
                                uniprotID,
                                geneID));
                    }
                }
                if(interaction == null){
                    throw new IllegalStateException();
                }
                interactionInterfaces.put(interaction, interactionInterface);
            } catch (NullPointerException npe) {
                logger.warn(String.format("Error finding interaction interface. %s: %s",
                        npe.getMessage(),
                        Arrays.toString(npe.getStackTrace())),
                        npe);
            }
        }
        return interactionInterfaces;
    }

    /**
     * This method...
     *
     * @throws Exception
     */
    private Set<String> extractInteractionsFromReactions(List<GKInstance> reactions) throws Exception {
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        Set<String> totalFIs = new HashSet<>();
        Set<String> fis;
        for (GKInstance reaction : reactions) {
            fis = reactomeAnalyzer.generateTentativePPIsForReaction(reaction,
                    false);
            if (fis.size() == 0)
                continue;
            totalFIs.addAll(fis);
        }
        return totalFIs;
    }

    private Map<String, Set<GKInstance>> generateInteractionReactionSetMap(Map<String, String> uniprotIdToGene,
                                                                           List<GKInstance> reactions) throws Exception {
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        Map<String, Set<GKInstance>> interactionReactionSetMap = new HashMap<>();
        Set<String> uniprotInteractions;
        for (GKInstance reaction : reactions) {
            uniprotInteractions = reactomeAnalyzer.generateTentativePPIsForReaction(reaction,
                    false);
            if (uniprotInteractions.size() == 0)
                continue;
            for (String uniprotInteraction : uniprotInteractions) {
                //TODO: figure out a cuter way to do this...
                String gene1 = uniprotIdToGene.get(uniprotInteraction.split("\t")[0]);
                String gene2 = uniprotIdToGene.get(uniprotInteraction.split("\t")[1]);
                if (gene1 == null || gene2 == null) {
                    //TODO: figure out which proteins don't map to genes and why
//                System.err.println("Cannot find gene for " + uniprotID);
                    continue;
                } else {
                    String geneInteraction = String.format("%s\t%s", gene1, gene2);
                    Set<GKInstance> interactionReactions;
                    if (interactionReactionSetMap.containsKey(geneInteraction)) {
                        interactionReactions = interactionReactionSetMap.get(geneInteraction);
                    } else {
                        interactionReactions = new HashSet<>();
                    }
                    interactionReactions.add(reaction);
                    interactionReactionSetMap.put(geneInteraction, interactionReactions);
                }
            }
        }
        return interactionReactionSetMap;
    }

    @Test
    public void testMafGenesToMutations() throws IOException {
        String mafDir = "datasets/firehose_all/";
        String mafPattern = "^.+\\.maf\\.txt$";
        Map<String, Map<Integer, Set<MutationObservation>>> rtn = summarizeAllMafs(mafPattern,
                mafDir);
        System.out.println(rtn.size());
    }

    /**
     * This method...
     *
     * @throws Exception
     */
    private Map<String, Map<Integer, Set<MutationObservation>>> summarizeAllMafs(String mafFileNamePattern,
                                                                                 String mafDirectoryPath) throws IOException {
        MAFFileLoader mafFileLoader = new MAFFileLoader();
        Map<String, Map<Integer, Set<MutationObservation>>> allSamplesMutatedGeneMap = new HashMap<>();
        Pattern mafFnPattern = Pattern.compile(mafFileNamePattern);
        FilenameFilter filenameFilter = (dir, name) -> {
            if (mafFnPattern.matcher(name).matches()) {
                return true;
            }
            return false;
        };

        Path path = Paths.get(mafDirectoryPath);
        List<File> mafFiles = new ArrayList<>();
        Files.walkFileTree(path, new SimpleFileVisitor<Path>() {

            @Override
            public FileVisitResult visitFile(Path file,
                                             BasicFileAttributes attrs) throws IOException {
                if (attrs.isRegularFile() && mafFnPattern.matcher(file.toFile().getName()).matches())
                    mafFiles.add(file.toFile());
                return super.visitFile(file, attrs);
            }

        });
//        System.out.println("Total maf files: " + mafFiles.size());
//        mafFiles.forEach(System.out::println);

//        File[] mafFiles = new File(mafDirectoryPath).listFiles(filenameFilter);
        Map<String, Map<Integer, Set<MutationObservation>>> sampleMutatedGeneMap;
        for (File mafFile : mafFiles) {
            sampleMutatedGeneMap = mafFileLoader.loadSampleToGenes(mafFile.getPath());
            //TODO: we can write our own data structure for hashed sets later
            Set<String> sampleGeneKeys = sampleMutatedGeneMap.keySet();
            for (String geneKey : sampleGeneKeys) {
                Map<Integer, Set<MutationObservation>> mismo = sampleMutatedGeneMap.get(geneKey);
                if (allSamplesMutatedGeneMap.containsKey(geneKey)) {
                    Map<Integer, Set<MutationObservation>> allSamplesMismo = allSamplesMutatedGeneMap.get(geneKey);
                    Set<Integer> sampleCoordKeys = mismo.keySet();
                    for (Integer coordKey : sampleCoordKeys) {
                        if (allSamplesMismo.containsKey(coordKey)) {
                            Set<MutationObservation> allSamplesSmo = allSamplesMismo.get(coordKey);
                            if (allSamplesSmo == null) {
                                throw new IllegalStateException();
                            }
                            Set<MutationObservation> smo = mismo.get(coordKey);
                            if (smo == null) {
                                throw new IllegalStateException();
                            }
                            allSamplesSmo.addAll(smo);
                            allSamplesMismo.put(coordKey, allSamplesSmo);
                        } else {
                            Set<MutationObservation> smo = mismo.get(coordKey);
                            if (smo == null) {
                                throw new IllegalStateException();
                            }
                            allSamplesMismo.put(coordKey, smo);
                        }
                    }
                    if(allSamplesMismo == null){
                        throw new IllegalStateException();
                    }
                    allSamplesMutatedGeneMap.put(geneKey, allSamplesMismo);
                } else {
                    if(mismo == null){
                        throw new IllegalStateException();
                    }
                    allSamplesMutatedGeneMap.put(geneKey, mismo);
                }
            }
        }
        for (String geneName : allSamplesMutatedGeneMap.keySet()) {
            if(geneName.equals("IRS2")){
                int x = 1;
            }
            Map<Integer, Set<MutationObservation>> mismo = allSamplesMutatedGeneMap.get(geneName);
            if (mismo == null) {
                throw new IllegalStateException();
            }
            for (Integer coord : mismo.keySet()) {
                Set<MutationObservation> smo = mismo.get(coord);
                if (smo == null) {
                    throw new IllegalStateException();
                }
            }

        }
        return allSamplesMutatedGeneMap;
    }

    /**
     * This method...
     *
     * @throws Exception
     */
    private Map<String, Map<String, Map<Integer, Set<MutationObservation>>>> summarizeIndividualMafs(String mafFileNamePattern,
                                                                                                     String mafDirectoryPath) throws IOException {
        MAFFileLoader mafFileLoader = new MAFFileLoader();
        Map<String, Map<String, Map<Integer, Set<MutationObservation>>>> allSamplesGeneMap = new HashMap<>();
        Pattern mafFnPattern = Pattern.compile(mafFileNamePattern);
        FilenameFilter filenameFilter = (dir, name) -> {
            if (mafFnPattern.matcher(name).matches()) {
                return true;
            }
            return false;
        };

        Path path = Paths.get(mafDirectoryPath);
        List<File> mafFiles = new ArrayList<>();
        Files.walkFileTree(path, new SimpleFileVisitor<Path>() {

            @Override
            public FileVisitResult visitFile(Path file,
                                             BasicFileAttributes attrs) throws IOException {
                if (attrs.isRegularFile() && mafFnPattern.matcher(file.toFile().getName()).matches())
                    mafFiles.add(file.toFile());
                return super.visitFile(file, attrs);
            }

        });
//        System.out.println("Total maf files: " + mafFiles.size());
//        mafFiles.forEach(System.out::println);

//        File[] mafFiles = new File(mafDirectoryPath).listFiles(filenameFilter);
        Map<String, Map<Integer, Set<MutationObservation>>> sampleGeneMap;
        for (File mafFile : mafFiles) {
            allSamplesGeneMap.put(mafFile.getName(), mafFileLoader.loadSampleToGenes(mafFile.getPath()));
        }
        return allSamplesGeneMap;
    }

    /**
     * This method...
     *
     * @throws Exception
     */
    private Set<String> keepGenesWithMappedProteins(Set<String> totalProteins, Map<String, String> uniprotIdToGene) {
        Set<String> totalGenes = new HashSet<>();
        for (String protein : totalProteins) {
            String gene = uniprotIdToGene.get(protein);
            if (gene == null) {
                //TODO: figure out which proteins don't map to genes and why
//                System.err.println("Cannot find gene for " + protein);
                continue;
            } else
                totalGenes.add(gene);
        }
        return totalGenes;
    }

    /**
     * This method...
     *
     * @throws Exception
     */
    private Map<String, File> getInteraction2PDBMap(String pdbDirectoryPath, Map<String, String> uniprotIdToGene) throws IOException {
        Interactome3dAnalyzer interactomeAnalyser = new Interactome3dAnalyzer();
        //comes out like "<uniprot ID>\t<uniprot ID>"
        Map<String, File> fiToPDB = interactomeAnalyser.loadPPIToPDBFile(pdbDirectoryPath,
                false);
        Map<String, File> geneSymbolToPDB = new HashMap<>();
        for (String uniprotInteraction : fiToPDB.keySet()) {
            //TODO: figure out a cuter way to do this...
            String gene1 = uniprotIdToGene.get(uniprotInteraction.split("\t")[0]);
            String gene2 = uniprotIdToGene.get(uniprotInteraction.split("\t")[1]);
            if (gene1 == null || gene2 == null) {
                //TODO: figure out which proteins don't map to genes and why
//                System.err.println("Cannot find gene for " + uniprotID);
                continue;
            } else
                geneSymbolToPDB.put(String.format("%s\t%s", gene1, gene2),
                        fiToPDB.get(uniprotInteraction));
        }
        return geneSymbolToPDB;
    }

    /**
     * Use this method to perform a one-stop analysis based on a FI file for one reaction.
     * The method will check structure PPI data downloaded from interactome3d, find its interaction
     * interfaces, load mutations downloaded from cosmic, and perform enrichment analysis.
     *
     * @throws Exception
     */
    @Test
    public void performInterfaceAnalysisForReaction() throws Exception {
        String dirName = "datasets/interactome3d/2016_06/reactions_fdr_05/";
        Map<String, File> ppiToPDB = new Interactome3dAnalyzer().loadPPIToPDBFile(dirName, true);
        System.out.println("Total ppiToPDB: " + ppiToPDB.size());        // Genes and proteins
        Set<String> genes = new HashSet<String>();
        Map<String, String> proteinToGene = new HashMap<String, String>();

        // Load FIs for the reaction
        Long dbId = 5672965L;
        dbId = 69213L;
        needDetails = true;
//        dbId = 5617896L;
        // The following files were generated by ReactomeAnalyzer.java in the FINetworkBuild project.
        // The code there should be refactored.
        //String reactionFIFile = "results/DriverGenes/Drivers_0816/FIsInReaction5617896.txt";
        String reactionFIFile = "results/DriverGenes/Drivers_0816/FIsInReaction" + dbId + ".txt";

        FileUtility fu = new FileUtility();
        fu.setInput(reactionFIFile);
        String line = fu.readLine();
        List<File> pdbFiles = new ArrayList<File>();
        System.out.println("Check for Reaction " + dbId);
        System.out.println("Found interactome3d PDB files for FIs:");
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String fi = InteractionUtilities.generateFIFromGene(tokens[0], tokens[1]);
            File file = ppiToPDB.get(fi);
            if (file != null) {
                System.out.println(line + "\t" + file.getName());
                genes.add(tokens[2]);
                genes.add(tokens[3]);
                proteinToGene.put(tokens[0], tokens[2]);
                proteinToGene.put(tokens[1], tokens[3]);
                pdbFiles.add(file);
            }
        }
        fu.close();

        Set<String> driverGenes = new CancerDriverAnalyzer().getDriverGenes(null);
        System.out.println("Total genes: " + genes.size());
        System.out.println("Total proteins: " + proteinToGene.size());
        Set<String> driverGenesinReaction = InteractionUtilities.getShared(driverGenes, genes);
        System.out.println("Driver genes in reaction: " + driverGenesinReaction.size() + ": " + driverGenesinReaction);
        CosmicAnalyzer cosmicAnalyzer = new CosmicAnalyzer();
        Map<String, List<MutationObservation>> geneToCosmicEntries = cosmicAnalyzer.loadMutations(genes);
        System.out.println("Loaded mutations for genes: " + geneToCosmicEntries.size());

        checkInterfaces(proteinToGene.keySet(),
                proteinToGene,
                geneToCosmicEntries,
                pdbFiles);
    }

    @Test
    public void generateFIsForOtherCancerReactionsAndComplexes() throws IOException {
        String dirName = "results/DriverGenes/Drivers_0816/";

        // These FIs have been analyzed via the Interactome3d Web site
        String fdr01File = dirName + "FIsInSelectedForInteractome3d_091416.txt";
        Map<String, String> fiToLineFDR01 = loadFIsWithFeaturesFile(fdr01File);
        System.out.println("Total size for fdr 0.01: " + fiToLineFDR01.size());
        String fdr0501File = dirName + "FIsInSelectedForInteractome3d_FDR_05_01_091416.txt";
        Map<String, String> fiToLineFDR0501 = loadFIsWithFeaturesFile(fdr0501File);
        System.out.println("Total size for fdr 0.05 - 0.01: " + fiToLineFDR0501.size());
        Set<String> analyzedFIs = new HashSet<String>();
        analyzedFIs.addAll(fiToLineFDR01.keySet());
        analyzedFIs.addAll(fiToLineFDR0501.keySet());
        System.out.println("Total analyzed FIs: " + analyzedFIs.size());

        // FIs extracted from cancer driver related reactions and complexes
        String fiFileName = dirName + "FIsInCancerReactionsAndComplexes_122116.txt";
        Set<String> allCancerFIs = fu.loadInteractions(fiFileName);
        System.out.println("All cancer FIs: " + allCancerFIs.size());

        // Generate a new set of FIs to send to the interactome3d group for local analysis
        Set<String> shared = InteractionUtilities.getShared(allCancerFIs, analyzedFIs);
        System.out.println("Shared: " + shared.size());
        allCancerFIs.removeAll(analyzedFIs);
        System.out.println("Not shared: " + allCancerFIs.size());
        // Output
        String outFileName = dirName + "FIsInCancerReactionsAndComplexes_Filtered_122116.txt";
        fu.saveCollection(allCancerFIs, outFileName);
    }

    @Test
    public void generateFIsForFDRBetween05And01() throws IOException {
        String dirName = "results/DriverGenes/Drivers_0816/";
        String fdr01File = dirName + "FIsInSelectedForInteractome3d_091416.txt";
        Map<String, String> fiToLineFDR01 = loadFIsWithFeaturesFile(fdr01File);
        System.out.println("Total size for fdr 0.01: " + fiToLineFDR01.size());
        String fdr0501File = dirName + "FIsInSelectedForInteractome3d_FDR_05_01_091416.txt";
        Map<String, String> fiToLineFDR0501 = loadFIsWithFeaturesFile(fdr0501File);
        System.out.println("Total size for fdr 0.05 - 0.01: " + fiToLineFDR0501.size());
        fiToLineFDR0501.keySet().removeAll(fiToLineFDR01.keySet());
        System.out.println("After removing FIs in fdr 0.01: " + fiToLineFDR0501.size());
        // Output the filtered results
        String outFileName = dirName + "FIsInSelectedForInteractome3d_FDR_05_01_filtered_091416.txt";
        fu.setInput(fdr0501File);
        String line = fu.readLine();
        fu.close();
        fu.setOutput(outFileName);
        fu.printLine(line); // Print the header
        for (String fi : fiToLineFDR0501.keySet())
            fu.printLine(fiToLineFDR0501.get(fi));
        fu.close();
    }

    private Map<String, String> loadFIsWithFeaturesFile(String fileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, String> fiToLine = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            fiToLine.put(tokens[0] + "\t" + tokens[1],
                    line);
        }
        fu.close();
        return fiToLine;
    }

    //TODO: More thoroughly document/describe this method
    private Set<String> checkInterfaces(Collection<String> accessions,
                                        Map<String, String> accessionToGene,
                                        Map<String, List<MutationObservation>> geneToCosmicEntries,
                                        Collection<File> pdbFiles) throws Exception {
        if (accToSeq == null) {
            ProteinSequenceHandler sequenceHandler = new ProteinSequenceHandler();
            accToSeq = sequenceHandler.loadSwissProtSequences();
        }

        // Handle structure
        Map<String, Set<Integer>> geneToContacts = new HashMap<String, Set<Integer>>();
        // Get the lowest p-value for single aa mutation
        Map<String, Double> geneToSinglePValue = new HashMap<String, Double>();
        Interactome3dAnalyzer interactome3dAnalyzer = new Interactome3dAnalyzer();
        interactome3dAnalyzer.setUpAtomCache();
        for (File pdbFile : pdbFiles) {
            System.out.println(pdbFile.getName());
            Structure structure = StructureIO.getStructure(pdbFile.getAbsolutePath());
            Map<Chain, List<Integer>> chainToCoordinates = interactome3dAnalyzer.extractContacts(structure);
            String[] tokens = pdbFile.getName().split("-");
            String[] acces = new String[]{tokens[0], tokens[1]};
            Map<Chain, PDBUniProtMatch> chainToMatch = null;
            if (!tokens[2].equals("EXP"))
                chainToMatch = interactome3dAnalyzer.mapCoordinatesToUniProtInPDB(structure,
                        acces,
                        accToSeq,
                        accessionToGene);
            else
                chainToMatch = interactome3dAnalyzer.getMatchForExpStructure(structure,
                        acces,
                        accessionToGene);
            remapCoordinates(chainToMatch, chainToCoordinates);
            checkInterfaces(chainToCoordinates,
                    chainToMatch,
                    geneToCosmicEntries,
                    geneToContacts,
                    geneToSinglePValue);
            System.out.println();
        }

        if (geneToLength == null)
            geneToLength = new UniProtAnalyzer().loadGeneToProteinLength();

        // Check contacts from all interactions
        System.out.println("\nGene\tLength\tContacts\tRatio\tMutations\tMutationsInContacts\tRatio\tP_Value\tp-value(single_aa_in_interface, corrected)");
        CosmicAnalyzer cosmicAnalyzer = new CosmicAnalyzer();
        Set<String> significantGenes = new HashSet<String>();
        // Choose pvalue cutoff 0.01
        double pvalueCutoff = 0.01d;
        minPValue = 1.0d; // The largest
        for (String gene : geneToContacts.keySet()) {
            Set<Integer> contacts = geneToContacts.get(gene);
            Integer length = geneToLength.get(gene);
            List<MutationObservation> mutations = geneToCosmicEntries.get(gene);
            List<MutationObservation> mutationsInContacts = cosmicAnalyzer.filterEntries(mutations, contacts);
            double ratio = contacts.size() / (double) length;
            double mutationRatio = (double) mutationsInContacts.size() / mutations.size();
            double pvalue = MathUtilities.calculateBinomialPValue(ratio, mutations.size(), mutationsInContacts.size());
            System.out.println(gene + "\t" + length + "\t" +
                    contacts.size() + "\t" + ratio + "\t" +
                    mutations.size() + "\t" + mutationsInContacts.size() + "\t" +
                    mutationRatio + "\t" + pvalue + "\t" +
                    geneToSinglePValue.get(gene));
            if (pvalue <= pvalueCutoff || geneToSinglePValue.get(gene) <= pvalueCutoff)
                significantGenes.add(gene);
            if (minPValue > pvalue)
                minPValue = pvalue;
            if (minPValue > geneToSinglePValue.get(gene))
                minPValue = geneToSinglePValue.get(gene);
        }
        return significantGenes;
    }

    private void checkInterfaces(Map<Chain, List<Integer>> chainToContactCoordiantes,
                                 Map<Chain, PDBUniProtMatch> chainToMatch,
                                 Map<String, List<MutationObservation>> geneToMutations,
                                 Map<String, Set<Integer>> geneToContacts,
                                 Map<String, Double> geneToSinglePValue) {
        System.out.println("Chain\tGene\tUniProt\tContacts\tLength\tRatio\tMutations\tMutationsInChain\tMutationsInContacts\tRatio\tp-value\tp-value(single_aa_in_interface, corrected)");
        CosmicAnalyzer cosmicHelper = new CosmicAnalyzer();
        for (Chain chain : chainToContactCoordiantes.keySet()) {
            List<Integer> contactCoords = chainToContactCoordiantes.get(chain);
            PDBUniProtMatch match = chainToMatch.get(chain);
            //TODO: Figure out why match can be null
            if (match != null) {
                List<MutationObservation> mutations = geneToMutations.get(match.getGene());
                if (mutations == null) {
                    System.out.println(match.getGene() + " doesn't have an entry in COSMIC!");
                    continue;
                }
                List<Integer> remappedSeqNumbers = remapChainSeqNumbers(chain, match);
                List<MutationObservation> mutationsInChain = cosmicHelper.filterEntries(mutations, remappedSeqNumbers);
                if (mutationsInChain.size() == 0) {
                    System.out.println(match.getGene() + " doesn't have mutations in chain!");
                    continue;
                }
                List<MutationObservation> interfaceMutations = cosmicHelper.filterEntries(mutations, contactCoords);
                if (interfaceMutations.size() == 0) {
                    System.out.println(match.getGene() + " dosn't have mutations in interface!");
                    continue;
                }
                double mutationRatio = (double) interfaceMutations.size() / mutationsInChain.size();
                double contactRatio = (double) contactCoords.size() / chain.getAtomGroups().size();
                double pvalue = 1.0;
                if (contactRatio < 1.0d) { // Otherwise, we cannot use binomial test
                    pvalue = MathUtilities.calculateBinomialPValue(contactRatio,
                            mutationsInChain.size(),
                            interfaceMutations.size());
                }
                double singleAApValue = calculatePValueForPositionMutation(interfaceMutations,
                        mutationsInChain,
                        remappedSeqNumbers);
                Double oldPValue = geneToSinglePValue.get(match.getGene());
                if (oldPValue == null || oldPValue > singleAApValue)
                    geneToSinglePValue.put(match.getGene(), singleAApValue);
                System.out.println(chain.getChainID() + "\t" +
                        match.getGene() + "\t" +
                        match.getUniprot() + "\t" +
                        contactCoords.size() + "\t" +
                        chain.getAtomGroups().size() + "\t" +
                        contactRatio + "\t" +
                        mutations.size() + "\t" +
                        mutationsInChain.size() + "\t" +
                        interfaceMutations.size() + "\t" +
                        mutationRatio + "\t" +
                        pvalue + "\t" +
                        singleAApValue);
//            System.out.println(chain.getChainID() + ": ");
                // Get mutation in the interfaces for checking
                outputMutationsInInterface(interfaceMutations, match);

                Set<Integer> geneContacts = geneToContacts.get(match.getGene());
                if (geneContacts == null) {
                    geneContacts = new HashSet<Integer>();
                    geneToContacts.put(match.getGene(), geneContacts);
                }
                geneContacts.addAll(contactCoords);
            } else {
                logger.warn(String.format("match was null for chain %s", chain));
            }
        }
    }

    private void outputMutationsInInterface(List<MutationObservation> interfaceMutations, PDBUniProtMatch match) {
        StringBuilder builder = new StringBuilder();
        builder.append("Coordinate offset: ").append(match.getOffset()).append(": ");
        Collections.sort(interfaceMutations, new Comparator<MutationObservation>() {
            public int compare(MutationObservation entry1, MutationObservation entry2) {
                return entry1.getCoordinate().compareTo(entry2.getCoordinate());
            }
        });
        for (MutationObservation entry : interfaceMutations)
            builder.append(entry.getCoordinate()).append(",").append(entry.getMutation()).append(";");
        builder.deleteCharAt(builder.length() - 1);
        System.out.println(builder.toString());
    }

    private double calculatePValueForPositionMutation(List<MutationObservation> interfaceMutations,
                                                      List<MutationObservation> chainMutations,
                                                      List<Integer> chainNumbers) {
        Map<Integer, Integer> posToMutCount = new HashMap<Integer, Integer>();
        for (MutationObservation entry : interfaceMutations) {
            Integer count = posToMutCount.get(entry.getCoordinate());
            if (count == null)
                posToMutCount.put(entry.getCoordinate(), 1);
            else
                posToMutCount.put(entry.getCoordinate(), ++count);
        }
        double pvalue = 1.0d; // Largest p-value
        double ratio = 1.0d / chainNumbers.size();
        for (Integer pos : posToMutCount.keySet()) {
            Integer count = posToMutCount.get(pos);
            double tmpPValue = MathUtilities.calculateBinomialPValue(ratio,
                    chainMutations.size(),
                    count);
            if (tmpPValue < pvalue)
                pvalue = tmpPValue;
            if (needDetails)
                System.out.println(pos + "\t" + tmpPValue);
        }
        double rtn = pvalue * posToMutCount.size(); // Perform a Bonferroni correction
        if (rtn > 1.0d)
            rtn = 1.0d;
        return rtn;
    }


    private void remapCoordinates(Map<Chain, PDBUniProtMatch> chainToMatch,
                                  Map<Chain, List<Integer>> chainToCoordinates) {
        if(chainToCoordinates == null) {
            throw new IllegalStateException();
        }
        if(chainToCoordinates.keySet().size() < 2){
            throw new IllegalStateException();
        }
        for (Chain chain : chainToCoordinates.keySet()) {
            List<Integer> coordinates = chainToCoordinates.get(chain);
//            System.out.println(chain.getChainID() + ": " + coordinates);
            PDBUniProtMatch match = chainToMatch.get(chain);
            if(match == null){
                throw new IllegalStateException();
            }
            for (int i = 0; i < coordinates.size(); i++) {
                Integer coord = coordinates.get(i);
                //TODO: figure out why NullPointerException is thrown
                if (coordinates != null && match != null) {
                    coordinates.set(i, coord + match.getOffset());
                } else {
                    logger.warn(String.format("i=%d,coordinates=%s,match=%s",
                            i, coordinates, match));
                }
            }
        }
    }

    private List<Integer> remapChainSeqNumbers(Chain chain,
                                               PDBUniProtMatch match) {
        List<Integer> list = new ArrayList<Integer>();
        List<Group> groups = chain.getAtomGroups(GroupType.AMINOACID);
        for (Group group : groups) {
            Integer old = group.getResidueNumber().getSeqNum();
            list.add(old + match.getOffset());
        }
        return list;
    }

    private Set<Integer> getSharedInterfaces(Map<String, List<Integer>> accToContacts,
                                             String[] acces) {
        if (acces.length == 1) {
            return new HashSet<Integer>(accToContacts.get(acces[0]));
        }
        Set<Integer> shared = null;
        for (int i = 0; i < acces.length - 1; i++) {
            List<Integer> contacts1 = accToContacts.get(acces[i]);
            for (int j = i + 1; j < acces.length; j++) {
                List<Integer> contacts2 = accToContacts.get(acces[j]);
                Set<Integer> shared1 = new HashSet<Integer>(contacts2);
                shared1.retainAll(contacts1);
                Set<Integer> merged = new HashSet<Integer>(contacts2);
                merged.addAll(contacts1);
                System.out.println(acces[i] + "\t" + acces[j] + "\t" +
                        contacts1.size() + "\t" + contacts2.size() + "\t" +
                        shared1.size() + "\t" + merged.size());
                if (shared == null)
                    shared = shared1;
                else
                    shared.retainAll(shared1);
            }
        }
        return shared;
    }

    private Set<Integer> getMergedInterfaces(Map<String, List<Integer>> accToContacts,
                                             String[] acces) {
        Set<Integer> merged = new HashSet<Integer>();
        for (String acc : acces)
            merged.addAll(accToContacts.get(acc));
        return merged;
    }

    /**
     * Use this method to generate interface residues.
     *
     * @throws Exception
     */
    @Test
    public void generateContacts() throws Exception {
        String pdbDir = "/Users/gwu/Documents/EclipseWorkspace/caBigR3/results/DriverGenes/Drivers_0816/interactome3D/reaction_69213/representative/interactions_1/";
        File dir = new File(pdbDir);
        String focusProtein = "P11802"; // CDK4
        // All coordinates for P11802 is based on file P11802-P55273-MDL-1blx.pdb1-A-0-B-0.pdb
        // Which has longest length for P11802, 298.
        Map<String, Integer> accToOffset = new HashMap<String, Integer>();
        accToOffset.put("P38936", 13);
        accToOffset.put("P46527", 13);
        accToOffset.put("P42771", 3);
        accToOffset.put("P42772", 3);
        String[] cdkn1 = new String[]{"P38936", "P46527"}; // CDKN1A, CDKN1B
        String[] cdkn2 = new String[]{"P42771", "P42772", "P42773", "P55273"}; // CDKN2A, B, C, D

//        focusProtein = "P38936"; // CDKN1A
//        // Use P11802-P38936-MDL-1jsu.pdb1-A-0-C-0.pdb as the reference
//        accToOffset.clear();
//        accToOffset.put("P11802", 0);
//        accToOffset.put("P24385", 159);
//        accToOffset.put("Q00534", 282);
//        // Borrowed to avoid change the code
//        cdkn1 = new String[] {"P11802", "Q00534"}; // CDK4, and CDK6
//        cdkn2 = new String[] {"P24385"}; // CCND1

        Map<String, List<Integer>> accToContacts = new HashMap<String, List<Integer>>();
        // Cache for quick performance
        AtomCache cache = new AtomCache();
        cache.setUseMmCif(true);
//        StringBuilder builder = new StringBuilder();
        for (File file : dir.listFiles()) {
            String fileName = file.getName();
            if (!fileName.endsWith(".pdb") || !fileName.contains(focusProtein))
                continue;
            if (fileName.equals("P11802-P42771-MDL-1bi7.pdb1-A-1-B-0.pdb"))
                continue; // A wrong model chosen by interactome3d as representative. This has been fixed
            // See notes on September 9, 2016.
            System.out.println(fileName);
            String[] tokens = fileName.split("-");
            Integer offset = null;
            String neededChain = null;
            if (tokens[0].equals(focusProtein)) {
                offset = accToOffset.get(tokens[1]);
                neededChain = "A";
            } else {
                offset = accToOffset.get(tokens[0]);
                neededChain = "B";
            }
            if (offset == null)
                offset = 0;
            Structure structure = StructureIO.getStructure(file.getAbsolutePath());
            Interactome3dAnalyzer helper = new Interactome3dAnalyzer();
            Map<Chain, List<Integer>> chainToContacts = helper.extractContacts(structure);
            StringBuilder builder = new StringBuilder();
            for (Chain chain : chainToContacts.keySet()) {
                if (!chain.getChainID().equals(neededChain))
                    continue;
                List<Integer> contacts = chainToContacts.get(chain);
                if (offset != 0) {
                    for (int i = 0; i < contacts.size(); i++) {
                        contacts.set(i, contacts.get(i) + offset);
                    }
                }
                String chainId = chain.getChainID();
                System.out.println(chain.getChainID());
                System.out.println(chain.getAtomGroup(0).getResidueNumber().getSeqNum() + ": " + chain.getAtomSequence().length() + ": " + chain.getAtomSequence());
                for (Integer contact : contacts) {
                    builder.append(contact).append(".").append(chainId);
                    builder.append(",");
                }
                builder.deleteCharAt(builder.length() - 1);
                System.out.println(builder.toString());
                builder.setLength(0);
                if (chain.getChainID().equals("A"))
                    accToContacts.put(tokens[1], contacts);
                else
                    accToContacts.put(tokens[0], contacts);
            }
            System.out.println();
        }
//        if (true)
//            return;
        // Compare between same proteins and different proteins
        Set<Integer> cdkn1All = new HashSet<Integer>();
        Set<Integer> cdkn1Shared = getSharedInterfaces(accToContacts, cdkn1);
        System.out.println("CDKN1 Shared: " + cdkn1Shared.size());
        Set<Integer> cdkn1Merged = getMergedInterfaces(accToContacts, cdkn1);
        System.out.println("CDKN1 merged: " + cdkn1Merged.size());
        outputInterfaces(cdkn1Merged, "A");
        Set<Integer> cdkn2Shared = getSharedInterfaces(accToContacts, cdkn2);
        System.out.println("CDKN2 Shared: " + cdkn2Shared.size());
        Set<Integer> cdkn2Merged = getMergedInterfaces(accToContacts, cdkn2);
        System.out.println("CDKN2 merged: " + cdkn2Merged.size());
        outputInterfaces(cdkn2Merged, "A");
        Set<Integer> allShared = new HashSet<Integer>(cdkn2Merged);
        allShared.retainAll(cdkn1Merged);
        System.out.println("All shared: " + allShared.size());
        outputInterfaces(allShared, "A");
    }

    private void outputInterfaces(Collection<Integer> contacts, String chain) {
        StringBuilder builder = new StringBuilder();
        for (Integer c : contacts)
            builder.append(c).append(".").append(chain).append(",");
        builder.deleteCharAt(builder.length() - 1);
        System.out.println(builder.toString());
    }

}
