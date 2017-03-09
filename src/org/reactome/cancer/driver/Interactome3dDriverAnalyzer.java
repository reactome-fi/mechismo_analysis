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
import org.reactome.r3.CosmicAnalyzer.CosmicEntry;
import org.reactome.r3.Interactome3dAnalyzer.PDBUniProtMatch;
import org.reactome.r3.ProteinSequenceHandler.Sequence;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.MathUtilities;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.*;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

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
        Map<String, List<CosmicEntry>> geneToCosmicEntries = cosmicAnalyzer.loadMutations(matchedGenes);

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
     * This method takes a reaction directly from the Reactome database, expand is into a set
     * of interactions in genes, load structures using interactome3d data, and then perform mutation
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
            Set<String> fis = reactomeDataAnalyzer.generateAttentativePPIsForReaction(reaction,
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
        // "<gene symbol>" -> [<CosmicEntry>]
        CosmicAnalyzer cosmicAnalyzer = new CosmicAnalyzer();
        Map<String, List<CosmicEntry>> geneToCosmicEntries = cosmicAnalyzer.loadMutations(totalGenes);
        System.out.println("Total genes in cosmic: " + geneToCosmicEntries.size());

        // Check interactions one by one
        Map<GKInstance, Double> reactionToMinPValue = new HashMap<>();
        Map<GKInstance, Integer> reactionToFINumber = new HashMap<>();
        Map<GKInstance, Set<String>> reactionToSigGenes = new HashMap<>();
        for (GKInstance reaction : reactions) {
            // Store interactions in a set
            // "<uniprot ID>\t<uniprot ID>"
            Set<String> fis = reactomeDataAnalyzer.generateAttentativePPIsForReaction(reaction,
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
     * This method overloads {@link #findInteractionsWithMutatedInterfaces(CancerDriverReactomeAnalyzer,String,String,String,String)
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

    /**
     * This method...
     *
     * @throws Exception
     */
    public void findInteractionsWithMutatedInterfaces(CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer,
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
        Set<String> totalFIs = interactionSet(reactions);
        System.out.println(String.format("Total FIs: %d", totalFIs.size()));

        // Map uniprot ID's to gene symbols (COSMIC uses gene symbols)
        // "<uniprot ID>" -> "<gene symbol>"
        UniProtAnalyzer uniprotAnalyzer = new UniProtAnalyzer();
        Map<String, String> uniprotIdToGene = uniprotAnalyzer.getUniProtAccessionToGeneName();

        // Extract and store unique proteins from interactions in a set
        // "<uniprot ID>"
        Set<String> totalProteins = InteractionUtilities.grepIDsFromInteractions(totalFIs);
        System.out.println(String.format("Total proteins: %d", totalProteins.size()));

        //Keep only Reactome FI's containing >= 1 gene mutated in MAF
        // Lookup and store gene symbols for proteins in a set
        // "<gene symbol>"
        Set<String> totalGenes = mutatedFIs(totalProteins,uniprotIdToGene);
        System.out.println(String.format("Total genes: %d", totalGenes.size()));

        //MAF Mutations
        Map<String,Map<String,Integer>> allSamplesGeneMap = mafGeneToMutation(mafFileNamePattern,mafDirectoryPath);

        //TODO: ensure this works... caps? some genes have multiple names (Akt/PKB) etc.
        allSamplesGeneMap.keySet().retainAll(totalGenes);

        // Map all pdb files to gene names;
        // "<gene symbol>\t<uniprot ID>" -> "path/to/complex.pdb"
        Map<String,File>geneSymbolToPDB = geneSymbolPdbMap(pdbDirectoryPath,uniprotIdToGene);
        System.out.println(String.format("Total PDB's: %d", geneSymbolToPDB.size()));

        //Reactome FI's Interface Mutation Ratio
        for(String mutatedGene : allSamplesGeneMap.keySet()){
            //do something like checkInterfaces with allSamplesGeneMap.get(mutatedGene) and geneSymbolToPDB.get(mutatedGene)
        }
    }

    /**
     * This method...
     *
     * @throws Exception
     */
    private Set<String> interactionSet(List<GKInstance> reactions) throws Exception {
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        Set<String> totalFIs = new HashSet<>();
        Set<String> fis;
        for (GKInstance reaction : reactions) {
            fis = reactomeAnalyzer.generateAttentativePPIsForReaction(reaction,
                    false);
            if (fis.size() == 0)
                continue;
            totalFIs.addAll(fis);
        }
        return totalFIs;
    }

    /**
     * This method...
     *
     * @throws Exception
     */
    private Map<String,Map<String,Integer>> mafGeneToMutation(String mafFileNamePattern,String mafDirectoryPath) throws IOException {
        MAFFileLoader mafFileLoader = new MAFFileLoader();
        Map<String,Map<String,Integer>> allSamplesGeneMap = new HashMap<>();
        Pattern mafFnPattern = Pattern.compile(mafFileNamePattern);
        FilenameFilter filenameFilter = (dir, name) -> {
            if(mafFnPattern.matcher(name).matches()){
                return true;
            }
            return false;
        };
        File[] mafFiles = new File(mafDirectoryPath).listFiles(filenameFilter);
        Map<String, Map<String,Integer>> sampleGeneMap;
        for (File mafFile : mafFiles) {
            sampleGeneMap = mafFileLoader.loadSampleToGenes(mafFile.getPath());
            //TODO: we can write our own data structure for hashed sets later
            for (String geneKey : sampleGeneMap.keySet()) {
                if (allSamplesGeneMap.containsKey(geneKey)) {
                    for (String mutationKey : sampleGeneMap.get(geneKey).keySet()) {
                        if (allSamplesGeneMap.get(geneKey).containsKey(mutationKey)) {
                            Map<String, Integer> map = allSamplesGeneMap.get(geneKey);
                            map.put(mutationKey, map.get(mutationKey) + 1);
                            sampleGeneMap.put(geneKey, map);
                        } else {
                            Map<String, Integer> map = allSamplesGeneMap.get(geneKey);
                            map.put(mutationKey, 1);
                            sampleGeneMap.put(geneKey, map);
                        }
                    }
                } else {
                    allSamplesGeneMap.put(geneKey, sampleGeneMap.get(geneKey));
                }
            }
        }
        return allSamplesGeneMap;
    }

    /**
     * This method...
     *
     * @throws Exception
     */
    private Set<String> mutatedFIs(Set<String> totalProteins,Map<String, String> uniprotIdToGene){
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
    private Map<String,File> geneSymbolPdbMap(String pdbDirectoryPath,Map<String, String> uniprotIdToGene) throws IOException {
        Interactome3dAnalyzer interactomeAnalyser = new Interactome3dAnalyzer();
        Map<String, File> fiToPDB = interactomeAnalyser.loadPPIToPDBFile(pdbDirectoryPath,
                false);
        Map<String, File> geneSymbolToPDB = new HashMap<>();
        for (String uniprotID : fiToPDB.keySet()) {
            String gene = uniprotIdToGene.get(uniprotID);
            if (gene == null) {
                //TODO: figure out which proteins don't map to genes and why
//                System.err.println("Cannot find gene for " + uniprotID);
                continue;
            } else
                geneSymbolToPDB.put(gene,fiToPDB.get(uniprotID));
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
        Map<String, List<CosmicEntry>> geneToCosmicEntries = cosmicAnalyzer.loadMutations(genes);
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
                                        Map<String, List<CosmicEntry>> geneToCosmicEntries,
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
            List<CosmicEntry> mutations = geneToCosmicEntries.get(gene);
            List<CosmicEntry> mutationsInContacts = cosmicAnalyzer.filterEntries(mutations, contacts);
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
                                 Map<String, List<CosmicEntry>> geneToMutations,
                                 Map<String, Set<Integer>> geneToContacts,
                                 Map<String, Double> geneToSinglePValue) {
        System.out.println("Chain\tGene\tUniProt\tContacts\tLength\tRatio\tMutations\tMutationsInChain\tMutationsInContacts\tRatio\tp-value\tp-value(single_aa_in_interface, corrected)");
        CosmicAnalyzer cosmicHelper = new CosmicAnalyzer();
        for (Chain chain : chainToContactCoordiantes.keySet()) {
            List<Integer> contactCoords = chainToContactCoordiantes.get(chain);
            PDBUniProtMatch match = chainToMatch.get(chain);
            //TODO: Figure out why match can be null
            if (match != null) {
                List<CosmicEntry> mutations = geneToMutations.get(match.getGene());
                if (mutations == null) {
                    System.out.println(match.getGene() + " doesn't have an entry in COSMIC!");
                    continue;
                }
                List<Integer> remappedSeqNumbers = remapChainSeqNumbers(chain, match);
                List<CosmicEntry> mutationsInChain = cosmicHelper.filterEntries(mutations, remappedSeqNumbers);
                if (mutationsInChain.size() == 0) {
                    System.out.println(match.getGene() + " doesn't have mutations in chain!");
                    continue;
                }
                List<CosmicEntry> interfaceMutations = cosmicHelper.filterEntries(mutations, contactCoords);
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

    private void outputMutationsInInterface(List<CosmicEntry> interfaceMutations, PDBUniProtMatch match) {
        StringBuilder builder = new StringBuilder();
        builder.append("Coordinate offset: ").append(match.getOffset()).append(": ");
        Collections.sort(interfaceMutations, new Comparator<CosmicEntry>() {
            public int compare(CosmicEntry entry1, CosmicEntry entry2) {
                return entry1.getCoordinate().compareTo(entry2.getCoordinate());
            }
        });
        for (CosmicEntry entry : interfaceMutations)
            builder.append(entry.getCoordinate()).append(",").append(entry.getMutation()).append(";");
        builder.deleteCharAt(builder.length() - 1);
        System.out.println(builder.toString());
    }

    private double calculatePValueForPositionMutation(List<CosmicEntry> interfaceMutations,
                                                      List<CosmicEntry> chainMutations,
                                                      List<Integer> chainNumbers) {
        Map<Integer, Integer> posToMutCount = new HashMap<Integer, Integer>();
        for (CosmicEntry entry : interfaceMutations) {
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
        for (Chain chain : chainToCoordinates.keySet()) {
            List<Integer> coordinates = chainToCoordinates.get(chain);
//            System.out.println(chain.getChainID() + ": " + coordinates);
            PDBUniProtMatch match = chainToMatch.get(chain);
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
