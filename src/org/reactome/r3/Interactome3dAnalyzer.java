/*
 * Created on Sep 22, 2016
 *
 */
package org.reactome.r3;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.contact.AtomContactSet;
import org.biojava.nbio.structure.contact.GroupContact;
import org.biojava.nbio.structure.contact.GroupContactSet;
import org.biojava.nbio.structure.contact.Pair;
import org.gk.util.GKApplicationUtilities;
import org.junit.Test;
import org.reactome.r3.ProteinSequenceHandler.Sequence;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.structure.model.PDBUniProtMatch;

/**
 * This class is used to handle interaction structure data downloaded from interactome3d.
 * @author gwu
 *
 */
public class Interactome3dAnalyzer {
    
    private final double INTERFACE_THRESHOLD = 5.0d; // The unit is Angstrom (A)
    private FileUtility fu = new FileUtility();
    private final static Logger logger = Logger.getLogger(Interactome3dAnalyzer.class);
    
    /**
     * Default constructor.
     */
    public Interactome3dAnalyzer() {
    }
    
    /**
     * If StructureIO.getSturcture() cannot work, check temp dir used by this method
     *  private String initPdbFilePath() in biojava class, UserConfiguration, and then
     *  deleted these two folders there: data and cheminfo.
     * @throws Exception
     */
    @Test
    public void testBioJava() throws Exception {
        String fileName = "datasets/interactome3d/2017_01/representative/interactions_06/P00533-P01133-EXP-1nql.pdb1-A-0-B-0.pdb";
//        setUpAtomCache();
        Structure structure = StructureIO.getStructure(fileName);
//        Structure structure = StructureIO.getStructure("4HHB");
        System.out.println("Structure: " + structure.getIdentifier());
        List<Chain> chains = structure.getChains();
        chains.forEach(chain -> System.out.println(chain.getChainID() + ": " + chain.getAtomSequence()));
    }
    
    /**
     * Compare if there are overlapping structures for FIs in two interaction3d queries.
     * @throws IOException
     */
    @Test
    public void compareTwoSetsOfInteractions() throws IOException {
        String dirName = "datasets/interactome3d/2016_06/reactions_fdr_01/";
        List<File> files01 = getInteractionStructureFiles(dirName + "complete/interactions",
                                                                true);
        System.out.println(dirName);
        System.out.println("Total files: " + files01.size());
        Set<String> fis01 = new HashSet<String>();
        for (File file : files01)
            fis01.add(extractFIFromFileName(file));
        System.out.println("Total FIs: " + fis01.size());
        dirName = "datasets/interactome3d/2016_06/reactions_fdr_05_01/";
        List<File> files05_01 = getInteractionStructureFiles(dirName + "complete/interactions",
                                                             true);
        System.out.println(dirName);
        System.out.println("Total fils: " + files05_01.size());
        Set<String> fis05_01 = new HashSet<String>();
        for (File file : files05_01)
            fis05_01.add(extractFIFromFileName(file));
        System.out.println("Total FIs: " + fis05_01.size());
        Set<String> shared = InteractionUtilities.getShared(fis01, fis05_01);
        System.out.println("Shared interactions: " + shared.size());
    }
    
    private Map<String, Integer> loadRank(String dir) throws IOException {
        String fileName = dir + File.separator + "interactions.dat";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, Integer> fileToRank = new HashMap<String, Integer>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            fileToRank.put(tokens[tokens.length - 1],
                           new Integer(tokens[2]));
        }
        fu.close();
        return fileToRank;
    }
    
    @Test
    public void pickUpBestStructuresBasedOnContacts() throws Exception {
        // To hold all picked PDBs
        String targetDir = "datasets/interactome3d/2016_06/reactions_fdr_05/";
        // Source dir
        String dirName = "datasets/interactome3d/2016_06/reactions_fdr_01/";
        dirName = "datasets/interactome3d/2016_06/reactions_fdr_05_01/";
//        String representDir = dirName + "representative";
//        List<File> representFiles = getInteractionStructureFiles(representDir, 
//                                                                 true);
//        System.out.println("Total interactions in " + dirName + ": ");
//        System.out.println("\trepresentative: " + representFiles.size());
        String completeDir = dirName + "complete";
        Map<String, Integer> fileToRank = loadRank(completeDir);
        List<File> completeFiles = getInteractionStructureFiles(completeDir + File.separator + "interactions", 
                                                                true);
        System.out.println("\tcomplete: " + completeFiles.size());
        
        List<File> allFiles = new ArrayList<File>();
//        allFiles.addAll(representFiles);
        allFiles.addAll(completeFiles);
        Map<String, List<InteractionPDBScore>> fiToScores = new HashMap<String, List<InteractionPDBScore>>();
        for (File file : allFiles) {
            InteractionPDBScore score = calculateInteractionScore(file,
                                                                  fileToRank);
            String fi = extractFIFromFileName(file);
            List<InteractionPDBScore> scores = fiToScores.get(fi);
            if (scores == null) {
                scores = new ArrayList<InteractionPDBScore>();
                fiToScores.put(fi, scores);
            }
            scores.add(score);
        }
        System.out.println("Total FIs: " + fiToScores.size());
        StringBuilder builder = new StringBuilder();
        for (String fi : fiToScores.keySet()) {
            List<InteractionPDBScore> scores = fiToScores.get(fi);
            Collections.sort(scores);
            builder.append(fi).append("\t");
            for (InteractionPDBScore score : scores) {
                builder.append(score.getShortFileName()).append("\t");
//                builder.append(score.isInRepresent).append("\t");
                builder.append(score.interactResidues).append("\t");
            }
            builder.deleteCharAt(builder.length() - 1);
            System.out.println(builder.toString());
            builder.setLength(0);
            // Copy the top one
            InteractionPDBScore topScore = scores.get(0);
            File targetFile = new File(targetDir, topScore.getShortFileName());
            GKApplicationUtilities.copy(new File(topScore.fileName),
                                        targetFile);
        }
    }
    
    public void setUpAtomCache() {
        AtomCache cache = new AtomCache();
        cache.setUseMmCif(true);
        StructureIO.setAtomCache(cache);
    }
    
    private InteractionPDBScore calculateInteractionScore(File file,
                                                          Map<String, Integer> fileToRank) throws Exception {
        Structure structure = StructureIO.getStructure(file.getAbsolutePath());
        Map<Chain, List<Integer>> chainToResidues = extractContacts(structure);
        int total = 0;
        for (Chain chain : chainToResidues.keySet()) {
            List<Integer> list = chainToResidues.get(chain);
            total += list.size();
        }
        // Get the total residues in the interface
        String fi = extractFIFromFileName(file);
        InteractionPDBScore score = new InteractionPDBScore();
        score.fileName = file.getAbsolutePath();
        score.interactResidues = total;
        score.isInRepresent = file.getAbsolutePath().contains("representative");
        score.rank = fileToRank.get(file.getName());
        return score;
    }
    
    @Test
    public void testLoadPPIToPDBFile() throws IOException {
        String dirName = "datasets/interactome3d/2016_06/prebuilt/representative/";
        Map<String, File> ppiToPDB = loadPPIToPDBFile(dirName, false);
        System.out.println("Total PPIToPDB: " + ppiToPDB.size());
        int count = 0;
        for (String ppi : ppiToPDB.keySet()) {
            File pdb = ppiToPDB.get(ppi);
            System.out.println(ppi + ": " + pdb.getName());
            count ++;
            if (count == 10)
                break;
        }
    }
    
    private List<File> getInteractionStructureFiles(String dirName,
                                                    boolean recursive) throws IOException {
        return getPDBFiles(dirName, recursive);
    }
    
    private List<File> getPDBFiles(String dirName,
                                   boolean recursive) throws IOException {
        List<File> interactionPDBFiles = new ArrayList<File>();
        List<File> dirs = new ArrayList<File>();
        dirs.add(new File(dirName));
        List<File> next = new ArrayList<File>();
        // At least once
        while (dirs.size() > 0) {
            for (File dir : dirs) {
                File[] files = dir.listFiles();
                if (files == null || files.length == 0)
                    continue;
                for (File file : files) {
                    if (file.isDirectory()) {
                        next.add(file);
                        continue;
                    }
                    if (!file.getName().endsWith(".pdb"))
                        continue;
                    interactionPDBFiles.add(file);
                }
            }
            if (!recursive)
                break; // Perform once only
            dirs.clear();
            dirs.addAll(next);
            next.clear();
        }
        return interactionPDBFiles;
    }
    
    public Map<String, File> loadPPIToPDBFile(String dirName,
                                              boolean recursive) throws IOException {
        Map<String, File> ppiToPDB = new HashMap<String, File>();
        List<File> interactionPDBFiles = getInteractionStructureFiles(dirName,
                                                                      recursive);
        for (File file : interactionPDBFiles) {
            // Do a simple parsing
            String fi = extractFIFromFileName(file);
            String[] proteins = fi.split("\t");
            assert proteins.length == 2;
            if(!proteins[0].equals("EXP") && !proteins[0].equals("MDL") &&
                    !proteins[1].equals("EXP") && !proteins[1].equals("MDL")) {
                if (ppiToPDB.containsKey(fi))
                    // TODO: simple parsing isn't good enough
                    // we end up with a bunch of "MDL" and "EXP" entries that
                    // shouldn't exist... should only contain uniprot ID's
                    // we'll fix this later with a regex like FIOncoNet
                    //throw new IllegalStateException("Duplicated PDB for " + fi);
                    logger.warn(String.format("ppiToPDB already contains %s", fi));
                ppiToPDB.put(fi, file);
            }else{
                //skip single-protein EXP and MDL structures
            }
        }
        return ppiToPDB;
    }
    
    /**
     * In the complete data set, one PPI may have multiple PDB files. It has been found
     * that the chosen PDF files in the representative may not be good. Therefore, it will
     * be better to perform analysis by loading all PDB files for a PPI.
     * @param dirName
     * @param recursive
     * @return
     * @throws IOException
     */
    public Map<String, List<File>> loadPPIToCompletePDBFiles(String dirName,
                                                             boolean recursive) throws IOException {
        Map<String, List<File>> ppiToPDBs = new HashMap<>();
        List<File> interactionPDBFiles = getInteractionStructureFiles(dirName,
                                                                      recursive);
        interactionPDBFiles.forEach(file -> {
            String fi = extractFIFromFileName(file);
            ppiToPDBs.compute(fi, (key, list) -> {
                if (list == null)
                    list = new ArrayList<>();
                list.add(file);
                return list;
            });
        });
        return ppiToPDBs;
    }
    
    @Test
    public void testLoadPPIToCompletePDBFiles() throws IOException {
        String dirName = "datasets/interactome3d/2017_01/complete";
        Map<String, List<File>> ppiToPDBs = loadPPIToCompletePDBFiles(dirName, true);
        System.out.println("Total PPIs: " + ppiToPDBs.size());
//        ppiToPDBs.forEach((ppi, list) -> System.out.println(ppi + ": " + list.size()));
        
        // Compare with representative dataset
        dirName = "datasets/interactome3d/2017_01/representative";
        Map<String, List<File>> repPPIToPDBs = loadPPIToCompletePDBFiles(dirName, true);
        System.out.println("Total PPIs in representative: " + repPPIToPDBs.size());
        
       ppiToPDBs.forEach((ppi, list) -> {
           List<File> repPDBs = repPPIToPDBs.get(ppi);
           List<String> fileNames = list.stream().map(file -> file.getName()).collect(Collectors.toList());
           List<String> repFileNames = repPDBs.stream().map(file -> file.getName()).collect(Collectors.toList());
           // There is only one case: P42345    Q8N122, caused different letter case usage.
           if (!fileNames.containsAll(repFileNames)) {
               System.out.println(ppi + " is not contained in complete!");
               System.out.println("Complete: " + fileNames);
               System.out.println("Representative: " + repFileNames);
           }
       });
    }
    
    /**
     * Load the map from proteins to their structures. One protein may have multiple structures in
     * either representative or complete folder.
     * @param dirName
     * @return
     * @throws IOException
     */
    public Map<String, Set<String>> loadProteinToPDBFiles(String dirName) throws IOException {
        List<File> proteinPDBFiles = getPDBFiles(dirName + File.separator + "proteins", 
                                                 true);
        Map<String, Set<String>> proteinToPDBs = new HashMap<String, Set<String>>();
        for (File pdbFile : proteinPDBFiles) {
            String fileName = pdbFile.getName();
            int index = fileName.indexOf("-");
            String protein = fileName.substring(0, index);
            InteractionUtilities.addElementToSet(proteinToPDBs, protein, pdbFile.getAbsolutePath());
        }
        return proteinToPDBs;
    }

    private String extractFIFromFileName(File file) {
        String[] tokens = file.getName().split("-");
        String fi = InteractionUtilities.generateFIFromGene(tokens[0], tokens[1]);
        return fi;
    }
    
    @Test
    public void downloadData() throws Exception {
        String interactome3dDir = "datasets/interactome3d/2016_06/";
        interactome3dDir = "datasets/interactome3d/2017_01/";
        
        String dirName = interactome3dDir + "representative/";
        String hostUrl = "http://interactome3d.irbbarcelona.org/user_data/human/download/representative/";
        int largestInt = 23;
//        
        dirName = interactome3dDir + "complete/";
        hostUrl = "http://interactome3d.irbbarcelona.org/user_data/human/download/complete/";
        largestInt = 149;
        
        // Results for reactions having enrichment fdrs <= 0.01
        // Representative
//        String dirName = "datasets/interactome3d/2016_06/reactions_fdr_01/representative/";
//        String hostUrl = "http://interactome3d.irbbarcelona.org/user_data/EXggwX2u6XZQVKWMpuBQ/download/representative/";
//        int largestInt = 1;
//        int largestProtein = 2;
//        // Complete
//        dirName = "datasets/interactome3d/2016_06/reactions_fdr_01/complete/";
//        hostUrl = "http://interactome3d.irbbarcelona.org/user_data/EXggwX2u6XZQVKWMpuBQ/download/complete/";
//        largestInt = 2;
//        largestProtein = 26;
//        // Results for reactions having encirhment fdr [0.05, 0.01)
//        // representative
//        dirName = "datasets/interactome3d/2016_06/reactions_fdr_05_01/representative/";
//        hostUrl = "http://interactome3d.irbbarcelona.org/user_data/CEtGowWWxcDu6bxu5rm6/download/representative/";
//        largestInt = 2;
//        largestProtein = 2;
//        // Complete
//        dirName = "datasets/interactome3d/2016_06/reactions_fdr_05_01/complete/";
//        hostUrl = "http://interactome3d.irbbarcelona.org/user_data/CEtGowWWxcDu6bxu5rm6/download/complete/";
//        largestInt = 7;
//        largestProtein = 28;
        
        String fileName = "interactions.dat";
        download(dirName, hostUrl, fileName);
//        fileName = "proteins.dat";
//        download(dirName, hostUrl, fileName);
        downloadDataFiles(dirName, hostUrl, "interactions", largestInt);
//        downloadDataFiles(dirName, hostUrl, "proteins", largestProtein);
        // After downloading, use a simple script as following:
        // $ for f in *.tar; do tar xf $f; done (http://stackoverflow.com/questions/583889/how-can-you-untar-more-than-one-file-at-a-time)
        // Or use the following script (preferred) to keep the original folder
//        for f in *.tgz; do d=`basename $f .tgz`; mkdir $d; cd $d; tar xzvf ../$f; cd ..; done
    }

    private void downloadDataFiles(String dirName, 
                                   String hostUrl,
                                   String dataType,
                                   int largestInt) throws Exception {
        for (int i = 1; i <= largestInt; i++) {
            String inset = "";
            if (largestInt >= 100) {
                if (i < 10)
                    inset = "00";
                else if (i < 100)
                    inset = "0";
            }
            else if (largestInt >= 10) {
                if (i < 10)
                    inset = "0";
            }
            String fileName = dataType + "_" + inset + i + ".tgz";
            download(dirName, hostUrl, fileName);
        }
    }

    private void download(String dirName, String hostUrl, String fileName) throws Exception {
        System.out.println("Download: " + fileName);
        String dataUrl = hostUrl + fileName;
        URL url = new URL(dataUrl);
        InputStream is = url.openStream();
        FileOutputStream fos = new FileOutputStream(dirName + fileName);
        byte[] buffer = new byte[1024 * 100];
        int read = is.read(buffer, 0, buffer.length);
        while (read > 0) {
            fos.write(buffer, 0, read);
            read = is.read(buffer, 0, buffer.length);
        }
        is.close();
        fos.close();
    }
    
    @Test
    public void checkInterfaces() throws Exception {
//        String dir = "datasets/interactome3d/2016_06/reactions_fdr_01/complete/interactions_2/";
//        String fileName = dir + "P20393-Q07869-MDL-1a6y.pdb1-B-0-A-0.pdb";
        String fileName = "datasets/interactome3d/2017_01/representative/interactions_06/P00533-P01133-EXP-1nql.pdb1-A-0-B-0.pdb";
        setUpAtomCache();
//        fileName = "results/DriverGenes/Drivers_0816/interactome3D/reaction_69213/representative/interactions_1/P11802-P42771-MDL-1bi7.pdb1-A-0-B-0.pdb";
        Structure structure = StructureIO.getStructure(fileName);
        Map<Chain, List<Integer>> chainToResidues = extractContacts(structure);
        for (Chain chain : chainToResidues.keySet()) {
            System.out.println(chain.getChainID());
            Group firstGroup = chain.getAtomGroups().get(0);
            int startCoord = firstGroup.getResidueNumber().getSeqNum();
            System.out.println("First coordinate: " + startCoord);
            List<Integer> residues = chainToResidues.get(chain);
            System.out.println(residues);
        }
    }
    
    public Map<Chain, List<Integer>> extractContacts(Chain chainA,
                                                     Chain chainB) throws StructureException {
        AtomContactSet contacts = StructureTools.getAtomsInContact(chainA,
                                                                   chainB,
                                                                   INTERFACE_THRESHOLD,
                                                                   false);
        GroupContactSet groupSet = new GroupContactSet(contacts);
        Map<Chain, Set<Integer>> chainToCoordinates = new HashMap<>();
        Map<Chain, List<Integer>> rtn = new HashMap<>();
        if(groupSet.size() > 0) {
            for (GroupContact contact : groupSet) {
                Pair<Group> pair = contact.getPair();
                Set<Integer> coordinates = chainToCoordinates.get(pair.getFirst().getChain());
                if (coordinates == null) {
                    coordinates = new HashSet<>();
                    chainToCoordinates.put(pair.getFirst().getChain(),
                            coordinates);
                }
                coordinates.add(pair.getFirst().getResidueNumber().getSeqNum());
                coordinates = chainToCoordinates.get(pair.getSecond().getChain());
                if (coordinates == null) {
                    coordinates = new HashSet<>();
                    chainToCoordinates.put(pair.getSecond().getChain(),
                            coordinates);
                }
                coordinates.add(pair.getSecond().getResidueNumber().getSeqNum());
            }
            for (Chain chain : chainToCoordinates.keySet()) {
                Set<Integer> set = chainToCoordinates.get(chain);
                List<Integer> list = new ArrayList<>(set);
                Collections.sort(list);
                rtn.put(chain, list);
            }
        }else{
            logger.info(String.format("Detected structure without %.2fA interface.",
                    INTERFACE_THRESHOLD));
            rtn.put(chainA,new ArrayList<>());
            rtn.put(chainB,new ArrayList<>());
        }
        return rtn;
    }
    
    public Map<Chain, List<Integer>> extractContacts(Structure structure) throws StructureException {
        // Assume there are two chains, named as A and B. This seems always true for
        // PDBs downloaded from Interactome3D.
        Chain chainA = structure.getChainByPDB("A");
        Chain chainB = structure.getChainByPDB("B");
        if(chainA == null || chainB == null){
            throw new NullPointerException();
        }
        return extractContacts(chainA, chainB);
    }
    
    public Map<Chain, PDBUniProtMatch> getMatchForExpStructure(Structure structure,
                                                                String[] accessions,
                                                                Map<String, String> accToGene) throws StructureException {
        if (accessions.length < 2)
            throw new IllegalArgumentException("The passed accessions array should have at least 2 elements!");
        Map<Chain, PDBUniProtMatch> chainToMatch = new HashMap<Chain, PDBUniProtMatch>();
        for (int i = 0; i < 2; i++) {
            String chainId = (i == 0 ? "A" : "B");
            Chain chain = structure.getChainByPDB(chainId);
            PDBUniProtMatch match = new PDBUniProtMatch();
            match.setChainID(chainId);
            match.setChainSequence(chain.getAtomSequence());
            match.setUniprot(accessions[i]);
            match.setGene(accToGene.get(match.getUniprot()));
            match.setOffset(0);
            chainToMatch.put(chain, match);
        }
        return chainToMatch;
    }
    
    @Test
    public void testMapCoordinatesToUniProtInPDB() throws Exception {
        ProteinSequenceHandler sequenceHandler = new ProteinSequenceHandler();
        Map<String, Sequence> accToSeq = sequenceHandler.loadSwissProtSequences();
        Map<String, String> accToGene = new UniProtAnalyzer().getUniProtAccessionToGeneName();
        String interactomeDirName = "datasets/interactome3d/2016_06/prebuilt/representative/";
        File pdbFile = new File(interactomeDirName + "P01709-P0CG04-MDD-V-set-C1-set-5dur-D-3-111-L-125-212.pdb");
        Structure structure = StructureIO.getStructure(pdbFile.getAbsolutePath());
        String[] tokens = pdbFile.getName().split("-");
        String[] acces = new String[]{tokens[0], tokens[1]};
        Map<Chain, PDBUniProtMatch> chainToMatch = mapCoordinatesToUniProtInPDB(structure,
                                                                                acces, 
                                                                                accToSeq,
                                                                                accToGene);
        for (Chain chain : chainToMatch.keySet()) {
            PDBUniProtMatch match = chainToMatch.get(chain);
            System.out.println(chain.getChainID() + ": " + match.getOffset());
        }
    }

    public void remapCoordinates(Map<Chain, PDBUniProtMatch> chainToMatch,
                                 Map<Chain, List<Integer>> chainToCoordinates) {
        for (Chain chain : chainToCoordinates.keySet()) {
            List<Integer> coordinates = chainToCoordinates.get(chain);
            //System.out.println(chain.getChainID() + ": " + coordinates);
            PDBUniProtMatch match = chainToMatch.get(chain);
            if (match == null)
                throw new IllegalStateException("Cannot find match for chain: " + chain);
            for (int i = 0; i < coordinates.size(); i++) {
                Integer coord = coordinates.get(i);
                coordinates.set(i, coord + match.getOffset());
            }
        }
    }

    /**
     * Note: removePrintLines in the method. Don't add debug related control as the parameter. Think this is used as
     * a public API.
     * @param structure
     * @param uniprotIDs
     * @param uniprotToSeq
     * @param uniprotToGene
     * @param printLines
     * @return
     * @throws IOException
     */
    public Map<Chain, PDBUniProtMatch> mapCoordinatesToUniProtInPDB(Structure structure,
                                                                    String[] uniprotIDs,
                                                                    Map<String, Sequence> uniprotToSeq,
                                                                    Map<String, String> uniprotToGene) throws IOException {
        List<Chain> chains = structure.getChains();
        Iterator<Chain> chainIterator = chains.iterator();
        List<String> uniprotIDsList = new ArrayList<>(Arrays.asList(uniprotIDs));
        Iterator<String> uniprotIDiterator = uniprotIDsList.iterator();
        Map<Chain, PDBUniProtMatch> chainToMatch = new HashMap<>();
        while (chainIterator.hasNext()) {
            Chain chain = chainIterator.next();
            // Don't use printLines as a parameter in the public API method!
//            if(printLines) {
//                System.out.println(chain.getChainID() + ": " + chain.getAtomSequence());
//            }
            while (uniprotIDiterator.hasNext()) {
                String uniprotID = uniprotIDiterator.next();
                Sequence seq = uniprotToSeq.get("UniProt:" + uniprotID);
//                System.out.println(acc + ": " + seq.getSequence());
                int index = seq.getSequence().indexOf(chain.getAtomSequence());

                if (index >= 0) {
                    //try {
                    //    calculateGlobalAlignmentScore(chainID.getAtomSequence(), seq.getSequence());
                    //}catch(CompoundNotFoundException cnfe){
                        //do nothing
                    //}
//                    System.out.println(seq.getSequence());
//                    System.out.println(chainID.getAtomSequence());
                    PDBUniProtMatch match = new PDBUniProtMatch();
                    match.setChainID(chain.getChainID());
                    match.setChainSequence(chain.getAtomSequence());
                    match.setGene(uniprotToGene.get(uniprotID));
                    match.setUniprot(uniprotID);
                    match.setPdbStart(chain.getAtomGroup(0).getResidueNumber().getSeqNum());
                    match.setUniprotStart(index + 1);
                    match.setOffset(match.getUniprotStart() - match.getPdbStart());
                    chainToMatch.put(chain, match);

                    //no duplicates
                    uniprotIDiterator.remove();
                    break;
                }
            }
        }
        return chainToMatch;
    }

    private double calculateGlobalAlignmentScore(String s1, String s2) throws CompoundNotFoundException {
        List<ProteinSequence> psl = new ArrayList<>();
        psl.add(new ProteinSequence(s1));
        psl.add(new ProteinSequence(s2));
        SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getBlosum62();
        double[] pair =
                Alignments.getAllPairsScores(
                        psl,
                        Alignments.PairwiseSequenceScorerType.GLOBAL,
                        new SimpleGapPenalty(),
                        matrix);
        for(int i = 0; i < pair.length; i++) {
            System.out.println(String.format("%s vs %s: %f", s1, s2, pair[i]));
        }
        return -1.0;
    }

    private class InteractionPDBScore implements Comparable<InteractionPDBScore> {
        
        private boolean isInRepresent;
        private String fileName;
        private int interactResidues;
        private String fi;
        private int rank;
        
        @Override
        public int compareTo(InteractionPDBScore o) {
            if (interactResidues != o.interactResidues)
                return o.interactResidues - interactResidues;
            return o.rank - rank;
        }
        
        public String getShortFileName() {
            int index = fileName.lastIndexOf(File.separator);
            if (index >= 0)
                return fileName.substring(index + 1);
            return fileName;
        }
        
    }
    
}
