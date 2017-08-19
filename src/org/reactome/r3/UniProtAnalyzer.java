/*
 * Created on Jun 20, 2006
 *
 */
package org.reactome.r3;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.gk.database.util.ReferencePeptideSequenceAutoFiller;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.funcInt.Protein;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

public class UniProtAnalyzer {
    //public static final String UNI_DIR_NAME = "D:\\documents\\Stein_lab\\Reactome\\Data\\";
    //private final String UNI_SPROT_FILE_NAME =  UNI_DIR_NAME + "uniprot_sprot_human.dat";
    //private final String UNI_TREMBL_FILE_NAME = UNI_DIR_NAME + "uniprot_trembl_human.dat";
    private final String UNI_SPROT_FILE_NAME = R3Constants.UNIPROT_DIR + "uniprot_sprot_human.dat";
    private final String UNI_TREMBL_FILE_NAME = R3Constants.UNIPROT_DIR + "uniprot_trembl_human.dat";
    private MySQLAdaptor reactomeDba;
    private FileUtility fu = new FileUtility();
    
    public String getDescription(String uniProtId) throws Exception {
        int index = uniProtId.indexOf("-");
        if (index > 0)
            uniProtId = uniProtId.substring(0, index);
        if (!(uniProtId.startsWith("P") ||
              uniProtId.startsWith("Q") ||
              uniProtId.startsWith("O"))) {
            System.out.println(uniProtId + " is not UniProt ID");
            return null;
        }
        if (reactomeDba == null) {
            reactomeDba = new MySQLAdaptor("localhost",
                                           "gk_central_101606",
                                           "root",
                                           "macmysql01",
                                           3306);
        }
        Collection c = reactomeDba.fetchInstanceByAttribute(ReactomeJavaConstants.ReferencePeptideSequence, 
                                                            ReactomeJavaConstants.identifier, 
                                                            "=", 
                                                            uniProtId);
        GKInstance refPepSeq = null;
        if (c != null && c.size() > 0) 
            refPepSeq = (GKInstance) c.iterator().next();
        else {
            refPepSeq = new GKInstance();
            refPepSeq.setSchemaClass(reactomeDba.getSchema().getClassByName(ReactomeJavaConstants.ReferencePeptideSequence));
            refPepSeq.setDbAdaptor(reactomeDba);
            ReferencePeptideSequenceAutoFiller autoFiller = new ReferencePeptideSequenceAutoFiller();
            autoFiller.setPersistenceAdaptor(reactomeDba);
            refPepSeq.setAttributeValue(ReactomeJavaConstants.identifier, uniProtId);
            autoFiller.process(refPepSeq);
        }
        String desc = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.description);
        if (desc == null || desc.length() == 0) {
            desc = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.name);
        }
        return desc;
    }
    
    public Map<String, String> loadUniProtToRefSeqMap() throws IOException {
        String fileName = "datasets" + File.separator + "iproclass" + File.separator + "human.iproclass.tb";
        FileUtility fu = new FileUtility();
        String line = null;
        fu.setInput(fileName);
        Map<String, String> uniProt2RefSeq = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            uniProt2RefSeq.put(tokens[0],
                               tokens[3]);
        }
        fu.close();
        return uniProt2RefSeq;
    }
    
    @Test
    public void generateENSEBMLToUniMap() throws IOException {
        long time1 = System.currentTimeMillis();
        // Only SwissProt part is handled for the pilot project
        FileUtility fu = new FileUtility();
        fu.setInput(UNI_SPROT_FILE_NAME);
        // For output
        String outputFileName = R3Constants.RESULT_DIR + "ensembl2uni.txt";
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outputFileName);
        String line = null;
        String ac = null;
        int index = 0;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("AC")) {
                index = line.indexOf(";");
                ac = line.substring(5);
            }
            else if (line.startsWith("DR   Ensembl;")) {
                // Escape the first ";"
                index = line.indexOf(";", 13);
                String ensemlId = line.substring(13, index).trim();
                String[] acs = ac.split("(;|; )");
                for (String a : acs) {
                    outFu.printLine(ensemlId + "\t" + a.trim());
                }
            }
        }
        // Extract TREMBL
        fu.setInput(UNI_TREMBL_FILE_NAME);
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("AC")) {
                index = line.indexOf(";");
                ac = line.substring(5);
            }
            else if (line.startsWith("DR   Ensembl;")) {
                // Escape the first ";"
                index = line.indexOf(";", 13);
                String ensemlId = line.substring(13, index).trim();
                String[] acs = ac.split("(;|; )");
                for (String a : acs) {
                    outFu.printLine(ensemlId + "\t" + a.trim());
                }
            }
        }
        // Need to close all writers and readers.
        fu.close();
        outFu.close();
        long time2 = System.currentTimeMillis();
        System.out.println("Time for generateEMBLToUniMap: " + (time2 - time1));
    }
    
    /**
     * Generate a map file from EMBL ids to UniProt accession numbers.
     * @throws IOException
     */
    public void generateEMBLToUniMap() throws IOException {
        long time1 = System.currentTimeMillis();
        // Only SwissProt part is handled for the pilot project
        String uniFileName = "/Users/wgm/Documents/caBIG_R3/datasets/UniProt/uniprot_sprot_human.dat";
        FileReader fileReader = new FileReader(uniFileName);
        BufferedReader reader = new BufferedReader(fileReader);
        // For output
        String outputFileName = "results/embl2uni.txt";
        FileWriter fileWriter = new FileWriter(outputFileName);
        PrintWriter printWriter = new PrintWriter(fileWriter);
        String line = null;
        String ac = null;
        String emblString = null;
        String[] emblIds = null;
        int index = 0;
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("AC")) {
                index = line.indexOf(";");
                ac = line.substring(5, index);
            }
            else if (line.startsWith("DR   EMBL")) {
                emblString = line.substring(11); // Escape "; " after EMBL
                emblIds = emblString.split("; ");
                // Sometimes there is no "-"
                int c = emblIds.length;
                int lastC = 0;
                // Have to find "-" since more tokens might exist
                for (int i = 0; i < c; i++) {
                    if (emblIds[i].equals("-") ||
                        emblIds[i].equals("JOINED")) {
                        lastC = i;
                        break;
                    }
                }
                if (lastC == 0)
                    lastC = c - 1;
                for (int i = 0; i < lastC; i++) {
                    printWriter.println(emblIds[i] + "\t" + ac);
                }
            }
        }
        // Need to close all writers and readers.
        reader.close();
        fileReader.close();
        printWriter.close();
        fileWriter.close();
        long time2 = System.currentTimeMillis();
        System.out.println("Time for generateEMBLToUniMap: " + (time2 - time1));
    }
    
    /**
     * Several UniProt Accession identifiers can be the same. This method is used 
     * to generate a map for several accession identifiers to the first AC id.
     * @throws IOException
     */
    public void generateUniProtIDsMap() throws IOException {
        // Run the following twice
        // Want to get both Swiss and Trembl files for all UniProt accession numbers
        String files[] = new String[] {
                UNI_SPROT_FILE_NAME,
                UNI_TREMBL_FILE_NAME
        };
        Map<String, String> idMaps = processUniProtIds(Arrays.asList(files));
        FileUtility fu = new FileUtility();
        String mapFileName = R3Constants.UNIPROT_DIR + "ACIDMap.txt";
        fu.exportMap(idMaps, mapFileName);
        
        // Want to use SwingProt only
        files = new String[] {
                UNI_SPROT_FILE_NAME,
        };
        idMaps = processUniProtIds(Arrays.asList(files));
        mapFileName = R3Constants.UNIPROT_DIR + "SwissProtACIDMap.txt";
        fu.exportMap(idMaps, mapFileName);
    }
    
    public void countUniProtEntries() throws IOException {
        FileUtility fu = new FileUtility();
        String dirName = R3Constants.UNIPROT_DIR;
        Map<String, String> map = fu.importMap(dirName + "SwissProtACIDMap.txt");
        Set<String> entries = new HashSet<String>(map.values());
        System.out.println("SwissProt Entries: " + entries.size());
        // Output file
        fu.outputSet(entries, dirName + "SwissProtIDs.txt");
        map = fu.importMap(dirName + "ACIDMap.txt");
        entries = new HashSet<String>(map.values());
        System.out.println("SwissProt and Trembl Entries: " + entries.size());
    }
    
    public Set<String> loadSwissProtIds() throws IOException {
        FileUtility fu = new FileUtility();
        return fu.loadSet(R3Constants.UNIPROT_DIR + "SwissProtIDs.txt");
    }
    
    /**
     * Load the pre-generated UniProt Accession identifier maps.
     * @return
     * @throws IOException
     * @see generateUniProtIDsMap()
     */
    public Map<String, String> loadUniProtIDsMap() throws IOException { 
        FileUtility fu = new FileUtility();
        String mapFileName = R3Constants.UNIPROT_DIR + "ACIDMap.txt";
        Map<String, String> acIdMap = fu.importMap(mapFileName);
        return acIdMap;
    }
    
    public boolean isHumanID(Set<String> uniProtId,
                             String id) {
        if (id.matches("^[0-9]+")) // id from NCBI or HPRD. No problems.
            return true;
        // id from IntAct (some non-human Ids) or from Reactome: some HIV
        // Have to make sure only accesion number not splice forms are used
        int index = id.indexOf("-");
        if (index > 0)
            id = id.substring(0, index);
        return uniProtId.contains(id);
    }
    
    /**
     * Load the pre-generated SwissProt Accession identifier maps.
     * @return
     * @throws IOException
     * @see generateUniProtIDsMap()
     */
    public Map<String, String> loadSwissProtIDsMap() throws IOException { 
        FileUtility fu = new FileUtility();
        //String mapFileName = "/Users/wgm/Documents/caBIG_R3/datasets/UniProt/SwissProtACIDMap.txt";
        String mapFileName = R3Constants.UNIPROT_DIR + "SwissProtACIDMap.txt";
        Map<String, String> acIdMap = fu.importMap(mapFileName);
        return acIdMap;
    }
    
    private Map<String, String> processUniProtIds(List<String> files) throws IOException {
        FileUtility fu = new FileUtility();
        Map<String, String> uniAccsToIdsMap = new HashMap<String, String>();
        List<String> acInEntry = new ArrayList<String>();
        for (String fileName : files) {
            fu.setInput(fileName);
            String line = null;
            String ids = null;
            while ((line = fu.readLine()) != null) {
                if (line.startsWith("AC")) {
                    ids = line.substring(2);
                    parseIds(ids, acInEntry);
                }
                else if (line.startsWith("//")) {
                    // Change to a new entry
                    String value = acInEntry.get(0);
                    for (String ac : acInEntry)
                        uniAccsToIdsMap.put(ac, value);
                    acInEntry.clear();
                }
            }
            fu.close();
        }
        System.out.println("total entries in uniprot: " + uniAccsToIdsMap.size());
        return uniAccsToIdsMap;
    }
    
    private void parseIds(String ids, List<String> acInEntry) {
        ids = ids.trim();
        String[] tmp = ids.split("(;| )+");
        for (String id : tmp) {
            acInEntry.add(id);
        }
    } 
    
    public Map<String, String> getUniProtAccessionToGeneName() throws IOException {
        Map<String, Set<String>> geneToAcc = generateGeneNameToUniAccesssMap(true);
        Map<String, String> accToGene = new HashMap<String, String>();
        for (String gene : geneToAcc.keySet()) {
            Set<String> acces = geneToAcc.get(gene);
            for (String acc : acces)
                accToGene.put(acc, gene);
        }
        return accToGene;
    }
    
    /**
     * Generate a map from gene name to UniProt accession numbers
     * @param needSwissProtOnly
     * @return
     * @throws IOException
     */
    public Map<String, Set<String>> generateGeneNameToUniAccesssMap(boolean needSwissProtOnly) throws IOException {
        String[] fileNames = new String[] {
                UNI_SPROT_FILE_NAME,
                UNI_TREMBL_FILE_NAME
        };
        String line = null;
        int index1, index2;
        String acc = null;
        String geneName;
        Map<String, Set<String>> nameToAcces = new HashMap<String, Set<String>>();
        for (String fileName : fileNames) {
            File file = new File(fileName);
            if (!file.exists()) {
                System.err.println(file.getAbsolutePath() + " doesn't exist!");
                continue;
            }
            fu.setInput(fileName);
            while ((line = fu.readLine()) != null) {
                if (acc == null && line.startsWith("AC")) { // Want to have the first accession. Acc lines many span more than one
                    index1 = line.indexOf(" ");
                    index2 = line.indexOf(";"); // Want the first accession number only
                    acc = line.substring(index1 + 1, index2).trim();
                }
                else if (line.startsWith("GN   Name=")) { // Sometimes there are other GN lines (e.g. P02768)
                    // Get gene name
                    geneName = extractGeneName(line);
                    InteractionUtilities.addElementToSet(nameToAcces, geneName, acc);
//                    // Check Synonyms
//                    index1 = line.indexOf("Synonyms=");
//                    if (index1 > 0) {
//                        index2 = line.indexOf(";", index1);
//                        if (index2 > 0) { // Some synonyms may be broken to more than one line. Just ignore them. This occurs in the trembl file only!
//                            String tmp = line.substring(index1 + "Synonyms=".length(), index2);
//                            String[] tokens = tmp.split(", ");
//                            for (String token : tokens) {
//                                InteractionUtilities.addElementToSet(nameToAcces, token, acc);
//                            }
//                        }
//                    }
                }
                else if (line.startsWith("//")) {
                    //rest
                    acc = null;
                }
            }
            fu.close();
            if (needSwissProtOnly)
                break;
        }
        return nameToAcces;
    }
    
    /**
     * Load the map from gene to the maximum length of its protein products. A gene may have multiple
     * protein products having different lengths. This method returns the longest length. If you want
     * to get the map from a gene to its all protein lengths, use loadGeneToProteinLengths().
     * @return
     * @throws IOException
     */
    public Map<String, Integer> loadGeneToProteinLength() throws IOException {
        Map<String, List<Integer>> geneToLengths = loadGeneToProteinLengths();
        Map<String, Integer> geneToLength = new HashMap<String, Integer>();
        for (String gene : geneToLengths.keySet()) {
            List<Integer> list = geneToLengths.get(gene);
            if (list.size() == 1)
                geneToLength.put(gene, list.get(0));
            else {
                // Find the longest
                Integer length = Integer.MIN_VALUE;
                for (Integer value : list) {
                    if (value > length)
                        length = value;
                }
                geneToLength.put(gene, length);
            }
        }
        return geneToLength;
    }
    
    public Map<String, List<Integer>> loadGeneToProteinLengths() throws IOException {
        String fileName = UNI_SPROT_FILE_NAME;
        Map<String, List<Integer>> geneToLengths = new HashMap<String, List<Integer>>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        Set<String> genes = new HashSet<String>();
        int index1, index2;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("GN   Name=")) { // Gene name
                String gene = extractGeneName(line);
                genes.add(gene);
                // A hard-coded fix
                if (gene.equals("STK26"))
                    genes.add("MST4");
                // Check Synonyms
                index1 = line.indexOf("Synonyms=");
                if (index1 > 0) {
                    index2 = line.indexOf(";", index1);
                    if (index2 > 0) { // Some synonyms may be broken to more than one line. Just ignore them. This occurs in the trembl file only!
                        String tmp = line.substring(index1 + "Synonyms=".length(), index2);
                        String[] tokens = tmp.split(", ");
                        for (String token : tokens) {
                            genes.add(token);
                        }
                    }
                }
            }
            else if (line.startsWith("SQ   SEQUENCE")) { // Sequence length
                if (genes.size() > 0) {
                    index1 = line.indexOf("AA");
                    String text = line.substring("SQ   SEQUENCE".length() + 1, index1).trim();
                    Integer length = new Integer(text);
                    for (String gene : genes) {
                        List<Integer> list = geneToLengths.get(gene);
                        if (list == null) {
                            list = new ArrayList<Integer>();
                            geneToLengths.put(gene, list);
                        }
                        list.add(length);
                    }
                }
            }
            if (line.startsWith("//"))
                genes.clear();
        }
        return geneToLengths;
    }
    
    @Test
    public void testLoadGeneToProteinLengths() throws IOException {
        Map<String, List<Integer>> geneToLengths = loadGeneToProteinLengths();
        int count = 0;
        for (String gene : geneToLengths.keySet()) {
            List<Integer> list = geneToLengths.get(gene);
            System.out.println(gene + "\t" + list);
            if (list.size() > 1)
                count ++;
        }
        System.out.println("Genes having more than one length: " + count);
    }
    
    private String extractGeneName(String line) {
        int index1;
        int index2;
        String geneName;
        // Get gene name
        index1 = line.indexOf("Name=");
        index2 = line.indexOf(";");
        // As of December, 2014, there are some lines as following:
//                    GN   Name=QRSL1 {ECO:0000255|HAMAP-Rule:MF_03150};
//                    GN   Name=GATB {ECO:0000255|HAMAP-Rule:MF_03147,
        if (index2 < 0) // For case: GN   Name=GATB {ECO:0000255|HAMAP-Rule:MF_03147,
            index2 = line.indexOf(" ", index1);
        geneName = line.substring(index1 + "Name=".length(), index2);
        if (geneName.contains("{")) { // For case: Name=QRSL1 {ECO:0000255|HAMAP-Rule:MF_03150}
            index2 = geneName.indexOf("{");
            geneName = geneName.substring(0, index2).trim();
        }
        return geneName;
    }
    
    public Map<String, Protein> generateUniAccToProteinMap() throws IOException {
        String[] fileNames = new String[] {
                UNI_SPROT_FILE_NAME,
                UNI_TREMBL_FILE_NAME
        };
        FileUtility fu = new FileUtility();
        int index = 0, index1 = 0;
        String ac = null, de = null, gn = null;
        boolean inNewEntry = true;
        Map<String, Protein> accToProtein = new HashMap<String, Protein>();
        for (String fileName : fileNames) {
            fu.setInput(fileName);
            String line = null;
            while ((line = fu.readLine()) != null) {
                if (inNewEntry && line.startsWith("AC")) {
                    inNewEntry = false;
                    if (ac != null) {
                        //System.out.println(ac + "\t" + de + "\t" + gn);
                        Protein protein = new Protein();
                        protein.setPrimaryDbName("UniProt");
                        protein.setPrimaryAccession(ac);
                        protein.setName(de);
                        protein.setShortName(gn);
                        accToProtein.put(ac, protein);
                    }
                    index = line.indexOf(";");
                    ac = line.substring(2, index).trim();
                    de = gn = null;
                }
                else if (line.startsWith("DE") && (de == null)) { // Want to check only one line
//                    index = line.indexOf("(");
//                    if (index > 0)
//                        de = line.substring(2, index).trim();
//                    else
//                        de = line.substring(2).trim();
//                    // Get rid of the last period
//                    if (de.endsWith("."))
//                        de = de.substring(0, de.length() - 1);
                    // For release 14.9
                    index = line.indexOf("Full=");
                    index1 = line.indexOf(";", index);
                    de = line.substring(index + "Full=".length(), index1);
                }
                else if (line.startsWith("GN") && (gn == null)) {
                    index = line.indexOf("Name=");
                    if (index > 0) { // Otherwise, no gene name
                        index1 = line.indexOf(";");
                        gn = line.substring(index + 5, index1);
                    }
                }
                else if (line.startsWith("//")) {
                    inNewEntry = true;
                }
            }
            fu.close();
        }
        return accToProtein;
    }
    
    @Test
    public void generateUniToPfamMap() throws IOException {
        String[] fileNames = new String[] {
                UNI_SPROT_FILE_NAME,
                UNI_TREMBL_FILE_NAME
        };
        Map<String, Set<String>> uni2Pfam = new HashMap<String, Set<String>>();
        processUni2Pfam(fileNames, uni2Pfam);
        FileUtility fu = new FileUtility();
        String outputFileName = R3Constants.UNIPROT_DIR + "Uni2Pfam.txt";
        fu.saveSetMap(uni2Pfam, outputFileName);
    }
    
    private void processUni2Pfam(String[] srcFileNames,
                                 Map<String, Set<String>> uni2PfamMap) throws IOException {
        String ac = null;
        String pfam = null;
        String line = null;
        Set<String> pFamSet = null;
        int index = 0;
        FileUtility fu = new FileUtility();
        for (String fileName : srcFileNames) {
            fu.setInput(fileName);
            while ((line = fu.readLine()) != null) {
                if (line.startsWith("AC")) {
                    index = line.indexOf(";");
                    ac = line.substring(5, index);
                    pFamSet = null; // reinitialization
                }
                else if (line.startsWith("DR   Pfam;")) {
                    // Starting from 11
                    index = line.indexOf(";", 11);
                    pfam = line.substring(11, index);
                    if (pFamSet == null)
                        pFamSet = new HashSet<String>();
                    pFamSet.add(pfam);
                }
                else if (line.startsWith("//")) {
                    // Done with one entry
                    if (pFamSet != null)
                        uni2PfamMap.put(ac, pFamSet);
                }
            }
            fu.close();
        }
    }
    
    public void generateEntrezGeneToUniProt() throws IOException {
        //String outputFile = "/Users/guanming/Documents/caBIG_R3/datasets/iproclass/EntrezToUniProt.txt";
        String inputFile = "/Users/gwu/Documents/caBIG_R3/datasets/iproclass/011612/h_sapiens.tb";
        String outputFile = R3Constants.ENTREZ_TO_UNIPROT_MAP_FILE_NAME;
        FileUtility fu = new FileUtility();
        fu.setInput(inputFile);
        String line = null;
        int c = 0;
        Set<String> outPutSet = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens[2].length() == 0)
                continue;
            c++;
            // Multiple Entrez Gene Ids may be used
            String[] entrezTokens = tokens[2].split("; ");
            for (String entrezToken : entrezTokens) {
                outPutSet.add(entrezToken + "\t" + tokens[0]);
            }
        }
        fu.saveInteractions(outPutSet, outputFile);
        fu.close();
        System.out.println("Counter: " + c);
    }
    
    private Map<String, Set<String>> loadEntrezIdToUniProts() throws IOException {
        String fileName = "/Users/wgm/Documents/caBIG_R3/datasets/iproclass/human.EntrezToUniProt.txt";
        FileUtility fu = new FileUtility();
        String line = null;
        Map<String, Set<String>> map = new HashMap<String, Set<String>>();
        int index = 0;
        fu.setInput(fileName);
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            String part1 = line.substring(0, index);
            String uniProt = line.substring(index + 1);
            String[] tokens = part1.split("[; ]+");
            for (String tmp : tokens) {
                Set<String> set = map.get(tmp);
                if (set == null) {
                    set = new HashSet<String>();
                    map.put(tmp, set);
                }
                set.add(uniProt);
            }
        }
        return map;
    }
    
    public Map<String, String> mapEntrezGeneIdToUniProt(Set<String> geneIds) throws IOException {
        Map<String, String> geneIdToUniProt = new HashMap<String, String>();
        Map<String, Set<String>> wholeMap = loadEntrezIdToUniProts();
        // Try SwissProt first. If all passed gene ids can be mapped to SwissProt,
        // don't try to use the whole set of UniProt to avoid overwriting swissprot
        // mapping
        Map<String, String> swissProtIdMap = loadSwissProtIDsMap();
        Set<String> left = new HashSet<String>();
        boolean isFound = false;
        for (String geneId : geneIds) {
            Set<String> set = wholeMap.get(geneId);
            if (set != null) {
                // Find swissProt 
                isFound = false;
                for (String tmp : set) {
                    String swissId = swissProtIdMap.get(tmp);
                    if (swissId != null) {
                        geneIdToUniProt.put(geneId, swissId);
                        isFound = true;
                        break;
                    }
                }
                if (!isFound)
                    left.add(geneId);
            }
        }
        if (left.size() == 0)
            return geneIdToUniProt;
        Map<String, String> uniprotIdMap = loadUniProtIDsMap();
        for (String geneId : left) {
            Set<String> set = wholeMap.get(geneId);
            if (set != null) {
                // Find swissProt 
                isFound = false;
                for (String tmp : set) {
                    String uniProtId = uniprotIdMap.get(tmp);
                    if (uniProtId != null) {
                        geneIdToUniProt.put(geneId, uniProtId);
                        isFound = true;
                        break;
                    }
                }
                if (!isFound)
                    left.add(geneId);
            }
        }
        return geneIdToUniProt;
    }
    
    /**
     * One gene name can be mapped to one uniprot only.
     * @return
     * @throws IOException
     * @deprecated @see UCSCDataAnalyzer.mapToUniProt(Collection<String>)
     */
    public Map<String, String> mapGeneNamesToUniProt(Collection<String> geneNames) throws IOException {
        String mapFile = "/Users/wgm/Documents/caBIG_R3/datasets/HPRD/GeneSymbolToUniProtFromHGNC.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(mapFile);
        Map<String, String> map = new HashMap<String, String>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 4) {
                //System.out.println("Line cannot be used: " + line);
                continue;
            }
            // tokens[1] gene name, tokens[3] uniprot id
            map.put(tokens[1], tokens[3]);
        }
        // Check if any dulpicate is found
        Map<String, String> results = new HashMap<String, String>();
        for (String name : geneNames) {
            String uni = map.get(name);
            if (uni != null)
                results.put(name, uni);
        }
        return results;
    }
    
    public Map<String, Set<String>> loadMimToUniProt() throws IOException {
        Map<String, Set<String>> mimToSp = new HashMap<String, Set<String>>();
        String fileName = R3Constants.UNIPROT_DIR + "mimtosp.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        Pattern mimPattern = Pattern.compile("^(\\d+):");
        Pattern spIdPattern = Pattern.compile("\\((.{6})\\)");
        int start = 0;
        String mimId = null;
        Set<String> spIds = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#"))
                continue;
            //System.out.println(line);
            Matcher matcher = mimPattern.matcher(line);
            if (matcher.find()) {
                mimId = matcher.group(1);
                spIds = new HashSet<String>();
                mimToSp.put(mimId, spIds);
                extractSpIds(line, spIdPattern, spIds);
            }
            else 
                extractSpIds(line, 
                             spIdPattern,
                             spIds);
        }
        fu.close();
        return mimToSp;
    }
    
    private void extractSpIds(String line,
                              Pattern pattern,
                              Set<String> ids) {
        Matcher matcher = pattern.matcher(line);
        int start = 0;
        while (matcher.find(start)) {
            String spId = matcher.group(1);
            ids.add(spId);
            start += spId.length();
        }
    }
    
    public void testLoadMimToUniProt() throws Exception {
        Map<String, Set<String>> mimToSp = loadMimToUniProt();
        System.out.println("Total MIM: " + mimToSp.size());
        // For some tests
        String[] mimIds = new String[] {
                "108770",
                "109100",
                "109195"
        };
        for (String mimId : mimIds) {
            Set<String> spIds = mimToSp.get(mimId);
            System.out.printf("%s: %s%n",
                              mimId,
                              spIds);
        }
    }
    
    @Deprecated
    public Map<String, String> loadGeneNameToUniProt() throws IOException {
        String mapFile = "/Users/wgm/Documents/caBIG_R3/datasets/HPRD/GeneSymbolToUniProtFromHGNC.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(mapFile);
        Map<String, String> map = new HashMap<String, String>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 4) {
                //System.out.println("Line cannot be used: " + line);
                continue;
            }
            // tokens[1] gene name, tokens[3] uniprot id
            map.put(tokens[1], tokens[3]);
        }
        return map;
    }
    
    /**
     * This method has been deprecated. Use method getUniProtAccessionToGeneName() instead.
     * @return
     * @throws IOException
     */
    @Deprecated
    public Map<String, String> loadUniProtToGeneNameMap() throws IOException {
        String mapFile = "/Users/wgm/Documents/caBIG_R3/datasets/HPRD/GeneSymbolToUniProtFromHGNC.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(mapFile);
        Map<String, String> map = new HashMap<String, String>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 4) {
                //System.out.println("Line cannot be used: " + line);
                continue;
            }
            // tokens[1] gene name, tokens[3] uniprot id
            map.put(tokens[3], tokens[1]);
        }
        return map;
    }
    
    public void filterIdMappingToHuman() throws IOException {
        // The file will be deleted after this method running to save space
        String inFileName = R3Constants.UNIPROT_DIR + "idmapping_selected.tab";
        fu.setInput(inFileName);
        FileUtility outFu = new FileUtility();
        String outFileName = R3Constants.UNIPROT_DIR + "human_idmapping_selected.tab";
        outFu.setOutput(outFileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // human ids
            if (tokens[1].endsWith("_HUMAN")) {
                outFu.printLine(line);
            }
        }
        fu.close();
        outFu.close();
    }
    
}
