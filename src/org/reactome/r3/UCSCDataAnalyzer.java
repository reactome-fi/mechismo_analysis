/*
 * Created on Jun 1, 2007
 *
 */
package org.reactome.r3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math.stat.inference.TestUtils;
import org.junit.Test;
import org.reactome.genome.Exon;
import org.reactome.genome.RefSeqInfo;
import org.reactome.genome.Segment;
import org.reactome.genome.Transcript;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

import weka.core.Utils;

public class UCSCDataAnalyzer {
    
    private final String DATA_DIR = "datasets/ucsc/080113/";
    // cache these two maps to increase performance
    private Map<String, String> geneToUniMap;
    private Map<String, String> geneToHprdMap;
    private FileUtility fu = new FileUtility();
    
    public UCSCDataAnalyzer() {
    }
    
    private Map<String, Double> calculateConservationScores(List<String> hubs,
                                                            List<String> nonHubs,
                                                            boolean isForProteins)
            throws Exception, MathException {
        Map<String, Double> hubScores = calculateConservationScores(hubs, isForProteins);
        Map<String, Double> nonHubScores = calculateConservationScores(nonHubs, isForProteins);
        System.out.println("Summary for hubs: " + hubs.size());
        summarize(hubScores);
        System.out.println("\nSummary for non-hubs: " + nonHubs.size());
        summarize(nonHubScores);
        // Do a ttest
        double[] hubArray = new double[hubScores.size()];
        int index = 0;
        for (String protein : hubScores.keySet()) {
            Double score = hubScores.get(protein);
            hubArray[index] = score;
            index ++;
        }
        index = 0;
        double[] nonHubArray = new double[nonHubScores.size()];
        for (String protein : nonHubScores.keySet()) {
            nonHubArray[index] = nonHubScores.get(protein);
            index ++;
        }
        double pvalue = TestUtils.tTest(hubArray, nonHubArray);
        System.out.println("\npvalue: " + pvalue);
        Map<String, Double> allScore = new HashMap<String, Double>();
        allScore.putAll(hubScores);
        allScore.putAll(nonHubScores);
        return allScore;
    }
    
    private void calculateCorrelationCoefficient(Map<String, Double> proteinToScore,
                                                 Map<String, Integer> proteinToDegree) {
        double[] scores = new double[proteinToScore.size()];
        double[] degrees = new double[proteinToScore.size()];
        int index = 0;
        for (String protein : proteinToScore.keySet()) {
            Double score = proteinToScore.get(protein);
            Integer degree = proteinToDegree.get(protein);
            scores[index] = score;
            degrees[index] = degree;
            index ++;
        }
        double cc = Utils.correlation(scores, degrees, scores.length);
        System.out.println("Correlation coefficient: " + cc);
    }
    
    private void summarize(Map<String, Double> nodeToScore) {
        DescriptiveStatistics stat = new DescriptiveStatistics();
        for (Double score : nodeToScore.values())
            stat.addValue(score);
        Double mean = stat.getMean();
        Double max = stat.getMax();
        Double min = stat.getMin();
        Double stderr = stat.getStandardDeviation();
        System.out.println("Mean: " + mean);
        System.out.println("std: " + stderr);
        System.out.println("Max: " + max);
        System.out.println("Min: " + min);
        System.out.println("Total data points: " + nodeToScore.size());
    }
    
    
    /**
     * Load gene names to coordinates.
     * @return String (gene name) -> List<String> (Chr, txStart, txEnd)
     * @throws IOException
     */
    public Map<String, List<String[]>> loadGeneNameToCoordinates() throws IOException {
        String fileName = R3Constants.UCSC_DIR + "HumanRefSeqGenes.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, List<String[]>> nameToCoordiantes = new HashMap<String, List<String[]>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String geneName = tokens[12];
            List<String[]> list = nameToCoordiantes.get(geneName);
            if (list == null) {
                list = new ArrayList<String[]>();
                nameToCoordiantes.put(geneName, list);
            }
            String[] coordinate = new String[3];
            coordinate[0] = tokens[2];
            coordinate[1] = tokens[4];
            coordinate[2] = tokens[5];
            list.add(coordinate);
        }
        fu.close();
        return nameToCoordiantes;
    }
    
    public Map<String, Integer> loadGeneNameToAverageExonLength() throws IOException {
        Map<String, List<Integer>> geneToLengths = loadGeneNameToTotalExonLengths();
        Map<String, Integer> rtn = new HashMap<String, Integer>();
        for (String gene : geneToLengths.keySet()) {
            List<Integer> list = geneToLengths.get(gene);
            if (list.size() == 1)
                rtn.put(gene, list.get(0));
            else {
                int total = 0;
                for (Integer i : list)
                    total += i;
                rtn.put(gene, total / list.size());
            }
        }
        return rtn;
    }
    
    public Map<String, List<Integer>> loadGeneNameToTotalExonLengths() throws IOException {
        String fileName = R3Constants.UCSC_DIR + "HumanRefSeqGenes.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, List<Integer>> geneToExonLengths = new HashMap<String, List<Integer>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String geneName = tokens[12];
            List<Integer> list = geneToExonLengths.get(geneName);
            if (list == null) {
                list = new ArrayList<Integer>();
                geneToExonLengths.put(geneName, list);
            }
            String[] starts = tokens[9].split(",");
            String[] ends = tokens[10].split(",");
            int total = 0;
            for (int i = 0; i < starts.length; i++) {
                total += (Integer.parseInt(ends[i]) - Integer.parseInt(starts[i]) + 1);
            }
            list.add(total);
        }
        fu.close();
        return geneToExonLengths;
    }
    
    @Test
    public void testLoadGeneNameToTotalExonLegnths() throws IOException {
        Map<String, List<Integer>> geneToLengths = loadGeneNameToTotalExonLengths();
        for (String gene : geneToLengths.keySet()) {
            List<Integer> list = geneToLengths.get(gene);
            if (list.size() > 1)
                System.out.println(gene + "\t" + list);
        }
    }
    
    private Map<String, List<String[]>> loadProteinToCoordinates() throws IOException {
        String fileName = R3Constants.UCSC_DIR + "knownGenes.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, List<String[]>> proteinToCoordinates = new HashMap<String, List<String[]>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String protein = tokens[10];
            List<String[]> list = proteinToCoordinates.get(protein);
            if (list == null) {
                list = new ArrayList<String[]>();
                proteinToCoordinates.put(protein, list);
            }
            String[] coordinate = new String[3];
            coordinate[0] = tokens[1];
            coordinate[1] = tokens[3];
            coordinate[2] = tokens[4];
            list.add(coordinate);
        }
        fu.close();
        return proteinToCoordinates;
    }
    
    /**
     * This method will provide you a map that you can choose only an interval 
     * of the scores, so that the speed is expected to be improved.
     * chr bin start end score will be organized as:
     * 
     * chr -> List<Integer> (tsStart, tsEnd, and Score)
     * @return
     * @throws Exception
     */
    public Map<String, List<int[]>> loadIntervalConservationScore() throws Exception{
        Map<String, List<int[]>> chrToCoordinatesScores = new HashMap<String, List<int[]>>();
        FileUtility fu = new FileUtility();
        String fileName = R3Constants.UCSC_DIR + "MammalPhastConsElements28wayPlacMammal.txt";
        fu.setInput(fileName);
        String line = fu.readLine(); // There is an extra line in the file as the header
        while((line = fu.readLine())!=null) {
            String[] tokens = line.split("\t");
            String chr = tokens[1];
            List<int[]> coordinatesAndScores = chrToCoordinatesScores.get(chr);
            if (coordinatesAndScores == null) {
                coordinatesAndScores = new ArrayList<int[]>();
                chrToCoordinatesScores.put(chr, coordinatesAndScores);
            }
            int[] coordinateAndScore = new int[3];
            coordinateAndScore[0] = Integer.parseInt(tokens[2]);
            coordinateAndScore[1] = Integer.parseInt(tokens[3]);
            coordinateAndScore[2] = Integer.parseInt(tokens[5]);
            coordinatesAndScores.add(coordinateAndScore);
        }
        return chrToCoordinatesScores;
    }
    
    private Map<String, Set<String>> loadProteinToUCSCNames() throws IOException {
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> unIdMap = uniAnalyzer.loadUniProtIDsMap();
        String fileName = R3Constants.UCSC_DIR + "knownGenes.txt";
        fu.setInput(fileName);
        Map<String, Set<String>> uniToUCSCNames = new HashMap<String, Set<String>>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String name = tokens[0];
            String protein = tokens[10];
            String uniprot = unIdMap.get(protein);
            if (uniprot == null)
                continue;
            Set<String> set = uniToUCSCNames.get(uniprot);
            if (set == null) {
                set = new HashSet<String>();
                uniToUCSCNames.put(uniprot, set);
            }
            set.add(name);
        }
        fu.close();
        return uniToUCSCNames;
    }
    
    private Map<String, Double> calculateConservationScores(List<String> nodes,
                                                     boolean isForProteins) throws Exception {
        Map<String, List<int[]>> scoreMap = loadIntervalConservationScore();
        Map<String, List<String[]>> nameToCoordinates = null;
        if (isForProteins)
            nameToCoordinates = loadProteinToCoordinates();
        else
            nameToCoordinates = loadGeneNameToCoordinates();
        double singleScore;
        double geneScore;
        Map<String, Double> nodeToScore = new HashMap<String, Double>();
        int c = 0;
        for(String gene : nodes) {
            // Map proteins UniProt accession names to UCSC names
            List<String[]> coordinates = nameToCoordinates.get(gene);
            if (coordinates == null) {
                //System.err.println(protein + " cannot be found in known gene file!");
                c ++;
                continue;
            }
            geneScore = 0;
            // Start the score analysis
            for(String[] coordinate : coordinates) {
                singleScore = 0;
                int start = Integer.parseInt(coordinate[1]);
                int end = Integer.parseInt(coordinate[2]);
                int geneLength = end-start;
                String chr = coordinate[0];
                List<int[]> coordinatesAndScores = scoreMap.get(chr);
                if (coordinatesAndScores == null) {
                    System.err.println(chr + " has no scores!");
                    continue;
                }
                singleScore = calculateScore(singleScore, 
                                             start, 
                                             end,
                                             geneLength, 
                                             coordinatesAndScores);
                geneScore += singleScore;
            }
            //Final score is averaged according to the mapping
            //result.add(protein + "\t" + proteinScore / names.size());
            nodeToScore.put(gene, geneScore / coordinates.size());
        }
        System.out.println(c + " out of " + nodes.size() + " are not in the known gene list!");
        return nodeToScore;
    }
    
    private double calculateScore(double singleScore, int start, int end,
                                  int geneLength,
                                  List<int[]> coordinatesAndScores) {
        for (int[] coordinateAndScore : coordinatesAndScores) {
            int subStart = coordinateAndScore[0];
            int subEnd = coordinateAndScore[1];
            int subScore = coordinateAndScore[2];
            int d1 = subStart - start;
            int d2 = end - subEnd;
            
            if(d1 >= 0 && d2 >=0){
                singleScore += (double) subScore * (subEnd-subStart) / geneLength; 
            }
            else if(d1 >= 0 && d2 < 0 && end > subStart){
                singleScore += ((double) subScore * (end-subStart) / (subEnd-subStart)) * (subEnd-subStart) / geneLength;
            }
            else if(d1 < 0 && d2 >= 0 && start < subEnd){
                //subscore / sublength * (subEnd - start) / (end - start)
                singleScore += ((double)subScore * (subEnd-start) / (subEnd-subStart)) * (subEnd-subStart) / geneLength;
            }
            else if(d1 < 0 && d2 < 0){
                singleScore +=  ((double)(end-start) * subScore / (subEnd-subStart)) * (subEnd-subStart) / geneLength;
            }
        }
        return singleScore;
    }
    
    @Test
    public void checkUniProtIdsInKnownGenes() throws IOException {
        String fileName = DATA_DIR + "knownGenes.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> proteinIds = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            proteinIds.add(tokens[10]);
        }
        fu.close();
        System.out.println("Total protein ids: " + proteinIds.size());
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> acIdMap = uniAnalyzer.loadUniProtIDsMap();
        Set<String> mapped = new HashSet<String>();
        int notMapped = 0;
        for (String knownId : proteinIds) {
            // Want to remove splice isoform
            int index = knownId.indexOf("-");
            if (index > 0)
                knownId = knownId.substring(0, index);
            String uniId = acIdMap.get(knownId);
            if (uniId == null) {
                System.out.println("Cannot be mapped: " + knownId);
                notMapped ++;
                continue;
            }
            mapped.add(uniId);
        }
        System.out.println("After mapping: " + mapped.size());
        System.out.println("Not mapped: " + notMapped);
    }
    
    public List<RefSeqInfo> loadHumanRefSeqInfos() throws IOException {
        String fileName = DATA_DIR + "HumanRefSeqGenes.txt";
        return loadRefSeqInfo(fileName);
    }
    
    @Test
    public void testGetHumanGeneToTranscript() throws IOException {
        Map<String, Transcript> geneToTx = getHumanGeneToTranscript();
        System.out.println("Size of geneToTx: " + geneToTx.size());
        System.out.println("Transcript for IHH: " + geneToTx.get("IHH"));
        System.out.println("Transcript for TP53: " + geneToTx.get("TP53"));
        int count = 0;
        for (String gene : geneToTx.keySet()) {
            System.out.println(gene + "\t" + geneToTx.get(gene).getName());
            if (count++ == 50)
                break;
        }
    }
    
    public Map<String, Transcript> getHumanGeneToTranscript() throws IOException {
        String fileName = "datasets/ucsc/HumanRefSeqGenes_hg19_072816.txt";
        List<RefSeqInfo> infoList = loadTranscripts(fileName);
        Map<String, Transcript> geneToTranscript = new HashMap<String, Transcript>();
        for (RefSeqInfo info : infoList) {
            Transcript tx = info.getTranscript();
            Transcript oldTx = geneToTranscript.get(tx.getName());
            if (oldTx == null || oldTx.getCds().getLength() < tx.getCds().getLength())
                geneToTranscript.put(tx.getName(), tx);
        }
        return geneToTranscript;
    }
    
    /**
     * Grep a list of RefSeqInfo objects in the provided Segment from the 
     * list. 
     * @param segment
     * @param refSeqInfos
     * @param inclusive if true, both start and end coordinates should be in the
     * segment. false either of coordinates may be not in the segment, but not both.
     * @return
     */
    public List<RefSeqInfo> grepRefSeqInfo(Segment segment,
                                           Map<String, List<RefSeqInfo>> chrToRefSeq,
                                           boolean inclusive) {
        List<RefSeqInfo> refSeqs = new ArrayList<RefSeqInfo>();
        List<RefSeqInfo> refSeqInfos = chrToRefSeq.get(segment.getChromosome());
        int start = segment.getStart();
        int end = segment.getEnd();
        if (inclusive) {
            for (RefSeqInfo refSeq : refSeqInfos) {
                //if (refSeq.containedBy(segment)) 
                //    refSeqs.add(refSeq);
                // Inline the above checking to avoid some overhaul (esp. chromosome comparsion)
                if (segment.getStart() <= refSeq.getStart() && 
                    segment.getEnd() >= refSeq.getEnd())
                    refSeqs.add(refSeq);
            }
        }
        else {
            for (RefSeqInfo refSeq : refSeqInfos) {
                if (refSeq.overlap(segment))
                    refSeqs.add(refSeq);
            }
        }
        return refSeqs;
    }
    
    public Map<String, List<RefSeqInfo>> convertRefSeqListToMap(List<RefSeqInfo> list) {
        Map<String, List<RefSeqInfo>> map = new HashMap<String, List<RefSeqInfo>>();
        for (RefSeqInfo refSeq : list) {
            String chr = refSeq.getChromosome();
            List<RefSeqInfo> chrList = map.get(chr);
            if (chrList == null) {
                chrList = new ArrayList<RefSeqInfo>();
                map.put(chr, chrList);
            }
            chrList.add(refSeq);
        }
        return map;
    }
    
    public void sortRefSeqsInChr(Map<String, List<RefSeqInfo>> chrToRefSeqs) {
        Comparator<RefSeqInfo> sorter = new Comparator<RefSeqInfo>() {
            public int compare(RefSeqInfo info1, RefSeqInfo info2) {
                int start = info1.getStart();
                int end = info2.getEnd();
                return start - end;
            }
        };
        for (String chr : chrToRefSeqs.keySet()) {
            List<RefSeqInfo> list = chrToRefSeqs.get(chr);
            Collections.sort(list, sorter);
        }
    }
    
    public Map<String, List<String>> convertRefSeqToGeneInOrder(Map<String, List<RefSeqInfo>> chrToRefSeqs) {
        Map<String, List<String>> chrToGeneList = new HashMap<String, List<String>>();
        for (String chr : chrToRefSeqs.keySet()) {
            List<RefSeqInfo> refSeqs = chrToRefSeqs.get(chr);
            List<String> genes = new ArrayList<String>();
            // There are many repeat refseqs
            for (RefSeqInfo refseq : refSeqs) {
                if (!genes.contains(refseq.getGeneName())) {
                    genes.add(refseq.getGeneName());
                }
            }
            chrToGeneList.put(chr, genes);
        }
        return chrToGeneList;
    }
    
    public List<RefSeqInfo> loadTranscripts(String refSeqInfoFile) throws IOException {
        List<RefSeqInfo> list = new ArrayList<RefSeqInfo>();
        FileUtility fu = new FileUtility();
        fu.setInput(refSeqInfoFile);
        String line = fu.readLine(); // headers
        // for some checking
        int count = 0;
        int error = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            RefSeqInfo info = new RefSeqInfo();
            info.setName(tokens[1]);
            count ++;
//            System.out.println("\n" + count + ": " + info.getName());
            if (info.getName().startsWith("NR_"))
                continue; // Escape some weird entries
            // Want to remove 'chr' from the chromosome name
            info.setChromosome(tokens[2].substring(3));
            info.setStrand(tokens[3]);
            // Start coordinates should add 1 in the UCSC genome table browser
            info.setStart(Integer.parseInt(tokens[4]) + 1);
            info.setEnd(Integer.parseInt(tokens[5]));
            info.setGeneName(tokens[12]);
            Transcript transcript = new Transcript();
            info.setTranscript(transcript);
            transcript.setName(info.getGeneName());
            transcript.setStart(info.getStart());
            transcript.setEnd(info.getEnd());
            transcript.setExons(extractExons(tokens));
//            System.out.println("Gene: " + info.getGeneName());
//            System.out.println("Total genomic length: " + info.getLength());
//            System.out.println("Total cDNA length: " + transcript.getLength());
            try {
                transcript.setCDSCoordinates(Integer.parseInt(tokens[6]) + 1,
                                             Integer.parseInt(tokens[7]));
            }
            catch(IllegalArgumentException e) {
                System.out.println("Parsing " + info.getName() + ": " + e.getMessage());
                error ++;
                continue;
            }
//            System.out.println("Coding region length: " + transcript.getCds().getLength());
            list.add(info);
//            if (count ++ == 10)
//                break;
        }
        fu.close();
        System.out.println("Total parsing errors: " + error);
        return list;
    }
    
    private List<Exon> extractExons(String[] tokens) {
        List<Exon> exons = new ArrayList<Exon>();
        String[] exonsStart = tokens[9].split(",");
        String[] exonsEnd = tokens[10].split(",");
        for (int i = 0; i < exonsStart.length; i++) {
            Exon exon = new Exon();
            exon.setStart(Integer.parseInt(exonsStart[i]) + 1);
            exon.setEnd(Integer.parseInt(exonsEnd[i]));
            exons.add(exon);
        }
        return exons;
    }
    
    public List<RefSeqInfo> loadRefSeqInfo(String refSeqInfoFile) throws IOException {
        List<RefSeqInfo> list = new ArrayList<RefSeqInfo>();
        FileUtility fu = new FileUtility();
        fu.setInput(refSeqInfoFile);
        String line = fu.readLine(); // headers
        // for some checking
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            RefSeqInfo info = new RefSeqInfo();
            info.setName(tokens[1]);
            // Want to remove 'chr' from the chromosome name
            info.setChromosome(tokens[2].substring(3));
            info.setStrand(tokens[3]);
            info.setStart(Integer.parseInt(tokens[4]));
            info.setEnd(Integer.parseInt(tokens[5]));
            info.setGeneName(tokens[12]);
            list.add(info);
        }
        fu.close();
        return list;
    }
    
    public Map<String, Integer> loadChromosomeSizes() throws IOException {
        Map<String, Integer> maps = new HashMap<String, Integer>();
        String fileName = DATA_DIR + "HumanChromosomeStatistics.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        boolean isInTable = false;
        while ((line = fu.readLine()) != null) {
            // There is an extra space
            if (line.startsWith(" ------")) {
                isInTable = !isInTable;
            }
            else if (isInTable) {
                String[] tokens = line.trim().split("( )+");
                maps.put(tokens[0], new Integer(tokens[1]));
            }
        }
        fu.close();
        return maps;
    }
    
    @Test
    public void testMethod() throws IOException {
        Set<String> geneNames = loadHumanGeneNames();
        System.out.println("Gene Names: " + geneNames.size());
        List<RefSeqInfo> refSeqs = loadHumanRefSeqInfos();
        Set<String> geneNames1 = new HashSet<String>();
        for (RefSeqInfo refSeq : refSeqs)
            geneNames1.add(refSeq.getGeneName());
        System.out.println("Gene Names from RefSeq: " + geneNames1.size());
    }
    
    public Map<String, String> mapToUniProt(Collection<String> geneNames) throws IOException {
        if (geneToUniMap == null) {
            UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
            geneToUniMap = uniAnalyzer.loadGeneNameToUniProt();
        }
        Map<String, String> nameToUniMap = new HashMap<String, String>();
        for (String name : geneNames) {
            String uni = geneToUniMap.get(name);
            if (uni != null)
                nameToUniMap.put(name, uni);
        }
        //System.out.println("Total Mapped: " + nameToUniMap.size());
        for (String name : geneNames) {
            if (nameToUniMap.containsKey(name))
                continue;
            String hprdId = geneToHprdMap.get(name);
            if (hprdId != null)
                nameToUniMap.put(name, hprdId);
            //else
            //    System.out.println("Cannot mapping to hprd: " + name);
        }
        //System.out.println("After HPRD mappping: " + nameToUniMap.size());
        return nameToUniMap;
    }
    
    public Set<String> loadHumanGeneNames() throws IOException {
        String fileName = DATA_DIR + "HumanRefSeqGenes.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> geneNames = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            geneNames.add(tokens[12]);
        }
        fu.close();
        return geneNames;
    }
    
}
