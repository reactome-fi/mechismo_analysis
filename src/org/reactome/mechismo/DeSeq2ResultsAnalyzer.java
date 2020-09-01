package org.reactome.mechismo;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.apache.commons.math3.stat.inference.OneWayAnova;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;

/**
 * This class is used to analyze DESEQ2 different gene expression results produced by Francesco.
 * @author wug
 *
 */
public class DeSeq2ResultsAnalyzer {
    private FileUtility fu = new FileUtility();
    private String dirName = "datasets/Mechismo/FrancescoResults/DESEQ2";
    
    public DeSeq2ResultsAnalyzer() {
    }
    
    /**
     * This method is used to check GSEA significant pathways for THCA and SKCM NRAS
     * results. The GSEA analysis was conducted using the local set up with 100 permutations,
     * minSize = 5, and maxSize = 1000, GMT file = ReactomePathways_122718.gmt
     * @throws IOException
     */
    @Test
    public void checkOverlapGSEAResults() throws Exception {
        double fdrCutoff = 0.01d;
        String dirName = "results/DESEQ2";
        String fileName = "GSEA_THCA_NRAS_RASAL2_Ranked.txt";
        Set<String> thcaPathways = loadGSEAPathways(dirName, fileName, fdrCutoff);
        System.out.println("FDR cutoff: " + fdrCutoff);
        System.out.println("Total THCA pathways: " + thcaPathways.size());
        fileName = "GSEA_SKCM_NRAS_RASAL2_Ranked.txt";
        Set<String> skcmPathways = loadGSEAPathways(dirName, fileName, fdrCutoff);
        System.out.println("Total SKCM pathways: " + skcmPathways.size());
        Set<String> shared = InteractionUtilities.getShared(thcaPathways, skcmPathways);
        System.out.println("Shared: " + shared.size());
        int total = 1534 * 2; // or 1535. Different cancers have a little bit of different total
                              // Time 2 for two directions
        double pvalue = MathUtilities.calculateHypergeometricPValue(total,
                                                                    thcaPathways.size(),
                                                                    skcmPathways.size(),
                                                                    shared.size());
        System.out.println("p-value: " + pvalue);
        if (fdrCutoff <= 0.01d)
            shared.stream().sorted(Comparator.naturalOrder()).forEach(System.out::println);
    }
    
    private Set<String> loadGSEAPathways(String dirName, String fileName, double fdrCutoff) throws IOException {
        return Files.lines(Paths.get(dirName, fileName))
                .skip(1)
                .map(line -> line.split("\t"))
                .filter(tokens -> new Double(tokens[4]) <= fdrCutoff)
                .map(tokens -> tokens[0] + "\t" + tokens[5]) // For two directions
                .collect(Collectors.toSet());
    }
    
    @Test
    public void generateRankedGeneList() throws IOException {
//        String fileName = "THCA_NRAS_RASAL2__results.csv";
//        String outFileName = "results/DESEQ2/THCA_NRAS_RASAL2_Ranked.txt";
        
        String fileName = "SKCM_NRAS_RASAL2__results.csv";
        String outFileName = "results/DESEQ2/SKCM_NRAS_RASAL2_Ranked.txt";
        
        fu.setInput(dirName + "/" + fileName);
        fu.setOutput(outFileName);
        String line = fu.readLine();
        fu.printLine("Gene\tStat");
        Map<String, Set<Double>> geneToStats = new HashMap<>();
        while ((line = fu.readLine()) != null) {
            if (line.contains("c(\""))
                continue;   // Escape an id can be mapped to multiple gene, which brings us some parsing difficult
            String[] tokens = line.split(",");
            String gene = tokens[tokens.length - 1];
            if (gene.equals("NA"))
                continue; // Escape NA
            geneToStats.compute(gene, (key, set) -> {
                if (set == null)
                    set = new HashSet<>();
                set.add(new Double(tokens[4])); // the stat column
                return set;
            });
        }
        List<String> geneList = new ArrayList<String>(geneToStats.keySet());
        Collections.sort(geneList);
        for (String gene : geneList) {
            Set<Double> stats = geneToStats.get(gene);
            double mean = getAverage(stats);
            fu.printLine(gene + "\t" + mean);
        }
        fu.close();
    }
    
    private double getAverage(Set<Double> stats) {
        int counter = 0;
        double total = 0.0d;
        for (Double stat : stats) {
            counter ++;
            total += stat;
        }
        return total / counter;
    }
    
    /**
     * Check overlapping of FIs between cancers for one specific gene
     * @throws IOException
     */
    @Test
    public void checkFIOverlapping() throws IOException {
        String gene = "NRAS";
        String fileName = "results/DESEQ2/CancerFIsDist.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, Set<String>> cancerToFIs = new HashMap<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens[0].contains(gene)) {
                collectFIForCancer(cancerToFIs, tokens[0]);
            }
            if (tokens[1].contains(gene)) {
                collectFIForCancer(cancerToFIs, tokens[1]);
            }
        }
        fu.close();
        List<String> cancerList = new ArrayList<>(cancerToFIs.keySet());
        Collections.sort(cancerList);
        for (int i = 0; i < cancerList.size() - 1; i++) {
            String cancer1 = cancerList.get(i);
            Set<String> fis1 = cancerToFIs.get(cancer1);
            System.out.println(cancer1 + "\t" + fis1.size());
            for (int j = i + 1; j < cancerList.size(); j++) {
                String cancer2 = cancerList.get(j);
                Set<String> fis2 = cancerToFIs.get(cancer2);
                System.out.println(cancer2 + "\t" + fis2.size());
                Set<String> shared = InteractionUtilities.getShared(fis1, fis2);
                System.out.println("Shared: " + shared.size());
                shared.forEach(System.out::println);
            }
        }
    }

    private void collectFIForCancer(Map<String, Set<String>> cancerToFIs, 
                                    String token) {
        String[] tokens1 = token.split("_");
        cancerToFIs.compute(tokens1[0], (key, set) -> {
            if (set == null)
                set = new HashSet<>();
            // Need to normalize FIs
            String fi = InteractionUtilities.generateFIFromGene(tokens1[1], tokens1[2]);
            set.add(fi);
            return set;
        });
    }
    
    /**
     * Generate the distances generated for cancer FIs based on DESeq2 differentilal
     * analysis results.
     * @throws IOException
     */
    @Test
    public void checkCancerFIsDist() throws IOException {
        String fileName = "results/DESEQ2/CancerFIsDist.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> cancers = new HashSet<>();
        Set<String> genes = new HashSet<>();
        Set<String> fis = new HashSet<>();
        SummaryStatistics stat = new SummaryStatistics();
        List<String> lines = new ArrayList<>();
        while ((line = fu.readLine()) != null) {
            lines.add(line);
            String[] tokens = line.split("\t");
            for (int i = 0; i < tokens.length - 1; i++) {
                String[] tokens1 = tokens[i].split("_");
                cancers.add(tokens1[0]);
                genes.add(tokens1[1]);
                genes.add(tokens1[2]);
                fis.add(InteractionUtilities.generateFIFromGene(tokens1[1],
                                                                tokens1[2]));
            }
            stat.addValue(new Double(tokens[2]));
        }
        fu.close();
        System.out.println("cancers: " + cancers.size());
        System.out.println("genes: " + genes.size());
        System.out.println("fis: " + fis.size());
        Map<String, Set<String>> geneToFIs = InteractionUtilities.generateProteinToPartners(fis);
        List<String> geneList = new ArrayList<>(geneToFIs.keySet());
        geneList.sort((gene1, gene2) -> geneToFIs.get(gene2).size() - geneToFIs.get(gene1).size());
        for (String gene : geneList) {
            System.out.println(gene + "\t" + geneToFIs.get(gene).size());
        }
        System.out.println("Average distance: " + stat.getMean());
        System.out.println("SD: " + stat.getStandardDeviation());
        
        // Perform one-way anova for cancer for each gene
        System.out.println("\nPerform ANOVA for cancer: ");
        OneWayAnova anova = new OneWayAnova();
        for (String gene : geneList) {
            // Nothing is interesting
            if (geneToFIs.get(gene).size() < 3)
                break;
            Map<String, SummaryStatistics> cancerToStat = new HashMap<>();
            for (String line1 : lines) {
                String[] tokens = line1.split("\t");
                // Check within cancer
                String[] tokens1 = tokens[0].split("_");
                String[] tokens2 = tokens[1].split("_");
                if (!tokens1[0].equals(tokens2[0]))
                    continue;
                if (!tokens[0].contains(gene) || !tokens[1].contains(gene))
                    continue; // Make sure both have this gene
                SummaryStatistics stat1 = cancerToStat.get(tokens1[0]);
                if (stat1 == null) {
                    stat1 = new SummaryStatistics();
                    cancerToStat.put(tokens1[0], stat1);
                }
                stat1.addValue(new Double(tokens[2]));
            }
            if (cancerToStat.size() == 1) {
                System.out.println(gene + "\tone cancer only!");
                continue;
            }
            System.out.println(gene + "\t" + anova.anovaPValue(cancerToStat.values(), true));
            if (gene.equals("NRAS")) {
                for (String cancer : cancerToStat.keySet()) {
                    SummaryStatistics cancerStat = cancerToStat.get(cancer);
                    System.out.println(cancer + ": " + cancerStat);
                }
            }
        }
        
        // Perform a comparison test between within cancer distances and between cancer distances
        // for a specific gene
        System.out.println("\nPerform Mann Whitney U Test for between and within cancer:");
        MannWhitneyUTest utest = new MannWhitneyUTest();
        for (String gene : geneList) {
            if (geneToFIs.get(gene).size() < 3)
                break;
            List<Double> within = new ArrayList<>();
            List<Double> between = new ArrayList<>();
            for (String line1 : lines) {
                String[] tokens = line1.split("\t");
                // Check within cancer
                String[] tokens1 = tokens[0].split("_");
                String[] tokens2 = tokens[1].split("_");
                if (!tokens[0].contains(gene) || !tokens[1].contains(gene))
                    continue; // Make sure both have this gene
                if (tokens1[0].equals(tokens2[0]))
                    within.add(new Double(tokens[2]));
                else
                    between.add(new Double(tokens[2]));
            }
            if (within.size() < 3 || between.size() < 3) {
                System.out.println(gene + "\tNot enough distances!");
                continue;
            }
            double pvalue = utest.mannWhitneyUTest(converToArray(within),
                                                   converToArray(between));
            System.out.println(gene + "\t" + pvalue);
            if (gene.equals("NRAS") ) {
                SummaryStatistics distStat = new SummaryStatistics();
                within.forEach(dist -> distStat.addValue(dist));
                System.out.println("Within distances: average " + distStat.getMean() + 
                                   ", sd " + distStat.getStandardDeviation());
                distStat.clear();
                between.forEach(dist -> distStat.addValue(dist));
                System.out.println("Between distances: average " + distStat.getMean() + 
                                   ", sd " + distStat.getStandardDeviation());
            }
        }
    }
    
    private double[] converToArray(List<Double> list) {
        double[] rtn = new double[list.size()];
        for (int i = 0; i < list.size(); i++)
            rtn[i] = list.get(i);
        return rtn;
    }
    
    /**
     * Generate a matrix from significant impacted FIs and differential gene expression
     * scores for clustering analysis.
     * @throws IOException
     */
    @Test
    public void generateMatrixFIAndDiffGeneExp() throws IOException {
        String localDirName = dirName;
        Set<String> geneIds = new HashSet<>();
        File dir = new File(localDirName);
        String line = null;
        String[] tokens = null;
        for (File file : dir.listFiles()) {
            if (!file.getName().endsWith(".csv"))
                continue;
            fu.setInput(file.getAbsolutePath());
            line = fu.readLine();
            while ((line = fu.readLine()) != null) {
                // Escape rows that don't have genes
                tokens = line.split(",");
                if (tokens[tokens.length - 1].equals("NA"))
                    continue;
                geneIds.add(tokens[0]);
            }
            fu.close();
        }
        System.out.println("Total gene ids: " + geneIds.size());
        
        String outFileName = localDirName + "/FIDiffStatMatrix.txt";
        fu.setOutput(outFileName);
        List<String> geneIdList = new ArrayList<>(geneIds);
        Collections.sort(geneIdList);
        StringBuilder builder = new StringBuilder();
        builder.append("FI");
        geneIdList.stream()
                  .map(id -> id.subSequence(1, id.length() - 1)) // Remove quotation marks
                  .forEach(id -> builder.append("\t").append(id));
        fu.printLine(builder.toString());
        builder.setLength(0);
        Map<String, Double> geneIdToValue = new HashMap<>();
        FileUtility readFU = new FileUtility();
        int fileCounter = 0;
        for (File file : dir.listFiles()) {
            if (!file.getName().endsWith(".csv"))
                continue;
            readFU.setInput(file.getAbsolutePath());
            line = readFU.readLine();
            while ((line = readFU.readLine()) != null) {
                tokens = line.split(",");
                if (!geneIds.contains(tokens[0]))
                    continue;
                geneIdToValue.put(tokens[0],
                                  new Double(tokens[4])); // Use the stat column
            }
            readFU.close();
            // Get the FI name from the file name
            String fiName = file.getName();
            int index = fiName.indexOf("__");
            fiName = fiName.substring(0, index);
            builder.append(fiName);
            for (String geneId : geneIdList) {
                Double value = geneIdToValue.get(geneId);
                builder.append("\t");
                if (value == null)
                    builder.append("NA");
                else
                    builder.append(value);
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
            geneIdToValue.clear();
            fileCounter ++;
        }
        System.out.println("Total processed files: " + fileCounter);
        fu.close();
    }

}
