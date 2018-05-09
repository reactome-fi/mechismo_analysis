/*
 * Created on Jan 11, 2011
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class group a set of methods that are related to gene expression data analysis.
 * @author wgm
 *
 */
public class CancerGeneExpressionCommon {
    private FileUtility fu = new FileUtility();
    
    public CancerGeneExpressionCommon() {
    }
    
    public void generateSampleToGeneExpValue(Map<String, Double> sampleToValue,
                                             String output,
                                             String label) throws IOException {
        fu.setOutput(output);
        StringBuilder builder = new StringBuilder();
        builder.append("Sample").append("\t").append(label);
        fu.printLine(builder.toString());
        for (String sample : sampleToValue.keySet()) {
            builder.setLength(0);
            builder.append(sample);
            
            builder.append("\t");
            if (sampleToValue.get(sample) == null)
                builder.append("NA");
            else
                builder.append(sampleToValue.get(sample));
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    /**
     * Use this method to output a matrix with samples as rows, genes as columns.
     * @param geneToSampleToValue
     * @param output
     * @throws IOExcepiton
     */
    public void generateSampleToGeneToExpValue(Map<String, Map<String, Double>> geneToSampleToValue,
                                               String output) throws IOException {
        // Get a list of samples and sort them
        Set<String> samples = CancerAnalysisUtilitites.grepSamplesInGeneToSampleToValue(geneToSampleToValue);
        List<String> sampleList = new ArrayList<String>(samples);
        Collections.sort(sampleList);
        generateSampleToGeneToExpValue(geneToSampleToValue, sampleList, output);
    }

    /**
     * Generate a matrix with samples as rows, genes as colums for a list of samples passed.
     * @param geneToSampleToValue
     * @param sampleList
     * @param output
     * @throws IOException
     */
    public void generateSampleToGeneToExpValue(Map<String, Map<String, Double>> geneToSampleToValue,
                                               Collection<String> sampleList,
                                               String output) throws IOException {
        // Get a list of all genes and sort them for easy debugging.
        List<String> geneList = new ArrayList<String>(geneToSampleToValue.keySet());
        Collections.sort(geneList);
        fu.setOutput(output);
        // Generate a header for genes
        StringBuilder builder = new StringBuilder();
        builder.append("Sample");
        for (String gene : geneList) {
            builder.append("\t").append(gene);
        }
        fu.printLine(builder.toString());
        boolean hasData = false;
        for (String sample : sampleList) {
            builder.setLength(0);
            hasData = false;
            builder.append(sample);
            for (String gene : geneList) {
                Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
                Double value = sampleToValue.get(sample);
                builder.append("\t");
                if (value == null)
                    builder.append("NA");
                else {
                    builder.append(value);
                    hasData = true;
                }
            }
            if (hasData) // Output only lines with data
                fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    private List<String> grepAllSamples(Map<String, Map<String, Double>> geneToSampleToValue) {
        Set<String> samples = CancerAnalysisUtilitites.grepSamplesInGeneToSampleToValue(geneToSampleToValue);
        return new ArrayList<String>(samples);
    }
    
    /**
     * Generate a sample to gene expression file for a sub-set of genes in the original data set.
     * @param geneToSampleToValue
     * @param output
     * @param genes
     * @throws IOException
     */
    public void generateSampleToGeneExpSubset(Map<String, Map<String, Double>> geneToSampleToValue,
                                              String output,
                                              Collection<String> genes) throws IOException {
        List<String> samples = grepAllSamples(geneToSampleToValue);
        fu.setOutput(output);
        StringBuilder builder = new StringBuilder();
        builder.append("Sample");
        for (String gene : genes)
            builder.append("\t").append(gene);
        fu.printLine(builder.toString());
        for (String sample : samples) {
            builder.setLength(0);
            builder.append(sample);
            for (String gene : genes) {
                Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
                builder.append("\t");
                if (sampleToValue == null || sampleToValue.get(sample) == null)
                    builder.append("NA");
                else
                    builder.append(sampleToValue.get(sample));
            }
            fu.printLine(builder.toString());
        }
        fu.close();
    }   
    
    /**
     * Generate sample to module gene expression scores for a list of clusters.
     * @param geneToSampleToValue
     * @param clusters
     * @return
     */
    public List<Map<String, Double>> generateSampleToGeneExpForClusters(Map<String, Map<String, Double>> geneToSampleToValue,
                                                                        List<Set<String>> clusters) {
        List<Map<String, Double>> sampleToClusterValues = new ArrayList<Map<String,Double>>();
        List<String> samples = grepAllSamples(geneToSampleToValue);
        for (Set<String> cluster : clusters) {
            Map<String, Double> sampleToClusterValue = new HashMap<String, Double>();
            sampleToClusterValues.add(sampleToClusterValue);
            for (String sample : samples) {
                double total = 0.0d;
                int count = 0;
                for (String gene : cluster) {
                    Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
                    if (sampleToValue == null || sampleToValue.get(sample) == null)
                        continue;
                    total += sampleToValue.get(sample);
                    count ++;
                }
                sampleToClusterValue.put(sample, total / count);
            }
        }
        return sampleToClusterValues;
    }

    public void generateSampleToGeneExpClusters(Map<String, Map<String, Double>> geneToSampleToValue,
                                                String output,
                                                List<Set<String>> clusters) throws IOException {
        generateSampleToGeneExpClusters(geneToSampleToValue, output, clusters, null);
    }    
    
    public void generateSampleToGeneExpClusters(Map<String, Map<String, Double>> geneToSampleToValue,
                                                String output,
                                                List<Set<String>> clusters,
                                                List<String> clusterNames) throws IOException {
        List<String> samples = grepAllSamples(geneToSampleToValue);
        Collections.sort(samples);
        generateSampleToGeneExpClusters(geneToSampleToValue, samples, output,
                                        clusters, clusterNames);
    }

    public void generateSampleToGeneExpClusters(Map<String, Map<String, Double>> geneToSampleToValue,
                                                 List<String> samples,
                                                 String output,
                                                 List<Set<String>> clusters,
                                                 List<String> clusterNames) throws IOException {
        if (samples == null) {
            samples = grepAllSamples(geneToSampleToValue);
            Collections.sort(samples);
        }
        fu.setOutput(output);
        StringBuilder builder = new StringBuilder();
        builder.append("Sample");
        if (clusterNames == null) {
            for (int i = 0; i < clusters.size(); i++)
                builder.append("\tModule").append(i);
        }
        else {
            for (int i = 0; i < clusters.size(); i++)
                builder.append("\t").append(clusterNames.get(i));
        }
        fu.printLine(builder.toString());
        double total = 0.0d;
        int count = 0;
        for (String sample : samples) {
            builder.setLength(0);
            builder.append(sample);
            for (Set<String> cluster : clusters) {
                total = 0.0d;
                count = 0;
                for (String gene : cluster) {
                    Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
                    if (sampleToValue != null && sampleToValue.get(sample) != null) {
                        total += sampleToValue.get(sample);
                        count ++;
                    }
                }
//                if (count == 0)
//                    builder.append("\t").append(Double.MIN_VALUE);
//                else
                    builder.append("\t").append(total / count);
            }
            fu.printLine(builder.toString());
        }
        fu.close();
    }   
    
    /**
     * Load a cluster list from file exported from FI Cytoscape plug-in.
     * @param clusterAnnotFileName
     * @return
     * @throws IOException
     */
    public List<Set<String>> loadGeneExpClusters(String clusterAnnotFileName) throws IOException {
        fu.setInput(clusterAnnotFileName);
        String line = fu.readLine();
        List<Set<String>> clusters = new ArrayList<Set<String>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            tokens = tokens[4].split(",");
            Set<String> cluster = new HashSet<String>();
            for (String token : tokens)
                cluster.add(token);
            clusters.add(cluster);
        }
        return clusters;
    }
    
    private Map<String, Double> loadFIToValue(String fileName) throws IOException {
        fu.setInput(fileName);
        String line = null;
        Map<String, Double> fiToValue = new HashMap<String, Double>();
        int index = 0;
        while ((line = fu.readLine()) != null) {
            index = line.lastIndexOf("\t");
            String fi = line.substring(0, index);
            String value = line.substring(index + 1);
            fiToValue.put(fi, new Double(value));
        }
        fu.close();
        return fiToValue;
    }
    
    /**
     * Generate edge weights based on gene expression correlations.
     * @param geneExpFileName
     * @param output
     * @throws Exception
     */
    public void generateFIsWithGeneExpCorr(String geneExpFileName,
                                           String output,
                                           String sampleEscapePattern) throws IOException {
        Map<String, Map<String, Double>> geneToSampleToValue = loadGeneExp(geneExpFileName);
        generateFIsWithGeneExpCorr(geneToSampleToValue, 
                                   output, 
                                   sampleEscapePattern);
    }
    
    /**
     * This is the actual method to calculate Pearson correlation for FIs based on provided
     * gene gene data set in a Map.
     * @param geneToSampleToValue
     * @param fis
     * @param useAbsoluteValue
     * @param sampleEscapePattern
     * @return
     */
    public Set<String> calculateGeneExpCorrForFIs(Map<String, Map<String, Double>> geneToSampleToValue,
                                                  Set<String> fis,
                                                  Boolean useAbsoluteValue,
                                                  String sampleEscapePattern) {
        Set<String> fisWithCorrs = new HashSet<String>();
        int index = 0;
        String gene1, gene2;
        List<Double> values1 = new ArrayList<Double>();
        List<Double> values2 = new ArrayList<Double>();
//        long time1 = System.currentTimeMillis();
        for (String fi : fis) {
//            fisWithCorrs.add(fi + "\t1.0");
//            if (true)
//                continue;
            index = fi.indexOf("\t");
            gene1 = fi.substring(0, index);
            Map<String, Double> sampleToValue1 = geneToSampleToValue.get(gene1);
            if (sampleToValue1 == null)
                continue;
            gene2 = fi.substring(index + 1);
            Map<String, Double> sampleToValue2 = geneToSampleToValue.get(gene2);
            if (sampleToValue2 == null)
                continue;
            for (String sample : sampleToValue1.keySet()) {
                if (sampleEscapePattern != null && sample.contains(sampleEscapePattern))
                    continue; // Escape normal samples
                Double value1 = sampleToValue1.get(sample);
                Double value2 = sampleToValue2.get(sample);
                // Always escape if there is any null value
                if (value1 == null || value2 == null)
                    continue;
                values1.add(value1);
                values2.add(value2);
            }
            double corr = MathUtilities.calculatePearsonCorrelation(values1, values2);
            // The following check is very slow. Need to depends a pre-processing to 
            // remove any no correct values in the array data set.
//            // Just in case
//            if ((corr + "").equals(nan)) // This is much faster than Double.isNaN(double)
//                continue;
            if (useAbsoluteValue)
                fisWithCorrs.add(fi + "\t" + Math.abs(corr));
            else
                fisWithCorrs.add(fi + "\t" + corr);
            values1.clear();
            values2.clear();
        }
//        long time2 = System.currentTimeMillis();
//        System.out.println("Time for correaltion caluation: " + (time2 - time1));
        return fisWithCorrs;
    }

    /**
     * Generate edge weights based on gene expression correlation.
     * @param geneToSampleToValue load gene expression data set in a map.
     * @param output
     * @param fiFileName
     * @param sampleEscapePattern
     * @throws IOException
     */
    public void generateFIsWithGeneExpCorr(Map<String, Map<String, Double>> geneToSampleToValue,
                                           String output, 
                                           String sampleEscapePattern) throws IOException {
        // Always use the big graph component.
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        long time1 = System.currentTimeMillis();
        Set<String> fisWithCorrs = calculateGeneExpCorrForFIs(geneToSampleToValue, 
                                                              fis, 
                                                              true, 
                                                              sampleEscapePattern);
        fu.setOutput(output);
        for (String fiWithCorr : fisWithCorrs)
            fu.printLine(fiWithCorr);
        long time2 = System.currentTimeMillis();
        System.out.println("Time for correlations(): " + (time2 - time1));
        fu.close();
    }
    
    /**
     * Use this method to filter out genes having null values for some samples.
     * @param geneToSampleToValue
     */
    public void filterOutGenesWithNullValues(Map<String, Map<String, Double>> geneToSampleToValue) {
        List<String> allSamples = grepAllSamples(geneToSampleToValue);
        for (Iterator<String> it = geneToSampleToValue.keySet().iterator(); it.hasNext();) {
            Map<String, Double> sampleToValue = geneToSampleToValue.get(it.next());
            if (sampleToValue.size() < allSamples.size())
                it.remove();
        }
    }
    
    public Map<String, Map<String, Double>> loadGeneExp(String fileName) throws IOException {
        return loadGeneExp(fileName, true);
    }
    
    /**
     * This method is used to load a gene expression file in the format of gene to sample to value.
     * Use this method for a data file with multiple rows mapped to the same gene. The values for a
     * gene will be averaged if multiple rows exist for that gene.
     * @param fileName
     * @return
     * @throws IOException
     */
    public Map<String, Map<String, Double>> loadGeneExpAfterAverage(String fileName) throws IOException {
        Map<String, Map<String, Double>> geneToSampleToValue = new HashMap<String, Map<String,Double>>();
        // Used for counting for averaging
        Map<String, Map<String, Integer>> geneToSampleToCount = new HashMap<String, Map<String,Integer>>();
        fu.setInput(fileName);
        String line = fu.readLine();
        // Get sample list
        String[] samples = line.split("\t");
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[0];
            Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
            if (sampleToValue == null) {
                sampleToValue = new HashMap<String, Double>();
                geneToSampleToValue.put(gene, sampleToValue);
                geneToSampleToCount.put(gene, new HashMap<String, Integer>());
            }
            Map<String, Integer> sampleToCount = geneToSampleToCount.get(gene);
            for (int i = 1; i < tokens.length; i++) {
                if (tokens[i].equals("NA")) // Just escape NA
                    continue;
                String sample = samples[i];
                Double value = new Double(tokens[i]);
                Double current = sampleToValue.get(sample);
                Integer count = sampleToCount.get(sample);
                if (current == null) {
                    sampleToValue.put(sample, value);
                    sampleToCount.put(sample, 1);
                }
                else {
                    sampleToValue.put(sample, value + current);
                    sampleToCount.put(sample, ++count);
                }
            }
        }
        fu.close();
        // Do a average
        for (String gene : geneToSampleToValue.keySet()) {
            Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
            Map<String, Integer> sampleToCount = geneToSampleToCount.get(gene);
            for (String sample : sampleToCount.keySet()) {
                Integer count = sampleToCount.get(sample);
                if (count == 1)
                    continue;
                Double value = sampleToValue.get(sample);
                sampleToValue.put(sample, value / count);
            }
        }
        return geneToSampleToValue;
    }
    
    /**
     * Load gene expression data per gene. The value is stored as gene to sample to value.
     * @return
     * @throws IOException
     */
    public Map<String, Map<String, Double>> loadGeneExp(String fileName,
                                                        boolean escapeRowContainNA) throws IOException {
        Map<String, Map<String, Double>> geneToData = new HashMap<String, Map<String,Double>>();
        fu.setInput(fileName);
        // Sample list
        String line = fu.readLine();
        List<String> sampleList = new ArrayList<String>();
        String[] tokens = line.split("\t");
        for (String token : tokens) {
            String sample = token.replaceAll("\"", "");
            //if (sample.length() > 12)
            //    sample = sample.substring(0, 12);
            sampleList.add(sample);
        }
        // Starting parsing
        while ((line = fu.readLine()) != null) {
//            if (line.contains("NA"))
//                continue; // Don't want any genes containing "NA".
//
            int index = line.indexOf("\t");
            if (line.substring(index + 1).contains("NA") && escapeRowContainNA)
                continue; // Don't want any genes with values containing "NA".
//            System.out.println(line);
            tokens = line.split("\t");
            // The first one is gene name
            String gene = tokens[0].replace("\"", "");
            Map<String, Double> sampleToValue = new HashMap<String, Double>();
            for (int i = 1; i < tokens.length; i++) {
                if (tokens[i].equals("NA") || tokens[i].length() == 0)
                    continue; // Just escape these values
                String sample = sampleList.get(i); // The first sample has been checked.
                Double value = new Double(tokens[i]);
                sampleToValue.put(sample, value);
            }
            if (geneToData.containsKey(gene)) {
                System.out.println("Duplicated gene: " + gene);
            }
            geneToData.put(gene, sampleToValue);
        }
        fu.close();
        // These two genes having NA values in exp
        //geneToData.remove("C1ORF129");
        //geneToData.remove("LRRC50");
        return geneToData;
    }
    
    /**
     * This method is used to output gene expression data set in a Map.
     * @param geneToSampleToValue
     * @param outFileName
     * @throws IOException
     */
    public void outputGeneExp(Map<String, Map<String, Double>> geneToSampleToValue, 
                              String outFileName) throws IOException {
        // Sort samples
        List<String> sampleList = grepAllSamples(geneToSampleToValue);
        Collections.sort(sampleList);
        // Sort genes
        List<String> geneList = new ArrayList<String>(geneToSampleToValue.keySet());
        Collections.sort(geneList);
        fu.setOutput(outFileName);
        StringBuilder builder = new StringBuilder();
        builder.append("GENE");
        for (String sample : sampleList)
            builder.append("\t").append(sample);
        fu.printLine(builder.toString());
        builder.setLength(0);
        for (String gene : geneList) {
            Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
            builder.append(gene);
            for (String sample : sampleList) {
                Double value = sampleToValue.get(sample);
                builder.append("\t");
                if (value == null)
                    builder.append("NA");
                else
                    builder.append(value);
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    /**
     * Do a gene based transformation.
     * @param srcFileName
     * @param targetFileName
     * @throws IOException
     */
    public void zscoreTransformDataOnGene(String srcFileName,
                                          String targetFileName) throws IOException {
        fu.setInput(srcFileName);
        String line = fu.readLine();
        fu.setOutput(targetFileName);
        fu.printLine(line);
        StringBuilder builder = new StringBuilder();
        DescriptiveStatistics stat = new DescriptiveStatistics();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            builder.append(tokens[0]);
            for (int i = 1; i < tokens.length; i++) {
                String token = tokens[i];
                if (token.length() == 0 || token.equals("NA"))
                    continue;
                stat.addValue(new Double(token));
            }
            double sd = stat.getStandardDeviation();
            double mean = stat.getMean();
            for (int i = 1; i < tokens.length; i++) {
                String token = tokens[i];
                if (token.length() == 0 || token.equals("NA"))
                    builder.append("\t").append(token);
                else {
                    Double value = (new Double(token) - mean) / sd;
                    builder.append("\t").append(value);
                }
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
            stat.clear();
        }
        fu.close();
    }
    
    public void ztransformDataGlobally(String srcFileName,
                                       String targetFileName) throws Exception {
        ztransformDataGlobally(srcFileName, targetFileName, true);
    }
    
    /**
     * This method is used to do a global z-score transformation.
     * @param srcFileName
     * @param targetFileName
     * @throws Exception
     */
    public void ztransformDataGlobally(String srcFileName,
                                       String targetFileName,
                                       boolean escapeLineContainingNA) throws Exception {
        Map<String, Map<String, Double>> geneToSampleToExp = loadGeneExp(srcFileName, escapeLineContainingNA);
        DescriptiveStatistics stat = new DescriptiveStatistics();
        // Need to get all samples first in case some samples don't have values
        Set<String> sampleSet = new HashSet<String>();
        for (String gene : geneToSampleToExp.keySet()) {
            Map<String, Double> sampleToExp = geneToSampleToExp.get(gene);
            for (String sample : sampleToExp.keySet()) {
                Double value = sampleToExp.get(sample);
                stat.addValue(value);
                sampleSet.add(sample);
            }
        }
        double mean = stat.getMean();
        double sd = stat.getStandardDeviation();
        System.out.println("Mean: " + mean);
        System.out.println("SD: " + sd);
        
        StringBuilder builder = new StringBuilder();
        fu.setOutput(targetFileName);
        List<String> sampleList = new ArrayList<String>(sampleSet);
        // Sort the samples
        Collections.sort(sampleList);
        // Generate sample header
        builder.append("Gene");
        for (String sample : sampleList) {
            builder.append("\t").append(sample);
        }
        fu.printLine(builder.toString());
        builder.setLength(0);
        // Sort the gene list
        List<String> geneList = new ArrayList<String>(geneToSampleToExp.keySet());
        Collections.sort(geneList);
        for (String gene : geneList) {
            Map<String, Double> sampleToExp = geneToSampleToExp.get(gene);
            builder.append(gene);
            for (String sample : sampleList) {
                Double value = sampleToExp.get(sample);
                if (value == null)
                    builder.append("\t").append("NA");
                else {
                    value = (value - mean) / sd;
                    builder.append("\t").append(value);
                }
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    /**
     * This method is used to load a specific module from a file.
     * @param moduleFileName the file containing a list of modules
     * @param index module index starting from 0.
     * @return
     * @throws IOException
     */
    public Set<String> loadMCLModule(String moduleFileName, int index) throws IOException {
        List<Set<String>> modules = loadMCLModules(moduleFileName);
        return modules.get(index);
    }
    
    /**
     * Load a list of modules from a file formated as following:
     * ......
     * 19: [COL15A1, COL12A1, COL5A3, COL6A3, PTGIS, ISLR, COL6A1, COL14A1, SRPX]
     * 20: [ZNF317, ZNF117, ZNF41, ZNF471, ZNF354C, ZFP28, ZNF266, ZNF154]
     * ......
     * @param moduleFileName
     * @return
     * @throws IOException
     */
    public List<Set<String>> loadMCLModules(String moduleFileName) throws IOException {
        FileUtility fu1 = new FileUtility();
        fu1.setInput(moduleFileName);
        String line = fu1.readLine();
        String[] tokens = line.split("\t");
        if (tokens.length > 1 && tokens[0].equals("Module")) {
            fu1.close();
            return loadPlugInMCLModules(moduleFileName);
        }
        // Check if it is a plug-in file
        fu.setInput(moduleFileName);
        int count = 0;
        List<Set<String>> modules = new ArrayList<Set<String>>();
        while ((line = fu.readLine()) != null) {
            Set<String> module = new HashSet<String>();
            // Load module
            int index1 = line.indexOf("[");
            int index2 = line.lastIndexOf("]");
            String sub = line.substring(index1 + 1, index2);
            tokens = sub.split(", ");
            for (String token : tokens)
                module.add(token);
            modules.add(module);
        }
        fu.close();
        return modules;
    }
    
    private List<Set<String>> loadPlugInMCLModules(String fileName) throws IOException {
        List<Set<String>> modules = new ArrayList<Set<String>>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String last = tokens[tokens.length - 1];
            String[] tmp = last.split(",");
            Set<String> cluster = new HashSet<String>();
            for (String gene : tmp)
                cluster.add(gene);
            modules.add(cluster);
        }
        fu.close();
        return modules;
    }
}
