/*
 * Created on Jun 9, 2009
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.gk.util.StringUtils;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to do network based clustering analysis.
 * @author wgm
 *
 */
public class NetworkClusterAnalyzer {
    private FileUtility fu = new FileUtility();
    
    public NetworkClusterAnalyzer() {
    }
    
    /**
     * This method is used to check overlapping between two network clustering results.
     * @param modules1
     * @param modules2
     * @param size
     * @param genes
     */
    public void checkNetworkModuleOverlapping(List<Set<String>> modules1,
                                              List<Set<String>> modules2,
                                              int size,
                                              double pvalueCutff,
                                              Set<String> genes) throws Exception {
        System.out.println("Module1\tSize1\tModule2\tSize2\tShared\tP-Value");
        for (int i = 0; i < modules1.size(); i++) {
            Set<String> module1 = modules1.get(i);
            if (module1.size() < size)
                continue;
            for (int j = 0; j < modules2.size(); j++) {
                Set<String> module2 = modules2.get(j);
                if (module2.size() < size)
                    continue;
                Set<String> shared = InteractionUtilities.getShared(module1, module2);
                double pvalue = MathUtilities.calculateHypergeometricPValue(genes.size(),
                                                                            module1.size(),
                                                                            module2.size(),
                                                                            shared.size());
                if (pvalue <= pvalueCutff)
                    System.out.println(i + "\t" + module1.size() + "\t" +
                                       j + "\t" + module2.size() + "\t" +
                                       shared.size() + "\t" + pvalue);
            }
        }
    }
    
    public void permutationTestForSamplesInNetworkClusters(Map<String, Set<String>> sampleToAlteredGenes,
                                                           String clusterFileName,
                                                           List<Integer> targetClusters,
                                                           double target,
                                                           boolean permutateClusters) throws Exception {
        List<Set<String>> clusters = loadNetworkClusters(clusterFileName);
        permutationTestForSamplesInNetworkClusters(sampleToAlteredGenes,
                                                   clusters,
                                                   targetClusters, 
                                                   target,
                                                   permutateClusters);
    }

    public void permutationTestForSamplesInNetworkClusters(Map<String, Set<String>> sampleToAlteredGenes,
                                                           List<Set<String>> clusters,
                                                           List<Integer> targetClusters,
                                                           double target,
                                                           boolean permutateClusters)
            throws IOException {
        // Will focus on FI genes only since only we know these genes
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        //Set<String> fis = fu.loadInteractions(R3Constants.RESULT_DIR + "FIsInGene_Pathway_041709.txt");
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> totalGenes = new HashSet<String>();
        for (String sample : sampleToAlteredGenes.keySet()) {
            Set<String> genes = sampleToAlteredGenes.get(sample);
            genes.retainAll(fiGenes);
            //System.out.println(sample + ": " + genes.size());
            System.out.print(genes.size() + ", ");
        }
        System.out.println();
        for (Set<String> set : sampleToAlteredGenes.values())
            totalGenes.addAll(set);
        System.out.println("Total genes: " + totalGenes.size());
        
//        int index = 0;
//        for (Set<String> cluster : clusters) {
//            System.out.println(index + ": " + cluster.size());
//            index ++;
//        }
        int permutationNumber = 10000;
        List<Double> permutPercentages = new ArrayList<Double>();
        Map<String, Set<String>> sampleToPermuGenes = new HashMap<String, Set<String>>();
        RandomData randomizer = new RandomDataImpl();
        if (permutateClusters) {
            for (int j = 0; j < permutationNumber; j++) {
                // Do a quick check using random clusters
                List<Set<String>> randomClusters = generateRandomClusters(clusters, 
                                                                          totalGenes, 
                                                                          randomizer);
                int counter1 = 0;
                for (String sample : sampleToAlteredGenes.keySet()) {
                    Set<String> genes = sampleToAlteredGenes.get(sample);
                    // Want to get the list of clusters
                    List<Integer> clusterIds = new ArrayList<Integer>();
                    for (int i = 0; i < randomClusters.size(); i++) {
                        Set<String> set = randomClusters.get(i);
                        for (String gene : genes) {
                            if (set.contains(gene)) {
                                clusterIds.add(i);
                                break;
                            }
                        }
                    }
                    if (isSampleInClusters(clusterIds, 
                                           targetClusters))
                        counter1 ++;
                }
                double percent = (double) counter1 / sampleToAlteredGenes.size();
                permutPercentages.add(percent);
            }
        }
        else {
//            // The following is testing code
//            sampleToAlteredGenes = new CancerResequenceDataSetAnalyzer().getScienceGBMSampleToAlteredGenes();
//            for (String sample : sampleToAlteredGenes.keySet()) {
//                Set<String> set = sampleToAlteredGenes.get(sample);
//                set.retainAll(fiGenes);
//            }
            
            for (int j = 0; j < permutationNumber; j++) {
                // Need to generate a same number of genes
                for (String sample : sampleToAlteredGenes.keySet()) {
                    Set<String> genes = sampleToAlteredGenes.get(sample);
                    Set<String> permutated = null;
                    if (genes.size() == 0)
                        permutated = new HashSet<String>();
                    else
                        permutated = MathUtilities.randomSampling(totalGenes, 
                                                                  genes.size(),
                                                                  randomizer);
                    sampleToPermuGenes.put(sample, permutated);
                }
                // Check how many samples having genes in target clusters
                int counter = 0;
                for (String sample : sampleToPermuGenes.keySet()) {
                    Set<String> genes = sampleToPermuGenes.get(sample);
                    // Want to get the list of clusters
                    List<Integer> clusterIds = new ArrayList<Integer>();
                    for (int i = 0; i < clusters.size(); i++) {
                        Set<String> set = clusters.get(i);
                        for (String gene : genes) {
                            if (set.contains(gene)) {
                                clusterIds.add(i);
                                break;
                            }
                        }
                    }
                    if (isSampleInClusters(clusterIds, 
                                           targetClusters))
                        counter ++;
                }
                double percent = (double) counter / sampleToPermuGenes.size();
                permutPercentages.add(percent);
                sampleToPermuGenes.clear();
            }
        }
        Collections.sort(permutPercentages, new Comparator<Double>() {
            public int compare(Double v1, Double v2) {
                return v2.compareTo(v1);
            }
        });
//        for (int i = 0; i < permutPercentages.size(); i++)
//            System.out.println(i + ": " + permutPercentages.get(i));
        calculatePValue(target, 
                        permutationNumber, 
                        permutPercentages);
//        String outFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609Permutation10000_071609.txt";
//        fu.setOutput(outFileName);
//        for (Double d : permutPercentages)
//            fu.printLine(d + "");
//        fu.close();
    }
    
    private List<Set<String>> generateRandomClusters(List<Set<String>> clusters,
                                                     Set<String> totalGenes,
                                                     RandomData randomizer) {
        List<Set<String>> samples = new ArrayList<Set<String>>();
        Set<String> copy = new HashSet<String>(totalGenes);
        for (Set<String> cluster : clusters) {
            int size = cluster.size();
            Set<String> tmp = MathUtilities.randomSampling(copy, 
                                                           size,
                                                           randomizer);
            samples.add(tmp);
            copy.removeAll(tmp);
        }
        return samples;
    }
    
    private void calculatePValue(double target, 
                                 int permutationNumber,
                                 List<Double> permutPercentages) {
        int index = -1;
        for (int i = 0; i < permutPercentages.size(); i++) {
            double pert = permutPercentages.get(i);
            if (pert < target) {
                index = i;
                break;
            }
        }
        if (index == 0)
            System.out.println("pvalue < " + 1.0 / permutationNumber);
        else if (index == -1)
            System.out.println("pvalue = 1.0");
        else
            System.out.println("pvalue: " + (double) index / permutationNumber);
        double max = permutPercentages.get(0);
        System.out.println("Max value: " + max);
        System.out.println("Min value: " + permutPercentages.get(permutPercentages.size() - 1));
        System.out.println("Target value: " + target);
    }
    
    public double checkSamplesInClusters(Map<String, Set<String>> sampleToAlteredGenes,
                                         String clusterFileName,
                                         List<Integer> targetClusters) throws IOException {
        List<Set<String>> clusters = loadNetworkClusters(clusterFileName);
        System.out.println("Total clusters: " + clusters.size());
        // Check how many samples having 0 and 1 clusters
        int counter = 0;
        System.out.println("Sample\tClusters\tSharedGenes");
        List<String> posSamples = new ArrayList<String>();
        List<String> negSamples = new ArrayList<String>();
        for (String sample : sampleToAlteredGenes.keySet()) {
            Set<String> genes = sampleToAlteredGenes.get(sample);
            // Want to get the list of clusters
            List<Integer> clusterIds = new ArrayList<Integer>();
            List<String> sharedGenes = new ArrayList<String>();
            for (int i = 0; i < clusters.size(); i++) {
                Set<String> set = clusters.get(i);
                for (String gene : genes) {
                    if (set.contains(gene)) {
                        clusterIds.add(i);
                        break;
                    }
                }
                Set<String> shared = new HashSet<String>(set);
                shared.retainAll(genes);
                sharedGenes.add(StringUtils.join(":", new ArrayList<String>(shared)));
            }
            System.out.println(sample + "\t" + clusterIds + "\t" + sharedGenes);
            if (isSampleInClusters(clusterIds,
                                   targetClusters)) {
                counter ++;
                posSamples.add(sample);
            }
            else
                negSamples.add(sample);
        }
        double percent = (double) counter / sampleToAlteredGenes.size();
//        System.out.println("Total samples checked: " + sampleToAlteredGenes.size());
//        System.out.println("Total samples having 0 and 1: " + counter + "(" + percent + ")");
//        System.out.println("Total samples in cluster 0: " + cluster0);
//        System.out.println("Total samples in cluster 1: " + cluster1);
        return percent;
    }
    
    private boolean isSampleInClusters(List<Integer> clusterIds,
                                       List<Integer> targetClusters) {
        for (Integer target : targetClusters) {
            if (!clusterIds.contains(target))
                return false;
        }
        return true;
//        // Check for a sample has mutation in both cluster 0 and 1.
//        if (clusterIds.contains(0) && 
//            //clusterIds.contains(1) && 
//            clusterIds.contains(2))
//            return true;
////         Check for a sample has a mutation in either cluster 0 or 1.
////        if (clusterIds.contains(0) || clusterIds.contains(1))
////            return true;
////        if (clusterIds.contains(0))
////            return true;
//        return false;
    }
    
    public Set<String> grepAllGenesInClusters(List<Set<String>> clusters) {
        Set<String> set = new HashSet<String>();
        for (Set<String> cluster : clusters)
            set.addAll(cluster);
        return set;
    }
    
    /**
     * Filter a list of network clusters based on a passed cluster size. After filtering,
     * clusters with sizes >= size are kept.
     * @param clusters
     * @param size
     */
    public void filterNetworkClusters(List<Set<String>> clusters,
                                      int size) {
        for (Iterator<Set<String>> it = clusters.iterator(); it.hasNext();) {
            Set<String> cluster = it.next();
            if (cluster.size() < size)
                it.remove();
        }
    }
    
    public List<Set<String>> loadNetworkClusters(String fileName) throws IOException {
        List<Set<String>> rtn = new ArrayList<Set<String>>();
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] headTokens = line.split("\t");
        if (line.equals("Module (class=java.lang.Integer)")) {
            while ((line = fu.readLine()) != null) {
                int index = line.indexOf("=");
                String gene = line.substring(0, index);
                int cluster = new Integer(line.substring(index + 1));
                Set<String> set = null;
                if (cluster < rtn.size()) {
                    set = rtn.get(cluster);
                }
                else {
                    set = new HashSet<String>();
                    rtn.add(set);
                }
                set.add(gene);
            }
        }
        else if ((headTokens.length == 6 || headTokens.length == 5) &&
                 headTokens[0].equals("Module") &&
                 headTokens[headTokens.length - 1].equals("Node List")) { // Format copied from Reactome FI plug-in
            while ((line = fu.readLine()) != null) {
                String[] tokens = line.split("\t");
                tokens = tokens[headTokens.length - 1].split(",");
                Set<String> cluster = new HashSet<String>();
                for (String token : tokens)
                    cluster.add(token);
                rtn.add(cluster);
            }
        }
        else {
            // Based on MCL or other types
            String[] tokens = line.split("\t");
            Set<String> cluster = new HashSet<String>();
            for (String token : tokens)
                cluster.add(token);
            rtn.add(cluster);
            while ((line = fu.readLine()) != null) {
                tokens = line.split("\t");
                cluster = new HashSet<String>();
                for (String token : tokens)
                    cluster.add(token);
                rtn.add(cluster);
            }
        }
        fu.close();
        return rtn;
    }
    
    public void outputNetworkClusters(List<Set<String>> clusterList,
                                      String outFileName) throws IOException {
        fu.setOutput(outFileName);
        // Export as attribute file
        fu.printLine("Module (class=java.lang.Integer)");
        for (int i = 0; i < clusterList.size(); i++) {
            Set<String> cluster = clusterList.get(i);
            for (String gene : cluster)
                fu.printLine(gene + "=" + i);
        }
        fu.close();
    }
    
    /**
     * This method is used to generate a file that maps samples to a lists of booleans.
     * @param sampleToGenes
     * @param selectedGenes
     * @param outFileName
     * @throws IOException
     */
    public void generateSampleToSelectedGenesMatrix(Map<String, Set<String>> sampleToGenes,
                                                   Set<String> selectedGenes,
                                                   String outFileName) throws IOException {
        // Make selected genes as list for easy viewing
        List<String> geneList = new ArrayList<String>(selectedGenes);
        Collections.sort(geneList);
        Map<String, List<Boolean>> sampleToVector = new HashMap<String, List<Boolean>>();
        for (String sample : sampleToGenes.keySet()) {
            Set<String> genes = sampleToGenes.get(sample);
            List<Boolean> values = new ArrayList<Boolean>();
            for (String gene : geneList) {
                values.add(genes.contains(gene));
            }
            sampleToVector.put(sample, values);
        }
        outputSampleToBooleanList(outFileName, 
                                  sampleToVector,
                                  geneList);
    }
    
    public void generateSampleToModulesMatrix(Map<String, Set<String>> sampleToGenes,
                                              String clusterFileName,
                                              String outFileName,
                                              int sizeCutoff) throws Exception {
        Map<String, List<Boolean>> sampleToVector = generateSampleToNetworkModuleVectors(sampleToGenes, 
                                                                                         clusterFileName,
                                                                                         sizeCutoff);
        outputSampleToBooleanList(outFileName, 
                                  sampleToVector,
                                  null);
    }
    
    private void outputSampleToBooleanList(String outFileName,
                                           Map<String, List<Boolean>> sampleToVector,
                                           List<String> headNames)
            throws IOException {
        // Peel the size of list
        List<Boolean> list = sampleToVector.values().iterator().next();
        fu.setOutput(outFileName);
        StringBuilder builder = new StringBuilder();
        // Total cluster
        int totalCluster = list.size();
        builder.append("Sample");
        if (headNames == null) {
            for (int i = 0; i < totalCluster; i++)
                builder.append("\tModule").append(i);
        }
        else {
            for (String name : headNames)
                builder.append("\t").append(name);
        }
        fu.printLine(builder.toString());
        for (String sample : sampleToVector.keySet()) {
            builder.setLength(0);
            builder.append(sample);
            List<Boolean> vector = sampleToVector.get(sample);
            for (Boolean b : vector) {
                builder.append("\t");
                if (b)
                    builder.append(1);
                else
                    builder.append(0);
            }
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    public void generateSampleToModulesMatrixInDegrees(Map<String, Set<String>> sampleToGenes,
                                                       String clusterFileName,
                                                       String outFileName,
                                                       int sizeCutoff) throws Exception {
        Map<String, List<Double>> sampleToVector = generateSampleToNetworkModuleVectorsInDegrees(sampleToGenes,
                                                                                                  clusterFileName,
                                                                                                  sizeCutoff);
        generateSampleToGeneSetDoubleMatrix(outFileName, 
                                            null,
                                            sampleToVector);
    }
    
    private void generateSampleToGeneSetDoubleMatrix(String outFileName,
                                                     List<String> headerNames,
                                                     Map<String, List<Double>> sampleToVector)
            throws IOException {
        // Peel the size of list
        List<Double> list = sampleToVector.values().iterator().next();
        fu.setOutput(outFileName);
        StringBuilder builder = new StringBuilder();
        // Total cluster
        int totalCluster = list.size();
        builder.append("Sample");
        if (headerNames == null)
            for (int i = 0; i < totalCluster; i++)
                builder.append("\tModule").append(i);
        else 
            for (String header : headerNames)
                builder.append("\t").append(header);
        fu.printLine(builder.toString());
        for (String sample : sampleToVector.keySet()) {
            builder.setLength(0);
            builder.append(sample);
            List<Double> vector = sampleToVector.get(sample);
            for (Double b : vector) {
                builder.append("\t").append(b);
            }
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    public void generateSampleToModulesMatrixInPValue(Map<String, Set<String>> sampleToGenes,
                                                      String clusterFileName,
                                                      String outFileName,
                                                      int sizeCutoff,
                                                      int totalGenes) throws Exception {
        Map<String, List<Double>> sampleToVector = generateSampleToNetworkModuleVectorsInPValue(sampleToGenes, 
                                                                                                 clusterFileName,
                                                                                                 sizeCutoff,
                                                                                                 totalGenes);
        generateSampleToGeneSetDoubleMatrix(outFileName,
                                            null,
                                            sampleToVector);
    }
    
    private Map<String, List<Boolean>> generateSampleToNetworkModuleVectors(Map<String, Set<String>> sampleToGenes,
                                                                            String clusterFileName,
                                                                            int sizeCutoff) throws IOException {
        List<Set<String>> clusters = loadNetworkClusters(clusterFileName);
        // Want to create a sample to a boolean vector map
        Map<String, List<Boolean>> sampleToVector = new HashMap<String, List<Boolean>>();
        for (String sample : sampleToGenes.keySet()) {
            Set<String> alteredGenes = sampleToGenes.get(sample);
            List<Boolean> vector = new ArrayList<Boolean>();
            for (Set<String> cluster : clusters) {
                if (cluster.size() < sizeCutoff)
                    continue;
                Set<String> shared = InteractionUtilities.getShared(cluster, alteredGenes);
                if (shared.size() > 0)
                    vector.add(Boolean.TRUE);
                else
                    vector.add(Boolean.FALSE);
                //vector.add(shared.size());
            }
            //System.out.println(sample + ": " + StringUtils.join(", ", vector));
            sampleToVector.put(sample, vector);
        }
        return sampleToVector;
    }
    
    private Map<String, List<Double>> generateSampleToNetworkModuleVectorsInDegrees(Map<String, Set<String>> sampleToGenes,
                                                                                     String clusterFileName,
                                                                                     int sizeCutoff) throws IOException {
        List<Set<String>> clusters = loadNetworkClusters(clusterFileName);
        return generateSampleToGeneSetVectorsInDegrees(sampleToGenes,
                                                       sizeCutoff, 
                                                       clusters,
                                                       null);
    }
    
    private Map<String, List<Double>> generateSampleToGeneSetVectorsInDegrees(Map<String, Set<String>> sampleToGenes,
                                                                              int sizeCutoff,
                                                                              List<Set<String>> clusters,
                                                                              List<Set<String>> fisInClusters) throws IOException {
        // Filter out clusters based on size
        for (java.util.Iterator<Set<String>> it = clusters.iterator(); it.hasNext();) {
            Set<String> cluster = it.next();
            if (cluster.size() < sizeCutoff)
                it.remove();
        }
        // Pre-calculate degrees for genes in clusters
        List<Map<String, Integer>> geneToDegreesInClusters = new ArrayList<Map<String,Integer>>();
        if (fisInClusters == null) {
            Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
            //Set<String> fis = fu.loadInteractions(R3Constants.RESULT_DIR + "FIsInGene_Pathway_041709.txt");
            for (Set<String> cluster : clusters) {
                Set<String> fisInCluster = InteractionUtilities.getFIs(cluster, fis);
                Map<String, Integer> geneToDegree = InteractionUtilities.generateProteinToDegree(fisInCluster);
                geneToDegreesInClusters.add(geneToDegree);
            }
        }
        else {
            for (int i = 0; i < clusters.size(); i++) {
                Set<String> fisInCluster = fisInClusters.get(i);
                Map<String, Integer> geneToDegree = InteractionUtilities.generateProteinToDegree(fisInCluster);
                geneToDegreesInClusters.add(geneToDegree);
            }
        }
        // Used to normalization
        List<Integer> totalDegrees = new ArrayList<Integer>();
        for (int i = 0; i < clusters.size(); i++) {
            int total = 0;
            Set<String> cluster = clusters.get(i);
            Map<String, Integer> geneToDegrees = geneToDegreesInClusters.get(i);
            for (String gene : geneToDegrees.keySet()) {
                total += geneToDegrees.get(gene);
            }
            totalDegrees.add(total);
        }
        // Want to create a sample to a boolean vector map
        Map<String, List<Double>> sampleToVector = new HashMap<String, List<Double>>();
        for (String sample : sampleToGenes.keySet()) {
            Set<String> alteredGenes = sampleToGenes.get(sample);
            List<Double> vector = new ArrayList<Double>();
            for (int i = 0; i < clusters.size(); i++) {
                Set<String> cluster = clusters.get(i);
                Set<String> shared = InteractionUtilities.getShared(cluster, alteredGenes);
                if (shared.size() == 0)
                    vector.add(0.0);
                else {
                    Map<String, Integer> geneToDegree = geneToDegreesInClusters.get(i);
                    int totalDegree = totalDegrees.get(i);
                    // Need to get the degree
                    int total = 0;
                    for (String gene : shared) {
                        Integer degree = geneToDegree.get(gene);
                        if (degree != null)
                            total += degree;
                    }
                    //vector.add((double)total);
                    vector.add((double) total / totalDegree);
                }
            }
            //System.out.println(sample + ": " + StringUtils.join(", ", vector));
            sampleToVector.put(sample, vector);
        }
        return sampleToVector;
    }
    
    private Map<String, List<Double>> generateSampleToNetworkModuleVectorsInPValue(Map<String, Set<String>> sampleToGenes,
                                                                                   String clusterFileName,
                                                                                   int sizeCutoff,
                                                                                   int totalGense) throws Exception {
        List<Set<String>> clusters = loadNetworkClusters(clusterFileName);
        // Want to create a sample to a boolean vector map
        Map<String, List<Double>> sampleToVector = new HashMap<String, List<Double>>();
        for (String sample : sampleToGenes.keySet()) {
            Set<String> alteredGenes = sampleToGenes.get(sample);
            List<Double> vector = new ArrayList<Double>();
            for (Set<String> cluster : clusters) {
                if (cluster.size() < sizeCutoff)
                    continue;
                Set<String> shared = InteractionUtilities.getShared(cluster, alteredGenes);
                Double pvalue = MathUtilities.calculateHypergeometricPValue(totalGense,
                                                                            alteredGenes.size(),
                                                                            cluster.size(), 
                                                                            shared.size());
                vector.add(pvalue);
            }
            sampleToVector.put(sample, vector);
        }
        return sampleToVector;
    }
    
}
