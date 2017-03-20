/*
 * Created on Dec 4, 2014
 *
 */
package org.reactome.cancer.driver;

import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.log4j.Logger;
import org.junit.Test;
import org.reactome.annotate.AnnotationHelper;
import org.reactome.cancer.MAFFileLoader;
import org.reactome.genome.Transcript;
import org.reactome.r3.UCSCDataAnalyzer;
import org.reactome.r3.UniProtAnalyzer;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

import weka.classifiers.AbstractClassifier;
import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.evaluation.output.prediction.AbstractOutput;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;

/**
 * This class is used to generate a WEKA instances for classifier.
 * @author gwu
 *
 */
public class CancerDriverInstancesGenerator {
//    private final String FI_FILE_NAME = "../FINetworkBuild/results/2015/FIsInGene_031516_BigComp.txt";
//    protected final String FI_FILE_NAME = "../FINetworkBuild/results/2015/FIsInGene_031516_NoUBC.txt";
    protected static final String FI_FILE_NAME = "../FINetworkBuild/results/2015/FIsInGene_081616_Reactome.txt";
//    private final String FI_FILE_NAME = R3Constants.GENE_FI_BIG_COMP_FILE_NAME;
    // Cached FIs for customization on-the-fly
    private Set<String> fis;

    private final static Logger logger = Logger.getLogger(CancerDriverInstancesGenerator.class);
    private NetworkFeature networkFeature;
    private boolean useRandomGenesForNegative; // Use a similar size of random genes as a negative data set.
    // Cache loaded data 
    private boolean cacheData = false;
    // A flag to turn off mutation normalization based on protein lengths
    private boolean noMutationNormalization;
    // For aa change position
    private String aaChangeColumnInMAF;
    private String fiColumnInMAF;
    // MutSig result file
    private String mutSigResultFile;
    private String fiPGMScoreFile;
    
    /**
     * Default constructor.
     */
    public CancerDriverInstancesGenerator() {
    }
    
    public String getFiPGMScoreFile() {
        return fiPGMScoreFile;
    }

    public void setFiPGMScoreFile(String fiPGMScoreFile) {
        this.fiPGMScoreFile = fiPGMScoreFile;
    }

    public String getMutSigResultFile() {
        return mutSigResultFile;
    }

    public void setMutSigResultFile(String mutSigResultFile) {
        this.mutSigResultFile = mutSigResultFile;
    }

    public String getFiColumnInMAF() {
        return fiColumnInMAF;
    }


    public void setFiColumnInMAF(String fiColumnInMAF) {
        this.fiColumnInMAF = fiColumnInMAF;
    }


    public String getAaChangeColumnInMAF() {
        return aaChangeColumnInMAF;
    }


    public void setAaChangeColumnInMAF(String aaChangeColumnInMAF) {
        this.aaChangeColumnInMAF = aaChangeColumnInMAF;
    }


    public boolean isNoMutationNormalization() {
        return noMutationNormalization;
    }

    public void setNoMutationNormalization(boolean noMutationNormalization) {
        this.noMutationNormalization = noMutationNormalization;
    }

    public void setCacheData(boolean cache) {
        this.cacheData = cache;
    }
    
    public void setNetworkFeature(NetworkFeature feature) {
        this.networkFeature = feature;
    }
    
    public NetworkFeature getNetworkFeature() {
        return this.networkFeature;
    }

    public boolean isUseRandomGenesForNegative() {
        return useRandomGenesForNegative;
    }

    public void setUseRandomGenesForNegative(boolean useRandomGenesForNegative) {
        this.useRandomGenesForNegative = useRandomGenesForNegative;
    }
    
    public Instances generateNetworkOnlyInstances(Set<String> driverGenes) throws IOException {
        if (networkFeature == null)
            throw new IllegalStateException("networkFeature has not be set up!");
        ArrayList<Attribute> attributes = new ArrayList<Attribute>();
        // We want to add the first attribute as an ID to later manipulation
        // This is a very weird syntax here. But it is required
        Attribute geneName = new Attribute("Gene", (ArrayList<String>)null);
        attributes.add(geneName);
        
        Attribute networkAtt = new Attribute("Network");
        attributes.add(networkAtt);
        Attribute driverRatio = new Attribute("DriverRatio");
        attributes.add(driverRatio);
        // Class types
        List<String> clsTypes = new ArrayList<String>();
        clsTypes.add(Boolean.TRUE + "");
        clsTypes.add(Boolean.FALSE + "");
        Attribute clsAtt = new Attribute("Driver", clsTypes);
        attributes.add(clsAtt);
        Instances instances = new NetworkMutableInstances("Cancer_Driver_TCGA_BRCA",
                                                          attributes,
                                                          10000);
        instances.setClassIndex(attributes.size() - 1); // The last attribute should be the class attribute.
        // Load genes in the FI network
        Set<String> fis = getFIs();
        Map<String, Set<String>> geneToPartners = InteractionUtilities.generateProteinToPartners(fis);
        if (instances instanceof NetworkMutableInstances) {
            ((NetworkMutableInstances)instances).setGeneToPartners(geneToPartners);
        }
        Map<String, String> geneToType = getGeneToType(fis, driverGenes);
        Map<String, Double> geneToDriverRatio = calculateGeneToDriverRatio(driverGenes, geneToPartners);
        for (String gene : geneToType.keySet()) {
            double[] values = new double[attributes.size()];
            // The first attribute is the gene name
            values[0] = attributes.get(0).addStringValue(gene);
            int index = 1;
            double networkFeature = calculateNetworkFeature(gene, 
                                                            driverGenes,
                                                            geneToPartners,
                                                            geneToDriverRatio);
            values[index++] = networkFeature;
            values[index++] = (double) geneToDriverRatio.get(gene);
            values[index++] = clsAtt.indexOfValue(geneToType.get(gene));
            DenseInstance instance = new DenseInstance(1.0d, values);
            instances.add(instance);
        }
        instances.compactify();
        logger.info("Total instances: " + instances.numInstances() + " with " + instances.numAttributes() + " attributes.");
        return instances;
    }
    
    /**
     * TODO: The implementation of this method has been deleted so that other source code can be compiled.
     * Call this method to generate a WEKA instances on the fly.
     * @return
     * @throws IOException
     * @throws ParseException 
     */
    public Instances generateInstances(Set<String> driverGenes) throws IOException, ParseException {
        return null;
    }

    private Set<String> getFIs() throws IOException {
        if (fis == null || fis.size() == 0)
            return new FileUtility().loadInteractions(FI_FILE_NAME);
        return fis;
    }
    
    public void setFIs(Set<String> fis) {
        this.fis = fis;
    }

    private Set<String> loadPathwayGenes(String pathway) throws IOException {
        AnnotationHelper annotationHelper = new AnnotationHelper();
        annotationHelper.setProteinNameToPathwayFile(R3Constants.RESULT_DIR + "ProteinNameToTopics121514.txt");
        Map<String, Set<String>> proteinToPathways = annotationHelper.loadProteinNameToPathwaysMap();
        Map<String, Set<String>> pathwayToProteins = InteractionUtilities.switchKeyValues(proteinToPathways);
        Set<String> pathwayGenes = pathwayToProteins.get(pathway);
        return pathwayGenes;
    }
    
    private void filterDriversToMutation(Set<String> driverGenes,
                                         Set<String> mutatedGenes) {
        logger.info("Driver genes before filtering to mutation: " + driverGenes.size());
        driverGenes.retainAll(mutatedGenes);
        logger.info("Driver genes after filtering to mutation: " + driverGenes.size());
    }
    
    private String calculateProbability(String gene,
                                        Map<String, Map<String, Float>> sampleToGeneToValue) {
        int count = 0;
        int sampleCount = 0;
        for (String sample : sampleToGeneToValue.keySet()) {
            Map<String, Float> geneToValue = sampleToGeneToValue.get(sample);
            Float fValue = geneToValue.get(gene);
            if (fValue == null)
                continue;
            Integer value = fValue.intValue();
            if (value != 0)
                count ++;
            sampleCount ++;
        }
        if (sampleCount == 0) {
//            logger.warn(gene + " cannot be found in a dat type!");
            //throw new IllegalArgumentException(gene + " cannot be found in a data type!");
            return "na";
        }
        return (double) count / sampleToGeneToValue.size() + "";
    }
    
    double calculateNetworkFeature(String gene,
                                   Set<String> driverGenes,
                                   Map<String, Set<String>> geneToPartners,
                                   Map<String, Double> geneToDriverRatio) {
        switch (networkFeature) {
            case ONE_HOP_SCORE :
                return _calculateNetworkFeature(gene, driverGenes, geneToPartners, true);
            case TWO_HOP_AVERAGE :
                return _calculateNetworkFeature(gene, driverGenes, geneToPartners, false);
            default :
                return _calculateNetworkFeature(gene, 
                                               geneToDriverRatio, 
                                               geneToPartners,
                                               networkFeature);
        }
    }
    
    private double _calculateNetworkFeature(String gene,
                                           Set<String> driverGenes,
                                           Map<String, Set<String>> geneToPartners,
                                           boolean useNeighborOnly) {
        Set<String> partners = getNeighbors(gene, geneToPartners, useNeighborOnly);
        int trueCounter = countTrueDrivers(partners, driverGenes);
        return (double) trueCounter / partners.size();
    }
    
    private Set<String> getNeighbors(String gene,
                                     Map<String, Set<String>> geneToPartners,
                                     boolean useNeighborOnly) {
        if (useNeighborOnly)
            return geneToPartners.get(gene);
        Set<String> partners = geneToPartners.get(gene);
        Set<String> allPartners = new HashSet<String>();
        allPartners.addAll(partners);
        for (String partner : partners) {
            allPartners.addAll(geneToPartners.get(partner));
        }
        allPartners.remove(gene); // Remove itself
        return allPartners;
    }
    
    private double _calculateNetworkFeature(String gene,
                                           Map<String, Double> geneToRatio,
                                           Map<String, Set<String>> geneToPartners,
                                           NetworkFeature feature) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        Double ratio = geneToRatio.get(gene);
        stats.addValue(ratio);
        Set<String> partners = getNeighbors(gene, geneToPartners, true);
        for (String partner : partners) {
            ratio = geneToRatio.get(partner);
            stats.addValue(ratio);
        }
        switch (feature) {
            case ONE_HOP_SCORE_AVERAGE : return stats.getMean();
            case ONE_HOP_SCORE_MAXIMUM : return stats.getMax();
            case ONE_HOP_SCORE_MINIMUM : return stats.getMin();
            case ONE_HOP_SCORE_MEDIAN : return stats.getPercentile(50.0d);
        }
        throw new IllegalArgumentException("NetworkFetaure is not supported: " + feature);
    }
    
    Map<String, Double> calculateGeneToDriverRatio(Set<String> driverGenes,
                                                   Map<String, Set<String>> geneToPartners) {
        Map<String, Double> geneToRatio = new HashMap<String, Double>();
        for (String gene : geneToPartners.keySet()) {
            Set<String> partners = geneToPartners.get(gene);
            int driverCounter = countTrueDrivers(partners, driverGenes);
            Double ratio = (double) driverCounter / partners.size();
            geneToRatio.put(gene, ratio);
        }
        return geneToRatio;
    }
    
    private int countTrueDrivers(Set<String> partners,
                                 Set<String> driverGenes) {
        int trueCounter = 0;
        for (String partner : partners) {
            if (driverGenes.contains(partner))
                trueCounter ++;
        }
        return trueCounter;
    }
    
    private Map<String, Double> loadGeneToMAFIScore(Set<String> hyperMutatedSamples) throws IOException {
        String mafFileName = getMAFFileName();
        MAFFileLoader fileLoader = new MAFFileLoader();
        return fileLoader.loadGeneToMAFIScore(mafFileName,
                                              fiColumnInMAF,
                                              hyperMutatedSamples);
    }
    
    //TODO: Change this to return a non-null file name
    private String getMAFFileName() {
        return null; 
    }
    
    private Map<String, Double> loadGeneToMutationClusterScore(Map<String, Integer> geneToLength,
                                                               Set<String> hyperMutatedSamples) throws IOException {
        String mafFileName = getMAFFileName();
        MAFFileLoader fileLoader = new MAFFileLoader();
        Map<String, int[]> geneToSiteMutFre = fileLoader.loadGeneToSiteMutationFrequency(mafFileName,
                                                                                         aaChangeColumnInMAF,
                                                                                         hyperMutatedSamples,
                                                                                         geneToLength);
        return calculateGeneToMutationClusterScore(geneToSiteMutFre, 5);
    }
    
    private Map<String, Double> loadGeneToMutationClusterScoreViaBP(Set<String> hyperMutatedSamples) throws IOException {
        Map<String, Transcript> geneToTx = new UCSCDataAnalyzer().getHumanGeneToTranscript();
        String mafFileName = getMAFFileName();
        MAFFileLoader fileLoader = new MAFFileLoader();
        Map<String, int[]> geneToSiteMutFre = fileLoader.loadGeneToNucleotideSiteMutationFrequency(mafFileName, hyperMutatedSamples, geneToTx);
        return calculateGeneToMutationClusterScore(geneToSiteMutFre, 15);
    }
    
    private Map<String, Double> loadGeneToFIPGMScore() throws IOException {
        if (fiPGMScoreFile == null)
            return null;
        Map<String, Double> geneToScore = new HashMap<String, Double>();
        FileUtility fu = new FileUtility();
        fu.setInput(fiPGMScoreFile);
        String line = fu.readLine(); // Escape header
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            geneToScore.put(tokens[0], new Double(tokens[1]));
        }
        fu.close();
        return geneToScore;
    }
    
    private Map<String, Double> loadGeneToMutSigScore() throws IOException {
        if (mutSigResultFile == null)
            return null;
        Map<String, Double> geneToPValue = new HashMap<String, Double>();
        FileUtility fu = new FileUtility();
        fu.setInput(mutSigResultFile);
        String line = fu.readLine();
        // The the smallest p-value that is great than 0 so that we can do log-transformation
        Double noZeroSmallest = Double.MAX_VALUE;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // We will use p-values to get a larger range for score
            double pvalue = Double.parseDouble(tokens[tokens.length - 2]);
            geneToPValue.put(tokens[0], pvalue);
            if (pvalue > 0.0 && pvalue < noZeroSmallest)
                noZeroSmallest = pvalue;
        }
        fu.close();
        // Perform log-transform
        Map<String, Double> geneToScore = new HashMap<String, Double>();
        double maxScore = -Math.log(noZeroSmallest / 10.0d);
        for (String gene : geneToPValue.keySet()) {
            Double pvalue = geneToPValue.get(gene);
            if (pvalue.equals(0.0d))
                geneToScore.put(gene, maxScore);
            else
                geneToScore.put(gene, -Math.log(pvalue));
        }
        return geneToScore;
    }
    
    private Set<String> getHyperMutatedSamples() throws IOException {
        // Do a hard-coded filtering to remove hyper-mutated samples
        MAFFileLoader fileLoader = new MAFFileLoader();
        String mafFileName = getMAFFileName();
        Map<String, Integer> sampleToMutationCount = fileLoader.countSampleMutationNumbers(mafFileName);
        Set<String> hyperSamples = new HashSet<String>();
        for (Iterator<String> it = sampleToMutationCount.keySet().iterator(); it.hasNext();) {
            String sample = it.next();
            Integer count = sampleToMutationCount.get(sample);
            if (count > 400) // This is an arbitrary number
                hyperSamples.add(sample);
        }
        return hyperSamples;
    }
    
    private Map<String, Double> loadGeneToMutationRate(Map<String, Integer> geneToLength,
                                                       Set<String> hyperMutatedSamples) throws IOException {
        String mafFileName = getMAFFileName();
        MAFFileLoader fileLoader = new MAFFileLoader();
        Map<String, Set<String>> sampleToGenes = fileLoader.loadSampleToGenes(mafFileName,
                                                                              aaChangeColumnInMAF,
                                                                              null,
                                                                              geneToLength);
        sampleToGenes.keySet().removeAll(hyperMutatedSamples);
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        Map<String, Double> geneToRate = new HashMap<String, Double>();
        for (String gene : geneToSamples.keySet()) {
            Integer geneLength = geneToLength.get(gene);
            if (geneLength == null)
                continue;
            double rate = geneToSamples.get(gene).size() / (double) sampleToGenes.size();
            geneToRate.put(gene, rate / geneLength);
        }
        return geneToRate;
    }
    
    @Test
    public void testLoadGeneToMutationClusterScore() throws IOException {
        Map<String, Integer> geneToLength = new UniProtAnalyzer().loadGeneToProteinLength();
        Set<String> hyperMutSamples = getHyperMutatedSamples();
//        Map<String, Double> geneToClusterScore = loadGeneToMutationClusterScore(geneToLength,
//                                                                                hyperMutSamples);
        Map<String, Double> geneToClusterScore = loadGeneToMutationClusterScoreViaBP(hyperMutSamples);
        String[] genes = new String[]{
                "TP53",
                "CASP8",
                "BRCA2",
                "CHEK2",
                "BAP1",
                "PALB2",
                "PBRM1",
                "BRIP1",
                "KRAS",
                "NCOR1",
                "TBL1XR1",
                "SF3B1",
                "ARID1A",
                "NF1",
                "VEZF1"
        };
        for (String gene : genes) {
            System.out.println(gene + "\t" + 
                               geneToLength.get(gene) + "\t" +
                               geneToClusterScore.get(gene));
        }
    }
    
    /**
     * This method is used to calculate a score for mutation clustering for each gene.
     * @param geneToSiteMutationFrequency
     * @return
     */
    private Map<String, Double> calculateGeneToMutationClusterScore(Map<String, int[]> geneToSiteMutationFrequency,
                                                                    int windowSize) {
        Map<String, Double> geneToScore = new HashMap<String, Double>();
        // A list of window sizes
        for (String gene : geneToSiteMutationFrequency.keySet()) {
            int[] siteMutationFreq = geneToSiteMutationFrequency.get(gene);
            // Get the total mutations
            int totalMutations = 0;
            for (int fre : siteMutationFreq)
                totalMutations += fre;
            double totalScore = 0.0d;
            for (int i = 0; i < siteMutationFreq.length; i++) {
                if (i + windowSize > siteMutationFreq.length)
                    break;
                int windowTotal = 0;
                for (int j = i; j < i + windowSize; j++) {
                    windowTotal += siteMutationFreq[j];
                }
                // Assume the minimum is 2 in order to make sure clustering is meaningful
                if (windowTotal < 2)
                    continue;
                double prior = 1.0d * windowSize / siteMutationFreq.length;
                double windowRatio = (double) windowTotal / totalMutations;
                if (windowRatio < prior)
                    continue; // Should not be counted as interesting
                double pvalue = MathUtilities.calculateBinomialPValue(prior, totalMutations, windowTotal);
                if (pvalue > 0.10d) // Hard-coded: following the selection from OncoDriverClust paper.
                    continue;
                double oddsRatio = windowRatio / prior;
                // Since we want to check clustering effect only, keep the positive odds ratio
                if (oddsRatio > 1.0d)
                    totalScore += Math.log(oddsRatio);
            }
            geneToScore.put(gene, totalScore / siteMutationFreq.length);
        }
        return geneToScore;
    }

    private Map<String, String> getGeneToType(Set<String> fis,
                                              Set<String> driverGenes) throws IOException {
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        logger.info("Total genes in the FI network: " + fiGenes.size());
        // Only want to work on drivers genes in the FI network
        driverGenes.retainAll(fiGenes);
        // Don't consider these known driver genes
        fiGenes.removeAll(driverGenes);
        logger.info("Driver genes in the FI network: " + driverGenes.size());
        logger.info("FI genes excluding the driver genes: " + fiGenes.size());
        
        Set<String> negativeGenes = null;
        if (useRandomGenesForNegative) {
            // Get the same number of random genes
            negativeGenes = MathUtilities.randomSampling(fiGenes, driverGenes.size() + 20);
        }
        else
            negativeGenes = new HashSet<String>(fiGenes); // Use all FI genes as the negative set
        
        Map<String, String> geneToType = new HashMap<String, String>();
        for (String gene : driverGenes)
            geneToType.put(gene, Boolean.TRUE.toString());
        for (String gene : negativeGenes)
            geneToType.put(gene, Boolean.FALSE.toString());
        return geneToType;
    }
    
    class NetworkMutableInstances extends Instances {
        // Want to keep the map from genes and their partners to update
        // network features during cross-validation.
        private Map<String, Set<String>> geneToPartners;
        // Record the network attribute
        private Attribute networkAtt;
        // Attribute for the ratio
        private Attribute driverRatioAtt;
        
        public NetworkMutableInstances(String name,
                                       ArrayList<Attribute> attributes,
                                       int size) {
            super(name, attributes, size);
            checkNetworkAttributes(attributes);
        }

        private void checkNetworkAttributes(ArrayList<Attribute> attributes) {
            for (Attribute att : attributes) {
                if (att.name().equalsIgnoreCase("network"))
                    networkAtt = att;
                else if (att.name().equalsIgnoreCase("DriverRatio"))
                    driverRatioAtt = att;
            }
        }

        public NetworkMutableInstances(Instances instances) {
            super(instances);
            ArrayList<Attribute> attributes = new ArrayList<Attribute>();
            for (int i = 0; i < instances.numAttributes(); i++) {
                Attribute att = instances.attribute(i);
                attributes.add(att);
            }
            checkNetworkAttributes(attributes);
        }
        
        public void setGeneToPartners(Map<String, Set<String>> geneToPartners) {
            this.geneToPartners = geneToPartners;
        }
        
        public Map<String, Set<String>> getGeneToPartners() {
            return this.geneToPartners;
        }

        @Override
        public Instances testCV(int numFolds, int numFold) {
            Instances rtn = super.testCV(numFolds, numFold);
            resetNetworkFeaturesForTest(rtn);
            return rtn;
        }
        
        private void resetNetworkFeaturesForTest(Instances subset) {
            if (geneToPartners == null || networkAtt == null)
                return; // Nothing to be done
            // We want to find driver genes in the train set
            Set<Instance> allInstanes = getInstances(this);
//            logger.info("All instances: " + allInstanes.size());
            Set<Instance> instances = getInstances(subset);
//            logger.info("Test instances: " + instances.size());
            allInstanes.removeAll(instances);
//            logger.info("Train instances: " + allInstanes.size());
            Set<String> driverGenes = getDriverGenes(allInstanes);
            resetNetworkFeatures(subset,
                                 driverGenes);
        }
        
        private Set<Instance> getInstances(Instances data) {
            Set<Instance> instances = new HashSet<Instance>();
            for (int i = 0; i < data.numInstances(); i++)
                instances.add(data.instance(i));
            return instances;
        }
        
        private void resetNetworkFeaturesForTrain(Instances subset) {
            if (geneToPartners == null || networkAtt == null)
                return;
//            logger.info("resetNewtorkFeature for an Instances: " + subset.numInstances() + " instances...");
            // Get all driver genes in this subset
            Set<Instance> instances = getInstances(subset);
            Set<String> driverGenes = getDriverGenes(instances);
            resetNetworkFeatures(subset, 
                                driverGenes);
        }

        private void resetNetworkFeatures(Instances subset,
                                         Set<String> driverGenes) {
            // Update network feature
            Map<String, Double> geneToRatio = calculateGeneToDriverRatio(driverGenes,
                                                                         geneToPartners);
            for (int i = 0; i < subset.numInstances(); i++) {
                Instance instance = subset.instance(i);
                // Get the gene name: it should be in the first attribute
                String geneName = instance.stringValue(0);
                double networkFeature = calculateNetworkFeature(geneName, 
                                                                driverGenes, 
                                                                geneToPartners,
                                                                geneToRatio);
                instance.setValue(networkAtt, networkFeature);
                instance.setValue(driverRatioAtt, geneToRatio.get(geneName));
            }
        }

        private Set<String> getDriverGenes(Set<Instance> instances) {
            Set<String> driverGenes = new HashSet<String>();
            for (Instance instance : instances) {
                int clsIndex = instance.classIndex();
                String cls = instance.stringValue(clsIndex);
                if (cls.equals(Boolean.TRUE.toString())) {
                    // Get the gene symbol
                    String geneName = instance.stringValue(0);
                    driverGenes.add(geneName);
                }
            }
            return driverGenes;
        }

        @Override
        public Instances trainCV(int numFolds, int numFold, Random random) {
            Instances rtn = super.trainCV(numFolds, numFold, random);
            resetNetworkFeaturesForTrain(rtn);
            return rtn;
        }
    }
    
    /**
     * We need a customized Evaluation so that we can use our own mutable network feature
     * in a customized Instances.
     * @author gwu
     *
     */
    class MutableDataEvalulation extends Evaluation {
        
        public MutableDataEvalulation(Instances data) throws Exception {
            super(data);
        }
        
        /**
         * Basically this is copied from the original implementation by commenting out the first line:
         * data = new Instances(data) to avoid copy, so that our implementation can be used.
         */
        @Override
        public void crossValidateModel(Classifier classifier, 
                                       Instances data,
                                       int numFolds, 
                                       Random random, 
                                       Object... forPredictionsPrinting) throws Exception {
            // Make a copy of the data we can reorder
            if (data instanceof NetworkMutableInstances) {
                Map<String, Set<String>> geneToPartners = ((NetworkMutableInstances)data).getGeneToPartners();
                // Make a new copy
                data = new NetworkMutableInstances(data);
                ((NetworkMutableInstances)data).setGeneToPartners(geneToPartners);
            }
            else
                data = new Instances(data);
            data.randomize(random);
            if (data.classAttribute().isNominal()) {
                data.stratify(numFolds);
            }
            
            // We assume that the first element is a
            // weka.classifiers.evaluation.output.prediction.AbstractOutput object
            AbstractOutput classificationOutput = null;
            if (forPredictionsPrinting.length > 0) {
                // print the header first
                classificationOutput = (AbstractOutput) forPredictionsPrinting[0];
                classificationOutput.setHeader(data);
                classificationOutput.printHeader();
            }
            
            // Do the folds
            for (int i = 0; i < numFolds; i++) {
                Instances train = data.trainCV(numFolds, i, random);
                setPriors(train);
                Classifier copiedClassifier = AbstractClassifier.makeCopy(classifier);
                copiedClassifier.buildClassifier(train);
                Instances test = data.testCV(numFolds, i);
                evaluateModel(copiedClassifier, test, forPredictionsPrinting);
            }
            m_NumFolds = numFolds;
            
            if (classificationOutput != null)
                classificationOutput.printFooter();
        }
    }
    
    static enum NetworkFeature {
        ONE_HOP_SCORE, // ratio of driver genes to one-hop partners
        TWO_HOP_AVERAGE, // ratio of driver genes to two-hop partners
        ONE_HOP_SCORE_AVERAGE, // average of one hop scores for one-hop neighbor genes
        ONE_HOP_SCORE_MAXIMUM, // maximum of one hop scores for one-hop neighbor genes
        ONE_HOP_SCORE_MEDIAN, // medium of one hop scores for one-hop neighbor genes
        ONE_HOP_SCORE_MINIMUM // minimum of one hop scores for one-hop neighbor genes
    }
    
}
