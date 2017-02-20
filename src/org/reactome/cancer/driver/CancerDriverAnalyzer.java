/*
 * Created on Nov 18, 2014
 *
 */
package org.reactome.cancer.driver;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.log4j.Logger;
import org.junit.Test;
import org.reactome.cancer.MATFileLoader;
import org.reactome.cancer.NetworkClusterAnalyzer;
import org.reactome.cancer.driver.CancerDriverInstancesGenerator.NetworkFeature;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

import weka.classifiers.AbstractClassifier;
import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.functions.Logistic;
import weka.classifiers.meta.FilteredClassifier;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Utils;
import weka.core.converters.ConverterUtils.DataSink;
import weka.filters.unsupervised.attribute.Remove;

/**
 * This class is used to analyze cancer driver genes.
 * @author gwu
 *
 */
public class CancerDriverAnalyzer {
    private static Logger logger = Logger.getLogger(CancerDriverAnalyzer.class);
//    public static final String RESULT_DIR = "results/DriverGenes/";
    public static final String RESULT_DIR = "datasets/ICGC/2016_04/Drivers/";
    
    private FileUtility fu = new FileUtility();
    
    public CancerDriverAnalyzer() {
    }
    
    @Test
    public void runNetworkFeaturesPermuation() throws Exception {
        Set<String> fis = new FileUtility().loadInteractions(CancerDriverInstancesGenerator.FI_FILE_NAME);
        Map<String, Set<String>> geneToPartners = InteractionUtilities.generateProteinToPartners(fis);
        System.out.println("Total genes: " + geneToPartners.size());
        
//        Set<String> drivers = getDriverGenes(null);
        Set<String> drivers = new FICancerDriverPredictor().loadCancerCensusGenes();
        
        System.out.println("Total driver genes: " + drivers.size());
        drivers.retainAll(geneToPartners.keySet());
        System.out.println("Filtered to the FI network: " + drivers.size());
        
        CancerDriverInstancesGenerator generator = getInstancesGenerator();
        
        List<Double> pvalues = performTTestForNetworkFeatures(geneToPartners, 
                                                              drivers,
                                                              generator);
        System.out.println("Data\tDriverRatio\tNetworkFeature");
        System.out.println("RealData\t" + pvalues.get(0) + "\t" + pvalues.get(1));
        int permutation = 1000;
        for (int i = 0; i < permutation; i++) {
            Set<String> randomFIs = InteractionUtilities.generateRandomFIsViaConfigurationModel(fis);
            Map<String, Set<String>> randomGeneToPartners = InteractionUtilities.generateProteinToPartners(randomFIs);
            pvalues = performTTestForNetworkFeatures(randomGeneToPartners, drivers, generator);
            System.out.println("Random" + i + "\t" + 
                    pvalues.get(0) + "\t" +
                    pvalues.get(1));
        }
    }

    private List<Double> performTTestForNetworkFeatures(Map<String, Set<String>> geneToPartners, 
                                                Set<String> drivers,
                                                CancerDriverInstancesGenerator generator) throws MathException {
        Map<String, Double> geneToDriverRatio = generator.calculateGeneToDriverRatio(drivers, 
                                                                                     geneToPartners);
        List<Double> driverValues1 = new ArrayList<Double>();
        List<Double> nonDriverValues1 = new ArrayList<Double>();
        List<Double> driverValues2 = new ArrayList<Double>();
        List<Double> nonDriverValues2 = new ArrayList<Double>();
        for (String gene : geneToDriverRatio.keySet()) {
            Double driverRatio = geneToDriverRatio.get(gene);
            Double networkFeature = generator.calculateNetworkFeature(gene, drivers, geneToPartners, geneToDriverRatio);
            if (drivers.contains(gene)) {
                driverValues1.add(driverRatio);
                driverValues2.add(networkFeature);
            }
            else {
                nonDriverValues1.add(driverRatio);
                nonDriverValues2.add(networkFeature);
            }
        }
        Double pvalue1 = MathUtilities.calculateTTest(driverValues1, nonDriverValues1);
        Double pvalue2 = MathUtilities.calculateTTest(driverValues2, nonDriverValues2);
        List<Double> rtn = new ArrayList<Double>();
        rtn.add(pvalue1);
        rtn.add(pvalue2);
        return rtn;
    }
    
    @Test
    public void runDriverGeneDistributionPermutation() throws Exception {
        Set<String> fis = new FileUtility().loadInteractions(CancerDriverInstancesGenerator.FI_FILE_NAME);
        Map<String, Set<String>> geneToPartners = InteractionUtilities.generateProteinToPartners(fis);
        System.out.println("Total genes: " + geneToPartners.size());
        
//        Set<String> drivers = getDriverGenes(null);
        Set<String> drivers = new FICancerDriverPredictor().loadCancerCensusGenes();
        
        System.out.println("Total driver genes: " + drivers.size());
        drivers.retainAll(geneToPartners.keySet());
        System.out.println("Filtered to the FI network: " + drivers.size());
        double ratio = (double) drivers.size() / geneToPartners.size();
        System.out.println("Ratio of drivers to all genes: " + ratio);
        System.out.println("Data\tFirstLevelTTest\tSecondLevelTTest");
        List<Double> pvalues = calculateTTest(geneToPartners, drivers, ratio);
        System.out.println("RealData\t" + pvalues.get(0) + "\t" + pvalues.get(1));
        int permutation = 1000;
        for (int i = 0; i < permutation; i++) {
            Set<String> randomFIs = InteractionUtilities.generateRandomFIsViaConfigurationModel(fis);
            Map<String, Set<String>> randomGeneToPartners = InteractionUtilities.generateProteinToPartners(randomFIs);
            pvalues = calculateTTest(randomGeneToPartners, drivers, ratio);
            System.out.println("Random" + i + "\t" + 
                               pvalues.get(0) + "\t" + 
                               pvalues.get(1));
        }
    }
    
    private List<Double> calculateTTest(Map<String, Set<String>> geneToPartners,
                                        Set<String> drivers,
                                        double ratio) throws Exception {
        List<Double> pvalueList1 = new ArrayList<Double>();
        List<Double> pvalueList2 = new ArrayList<Double>();
        List<Double> pvalueListL1 = new ArrayList<Double>();
        List<Double> pvalueListL2 = new ArrayList<Double>();
        for (String gene : geneToPartners.keySet()) {
            Set<String> parteners = getNeighbors(gene, geneToPartners, 1);
            if (parteners == null || parteners.size() == 0)
                continue; // Escape if no interaction 
            Set<String> shared = InteractionUtilities.getShared(parteners, drivers);
            Double pvalue = MathUtilities.calculateBinomialPValue(ratio, 
                                                                  parteners.size(), 
                                                                  shared.size());
            Set<String> partenersL2 = getNeighbors(gene, geneToPartners, 2);
            Set<String> sharedL2 = InteractionUtilities.getShared(partenersL2, drivers);
            Double pvalueL2 = MathUtilities.calculateBinomialPValue(ratio,
                                                                    partenersL2.size(),
                                                                    sharedL2.size());
            int isDriver = drivers.contains(gene) ? 1 : 0;
            if (isDriver == 1) {
                pvalueList1.add(Math.log10(pvalue));
                pvalueListL1.add(Math.log10(pvalueL2));
            }
            else {
                pvalueList2.add(Math.log10(pvalue));
                pvalueListL2.add(Math.log10(pvalueL2));
            }
        }
        Double tTest1 = MathUtilities.calculateTTest(pvalueList1, pvalueList2);
        Double tTest2 = MathUtilities.calculateTTest(pvalueListL1, pvalueListL2);
        List<Double> rtn = new ArrayList<Double>();
        rtn.add(tTest1);
        rtn.add(tTest2);
        return rtn;
    }
    
    @Test
    public void analyzeDriverGeneDistribution() throws Exception {
        Set<String> fis = new FileUtility().loadInteractions(CancerDriverInstancesGenerator.FI_FILE_NAME);
        Map<String, Set<String>> geneToPartners = InteractionUtilities.generateProteinToPartners(fis);
        System.out.println("Total genes: " + geneToPartners.size());
        
//        Set<String> drivers = getDriverGenes(null);
        Set<String> drivers = new FICancerDriverPredictor().loadCancerCensusGenes();
        
        System.out.println("Total driver genes: " + drivers.size());
        drivers.retainAll(geneToPartners.keySet());
        System.out.println("Filtered to the FI network: " + drivers.size());
        double ratio = (double) drivers.size() / geneToPartners.size();
        System.out.println("Ratio of drivers to all genes: " + ratio);
        System.out.println("Gene\tPartners\tDrivers\tPValue\tParntersL2\tDriversL2\tPValueL2\tIsDriver");
        List<Double> pvalueList1 = new ArrayList<Double>();
        List<Double> pvalueList2 = new ArrayList<Double>();
        List<Double> pvalueListL1 = new ArrayList<Double>();
        List<Double> pvalueListL2 = new ArrayList<Double>();
        for (String gene : geneToPartners.keySet()) {
            Set<String> parteners = getNeighbors(gene, geneToPartners, 1);
            Set<String> shared = InteractionUtilities.getShared(parteners, drivers);
            Double pvalue = MathUtilities.calculateBinomialPValue(ratio, 
                                                                  parteners.size(), 
                                                                  shared.size());
            Set<String> partenersL2 = getNeighbors(gene, geneToPartners, 2);
            Set<String> sharedL2 = InteractionUtilities.getShared(partenersL2, drivers);
            Double pvalueL2 = MathUtilities.calculateBinomialPValue(ratio,
                                                                    partenersL2.size(),
                                                                    sharedL2.size());
            int isDriver = drivers.contains(gene) ? 1 : 0;
            System.out.println(gene + "\t" + parteners.size() + "\t" + shared.size() + "\t" + pvalue + "\t" +
                               partenersL2.size() + "\t" + sharedL2.size() + "\t" + pvalueL2 + "\t" + isDriver);
            if (isDriver == 1) {
                pvalueList1.add(Math.log10(pvalue));
                pvalueListL1.add(Math.log10(pvalueL2));
            }
            else {
                pvalueList2.add(Math.log10(pvalue));
                pvalueListL2.add(Math.log10(pvalueL2));
            }
        }
        Double tTest1 = MathUtilities.calculateTTest(pvalueList1, pvalueList2);
        System.out.println("Level 1 t-test for p-values: " + tTest1);
        System.out.println("Average: " + MathUtilities.calculateMean(pvalueList1) + ", " + MathUtilities.calculateMean(pvalueList2));
        Double tTest2 = MathUtilities.calculateTTest(pvalueListL1, pvalueListL2);
        System.out.println("Level 2 t-test for p-values: " + tTest2);
        System.out.println("Average: " + MathUtilities.calculateMean(pvalueListL1) + ", " + MathUtilities.calculateMean(pvalueListL2));
    }
    
    private Set<String> getNeighbors(String gene,
                                     Map<String, Set<String>> geneToPartners,
                                     int level) {
        Set<String> neighbor = new HashSet<String>();
        int tempLevel = 0;
        neighbor.add(gene);
        while (tempLevel < level) {
            tempLevel ++;
            Set<String> tmpSet = new HashSet<String>(neighbor);
            for (String tmp : tmpSet) {
                neighbor.addAll(geneToPartners.get(tmp));
            }
        }
        neighbor.remove(gene); // Remove itself
        return neighbor;
    }
    
    @Test
    public void generatePredictedScores() throws Exception {
        FileUtility.initializeLogging();
        FICancerDriverPredictor predictor = new FICancerDriverPredictor();
        
        //              Set<String> trainDrivers = getDriverGenes("ovarian");
        //        Set<String> trainDrivers = getDriverGenes("colorectal");
//        Set<String> trainDrivers = getDriverGenes("breast");
        Set<String> trainDrivers = predictor.loadCancerCensusGenes();
        //      Set<String> testDriverGenes = predictor.loadScienceReportDriverGenes(true);
        System.out.println("Train drivers: " + trainDrivers.size());
//        Set<String> cgcTrainer = predictor.loadCancerCensusGenes("breast");
//        Set<String> natureTrainer = predictor.loadNatureDriverGenes("BRCA");
//        Set<String> natureTrainer = predictor.loadNatureGADDriverGenes("BRCA");
        // Load test genes
        Set<String> natureGenes = predictor.loadNatureDriverGenes();
        Set<String> cgcGenes = predictor.loadCancerCensusGenes();
        Set<String> testDrivers = new HashSet<String>();
        testDrivers.addAll(natureGenes);
        testDrivers.addAll(cgcGenes);
        //testDrivers.addAll(predictor.loadScienceReportDriverGenes(false));
        testDrivers.addAll(predictor.loadSRHCDriverNoFIGenes());
        testDrivers.removeAll(trainDrivers);
        System.out.println("Test Drivers: " + testDrivers.size());
        trainDrivers.addAll(testDrivers);
        
        // Output
//        String outputFileName = RESULT_DIR + "TCGABreastCancerScores_050115.txt";
//        String outputFileName = RESULT_DIR + "TCGABreastCancerScores_051915_1.txt";
//        String outputFileName = RESULT_DIR + "TCGABreastCancerScores_052115.txt";
//        String outputFileName = RESULT_DIR + "ICGCPanCancerScores_052115.txt";
//        String outputFileName = RESULT_DIR + "ICGCPanCancerScores_NoMutation_052815.txt";
//        String outputFileName = RESULT_DIR + "ICGCPanCancerScores_Network_Only_061115.txt";
//        String outputFileName = RESULT_DIR + "ICGCPanCancerScores_NoMutSig_061115.txt";
//        String outputFileName = "results/ICGC_PanCancer/ICGCPanCancerScores_NoMutation_052615.txt";
//        String outputFileName = RESULT_DIR + "TCGAPanCancerScores_052115.txt";
        String outputFileName = RESULT_DIR + "TCGAPanCancerScores_081116_2.txt";
        fu.setOutput(outputFileName);
        
        CancerDriverInstancesGenerator generator = getInstancesGenerator();
        
        Remove remove = new Remove();
//        remove.setAttributeIndices("1,2");
        remove.setAttributeIndices("1");
//        remove.setAttributeIndices("1,2,3,4,7"); // Network only features left
        
        Instances instances = predictor.calculateClassifierScore(trainDrivers, generator, remove);
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < instances.numAttributes(); i++)
            builder.append(instances.attribute(i).name()).append("\t");
        builder.append("TrainDrivers\tOtherDrivers");
        // TrainDrivers to indicate where it is from: 0: none, 1: nature, 2: cgc, 3: both
//        System.out.println(builder.toString());
        fu.printLine(builder.toString());
        for (Instance inst : instances) {
            builder.setLength(0);
            String gene = inst.stringValue(0);
            builder.append(gene).append("\t");
            for (int i = 1; i < instances.numAttributes(); i++)
                builder.append(inst.value(i)).append("\t");
            int train = 0;
            if (trainDrivers.contains(gene))
                train = 1;
//            if (natureTrainer.contains(gene))
//                train += 1;
//            if (cgcTrainer.contains(gene))
//                train += 2;
            builder.append(train).append("\t");
            builder.append((testDrivers.contains(gene) ? "1" : "0"));
//            System.out.println(builder.toString());
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    @Test
    public void generateSampleToDriversFileForCOADREAD() throws Exception {
        // Mutation file
        String mutationFile = "test_data/tcga_coadread/coadread.maf.txt";
        MATFileLoader loader = new MATFileLoader();
        Map<String, Set<String>> sampleToGenes = loader.loadSampleToGenes(mutationFile, false);
        NetworkClusterAnalyzer analyzer = new NetworkClusterAnalyzer();
        
        String driverGeneFile = RESULT_DIR + "TCGA_COADREAD_PreictedDriver_Known_001_031815.txt";
        
        String outputFile = RESULT_DIR + "Sample_TCGA_COADREAD_PreictedDriver_Known_001_031815.txt";
        
        Set<String> driverGenes = fu.loadInteractions(driverGeneFile);
        System.out.println("Total drivers: " + driverGenes.size());
        analyzer.generateSampleToSelectedGenesMatrix(sampleToGenes,
                                                     driverGenes, 
                                                     outputFile);
        
        // We want to compare prediction results with known drivers
//        driverGenes = getDriverGenes("colorectal");
//        System.out.println("Known drivers: " + driverGenes.size());
////        outputFile = RESULT_DIR + "Sample_COADREAD_Known_Drivers_022615.txt";
//        outputFile = RESULT_DIR + "Sample_COADREAD_Known_Drivers_030415.txt";      
//        analyzer.generateSampleToSelectedGenesMatrix(sampleToGenes,
//                                                     driverGenes, 
//                                                     outputFile);
    }
    
    @Test
    public void generateSampleToDriversFileForOV() throws Exception {
        // Mutation file
        String mutationFile = "test_data/tcga_ov/ov.maf.txt";
        MATFileLoader loader = new MATFileLoader();
        Map<String, Set<String>> sampleToGenes = loader.loadSampleToGenes(mutationFile, false);
        NetworkClusterAnalyzer analyzer = new NetworkClusterAnalyzer();
        
        String driverGeneFile = RESULT_DIR + "TCGA_OV_PreictedDriver_Known_00026_022615.txt";
        Set<String> driverGenes = fu.loadInteractions(driverGeneFile);
        String outputFile = RESULT_DIR + "Sample_TCGA_OV_PreictedDriver_Known_00026_022615.txt";
        
        System.out.println("Total drivers: " + driverGenes.size());
        analyzer.generateSampleToSelectedGenesMatrix(sampleToGenes,
                                                     driverGenes, 
                                                     outputFile);
        
        // We want to compare prediction results with known drivers
//        driverGenes = getDriverGenes("ovarian");
//        outputFile = RESULT_DIR + "Sample_TCGA_OV_Known_Drivers_022615.txt";        
//        analyzer.generateSampleToSelectedGenesMatrix(sampleToGenes,
//                                                     driverGenes, 
//                                                     outputFile);
    }
    
    @Test
    public void generateSampleToDriversFileForBRCA() throws Exception {
        // Mutation file
        String mutationFile = "test_data/tcga_brca/brca.maf.txt";
        MATFileLoader loader = new MATFileLoader();
        Map<String, Set<String>> sampleToGenes = loader.loadSampleToGenes(mutationFile, false);
        NetworkClusterAnalyzer analyzer = new NetworkClusterAnalyzer();
        
//        String driverGeneFile = RESULT_DIR + "TCGA_BRCA_PreictedDriver_Known_001_022615.txt";
//        String outputFile = RESULT_DIR + "Sample_TCGA_BRCA_PreictedDriver_Known_001_022615.txt";
        
        String driverGeneFile = RESULT_DIR + "TCGA_BRCA_PreictedDriver_Known_0009_031815_1.txt";
        String outputFile = RESULT_DIR + "Sample_TCGA_BRCA_PreictedDriver_Known_0009_031815_1.txt";
        
        Set<String> driverGenes = fu.loadInteractions(driverGeneFile);
        System.out.println("Total drivers: " + driverGenes.size());
        analyzer.generateSampleToSelectedGenesMatrix(sampleToGenes,
                                                     driverGenes, 
                                                     outputFile);
        
        // We want to compare prediction results with known drivers
//        driverGenes = getDriverGenes("ovarian");
//        outputFile = RESULT_DIR + "Sample_TCGA_OV_Known_Drivers_022615.txt";        
//        analyzer.generateSampleToSelectedGenesMatrix(sampleToGenes,
//                                                     driverGenes, 
//                                                     outputFile);
    }
    
    @Test
    public void checkPredictedDriverGenes() throws Exception {
        FileUtility.initializeLogging();
        FICancerDriverPredictor predictor = new FICancerDriverPredictor();
        
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
//        Set<String> trainDrivers = getDriverGenes("ovarian");
//        Set<String> trainDrivers = getDriverGenes("colorectal");
//        Set<String> trainDrivers = getDriverGenes("breast");
        Set<String> trainDrivers = predictor.loadNatureDriverGenes();
//        Set<String> testDriverGenes = predictor.loadScienceReportDriverGenes(true);
        Set<String> testDrivers = predictor.loadNatureDriverGenes();
        testDrivers.addAll(predictor.loadCancerCensusGenes());
        //testDrivers.addAll(predictor.loadScienceReportDriverGenes(false));
//        testDrivers.addAll(predictor.loadSRHCDriverNoFIGenes());
        testDrivers.removeAll(trainDrivers);
        System.out.println("Test Drivers: " + testDrivers.size());
        testDrivers.retainAll(fiGenes);
        System.out.println("in FI: " + testDrivers.size());
        System.out.println("Total FI genes: " + fiGenes.size());
        fiGenes.removeAll(trainDrivers);
        System.out.println("Remove train Drivers: " + fiGenes.size());
        double originalRatio = (double) testDrivers.size() / fiGenes.size();
        System.out.println("Ratio: " + originalRatio);
        
        // Try different threshold values: values for breast cancer
//        String cancerType = "BRCA";
//        double[] thresholds = new double[] {
//                0.01,
//                0.02,
//                0.03,
//                0.04,
//                0.05,
//                0.06
//        };
//        thresholds = new double[]{0.009};
//        
//        // Values for ovarian
//        String cancerType = "OV";
//        double[] thresholds = new double[]{0.002, 0.0023, 0.0026, 0.0029, 0.003};
//        thresholds = new double[]{0.0026};
        
//        String cancerType = "COADREAD";
//        double[] thresholds = new double[]{0.008, 0.009, 0.01, 0.015, 0.02};
//        thresholds = new double[]{0.01};
        
        // All driver genes
        String cancerType = "All";
        double[] thresholds = new double[] {
                0.1,
                0.2,
                0.3,
                0.4,
                0.5,
                0.6
        };
        
        List<Set<String>> geneSets = predictor.getPredictedDriverGenesInBatch(trainDrivers, 
                                                                              thresholds);
        for (int i = 0; i < thresholds.length; i++) {
            Set<String> predicted = geneSets.get(i);
            System.out.println("Threshold " + thresholds[i] + ": " + predicted.size());
            Set<String> sharedWithNature = InteractionUtilities.getShared(predicted, trainDrivers);
            System.out.println("Shared with train: " + sharedWithNature.size());
            predicted.removeAll(trainDrivers);
            System.out.println("Remove train: " + predicted.size());
//            System.out.println("Genes: " + predicted);
            Set<String> predictedInSrAll = InteractionUtilities.getShared(testDrivers, predicted);
            System.out.println("In test: " + predictedInSrAll.size());
//            System.out.println("Genes: " + predictedInSrAll);
            double ratio = (double) predictedInSrAll.size() / predicted.size();
            System.out.println("Ratio: " + ratio);
            double pvalue = MathUtilities.calculateBinomialPValue(originalRatio, predicted.size(), predictedInSrAll.size());
            System.out.println("Pvalue: " + pvalue);
            System.out.println();
            if (thresholds.length > 1)
                continue;
            String tmp = thresholds[i] + "";
            String fileName = "TCGA_" + cancerType + "_PreictedDriver_Known_" + tmp.replace(".", "") + "_031815.txt";
//            String fileName = "TCGA_" + cancerType + "_PreictedDriver_Known_" + tmp.replace(".", "") + "_022615.txt";
            Set<String> set = new HashSet<String>(geneSets.get(i));
            set.addAll(trainDrivers);
            List<String> list = new ArrayList<String>(set);
            Collections.sort(list);
            fu.saveCollection(list,
                              RESULT_DIR + fileName);
        }
    }
    
    Set<String> getDriverGenes(String cancerType) throws IOException {
        FICancerDriverPredictor helper = new FICancerDriverPredictor();
        Set<String> allDriverGenes = new HashSet<String>();
        if (cancerType == null) {
            allDriverGenes.addAll(helper.loadNatureDriverGenes());
            //      driverGenes = helper.loadScienceReportDriverGenes(false);
            allDriverGenes.addAll(helper.loadSRHCDriverNoFIGenes());
            allDriverGenes.addAll(helper.loadCancerCensusGenes());
        }
        else if (cancerType.equals("breast")) {
//            Set<String> driverGenes = helper.loadNatureDriverGenes("BRCA");
            Set<String> driverGenes = helper.loadNatureGADDriverGenes("BRCA");
            logger.info("Known drivers from Nautre: " + driverGenes.size());
            allDriverGenes.addAll(driverGenes);
            driverGenes = helper.loadCancerCensusGenes(cancerType);
            logger.info("Known drivers from CGC: " + driverGenes.size());
            allDriverGenes.addAll(driverGenes);
            logger.info("Total: " + allDriverGenes.size());
        }
        else if (cancerType.equals("ovarian")) {
            Set<String> driverGenes = helper.loadNatureDriverGenes("OV");
            logger.info("Known drivers from Nautre: " + driverGenes.size());
            allDriverGenes.addAll(driverGenes);
            driverGenes = helper.loadCancerCensusGenes(cancerType);
            logger.info("Known drivers from CGC: " + driverGenes.size());
            allDriverGenes.addAll(driverGenes);
            logger.info("Total: " + allDriverGenes.size());
        }
        else if (cancerType.equals("colorectal")) {
            Set<String> driverGenes = helper.loadNatureDriverGenes("COADREAD");
            logger.info("Known drivers from Nautre: " + driverGenes.size());
            allDriverGenes.addAll(driverGenes);
            driverGenes = helper.loadCancerCensusGenes(cancerType);
            logger.info("Known drivers from CGC: " + driverGenes.size());
            allDriverGenes.addAll(driverGenes);
            logger.info("Total: " + allDriverGenes.size());
        }
        return allDriverGenes;
    }
    
    @Test
    public void checkOverlapAmongPredictedGenes() throws Exception {
        String[] fileNames = new String[]{
                "TCGABreastCancerScores_030415.txt",
                "TCGAColorectalCancerScores_030415.txt",
                "TCGAOvarianCancerScores_030415.txt"
        };
        String[] cancers = new String[]{"breast", "colorectal", "ovarian"};
        List<List<String>> list = new ArrayList<List<String>>();
        for (String fileName : fileNames) {
            fu.setInput(RESULT_DIR + fileName);
            String line = fu.readLine();
            final Map<String, Double> geneToScore = new HashMap<String, Double>();
            while ((line = fu.readLine()) != null) {
                String[] tokens = line.split("\t");
                geneToScore.put(tokens[0], new Double(tokens[1]));
            }
            fu.close();
            List<String> geneList = new ArrayList<String>(geneToScore.keySet());
            Collections.sort(geneList, new Comparator<String>() {
                public int compare(String gene1, String gene2) {
                    Double score1 = geneToScore.get(gene1);
                    Double score2 = geneToScore.get(gene2);
                    return score2.compareTo(score1);
                }
            });
            list.add(geneList);
        }
        int top = 30;
        System.out.println("Number of top genes: " + top);
        for (int i = 0; i < cancers.length - 1; i++) {
            List<String> list1 = list.get(i);
            System.out.println(cancers[i] + ": " + list1.size());
            for (int j = i + 1; j < cancers.length; j++) {
                List<String> list2 = list.get(j);
                System.out.println("\t" + cancers[j] + ": " + list2.size());
                Set<String> totalShared = InteractionUtilities.getShared(list1, list2);
                System.out.println("\tTotal genes shared: " + totalShared.size());
                // We want to analyze overlap of top 30 genes only
                // Need to make a copy to avoid modify the original list
                List<String> subList1 = new ArrayList<String>(list1.subList(0, top));
                // Need to focus to shared genes only
                subList1.retainAll(totalShared);
                System.out.println("\tTop 30 genes: " + subList1.size());
                List<String> subList2 = new ArrayList<String>(list2.subList(0, top));
                subList2.retainAll(totalShared);
                System.out.println("\tTop 30 genes: " + subList2.size());
                Set<String> shared = InteractionUtilities.getShared(subList1, subList2);
                System.out.println("\tShared in top 30 genes: " + shared.size());
                double pvalue = MathUtilities.calculateHypergeometricPValue(totalShared.size(),
                                                                            subList1.size(),
                                                                            subList2.size(),
                                                                            shared.size());
                System.out.println("\tp-value: " + pvalue);
            }
        }
    }
    
    @Test
    public void testKnownDriverGenes() throws Exception {
        FileUtility.initializeLogging();
        
        Set<String> allKnownGenes = getDriverGenes(null);
        System.out.println("All known cancer genes: " + allKnownGenes.size());
        System.out.println("Gene\tDriver");
        for (String gene : allKnownGenes)
            System.out.println(gene + "\ttrue");
        if (true)
            return;
        
        System.out.println("Known Breast Cancers:");
        Set<String> breast = getDriverGenes("breast");
        System.out.println("Known ovarian Cancers:");
        Set<String> ovarian = getDriverGenes("ovarian");
        System.out.println(ovarian);
        System.out.println("Known colorectal Cancers:");
        Set<String> colorectal = getDriverGenes("colorectal");
        // Perform overlap analysis
        List<Set<String>> list = new ArrayList<Set<String>>();
        list.add(breast);
        list.add(colorectal);
        list.add(ovarian);
        String[] names = new String[]{"breast", "colorectal", "ovarian"};
        for (int i = 0; i < names.length - 1; i++) {
            Set<String> set1 = list.get(i);
            System.out.println(names[i] + ": " + set1.size());
            for (int j = i + 1; j < names.length; j++) {
                Set<String> set2 = list.get(j);
                System.out.println("\t" + names[j] + ": " + set2.size());
                Set<String> shared = InteractionUtilities.getShared(set1, set2);
                System.out.println("\tShared: " + shared.size());
                double pvalue = MathUtilities.calculateHypergeometricPValue(R3Constants.TOTAL_HUMAN_GENES,
                                                                            set1.size(),
                                                                            set2.size(),
                                                                            shared.size());
                System.out.println("\tp-value: " + pvalue);
            }
        }
    }
    
    @Test
    public void runPermutationTest() throws Exception {
        FileUtility.initializeLogging();
        CancerDriverInstancesGenerator generator = getInstancesGenerator();
        
        FICancerDriverPredictor helper = new FICancerDriverPredictor();
        Set<String> driverGenes = helper.loadCancerCensusGenes();
        logger.info("Loaded driver genes: " + driverGenes.size());
        
        Set<String> fis = fu.loadInteractions(generator.FI_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
//        Map<String, Set<String>> geneToPartners = InteractionUtilities.generateProteinToPartners(fis);
//        System.out.println("TP53: " + geneToPartners.get("TP53").size());
        driverGenes.retainAll(fiGenes);
        logger.info("Driver genes in the FI network: " + driverGenes.size());
        
        FilteredClassifier fc = getClassifier();
        Classifier copy = AbstractClassifier.makeCopy(fc);
        // For output classifier information. In case the data be changed,
        // we want to use this original dataset first.
        Instances data = generator.generateNetworkOnlyInstances(driverGenes);
//        Instances data = generator.generateInstances(driverGenes);
        
        copy.buildClassifier(data);
        Evaluation eval = generator.new MutableDataEvalulation(data);
        eval.crossValidateModel(copy, data, 5, new Random(1));
        System.out.println("Data\tPrecision\tRecall\tAUC_ROC");
        System.out.println("RealData" + "\t" + 
                           eval.precision(0) + "\t" +
                           eval.recall(0) + "\t" +
                           eval.areaUnderROC(0));
        int permutation = 100;
        for (int i = 0; i < permutation; i++) {
//            Set<String> randomGenes = MathUtilities.randomSampling(fiGenes, driverGenes.size());
//            Instances randomData = generator.generateNetworkOnlyInstances(randomGenes);
//            Set<String> randomFIs = InteractionUtilities.generateRandomFIsViaSwitch(fis);
            Set<String> randomFIs = InteractionUtilities.generateRandomFIsViaConfigurationModel(fis);
            generator.setFIs(randomFIs);
//            geneToPartners = InteractionUtilities.generateProteinToPartners(randomFIs);
//            System.out.println("TP53: " + geneToPartners.get("TP53").size());
            Instances randomData = generator.generateNetworkOnlyInstances(driverGenes);
//            Instances randomData = generator.generateInstances(driverGenes);
            
            copy = AbstractClassifier.makeCopy(fc);
            copy.buildClassifier(randomData);
            eval = generator.new MutableDataEvalulation(randomData);
            eval.crossValidateModel(copy, randomData, 5, new Random(1));
            System.out.println("RandomData" + i + "\t" + 
                    eval.precision(0) + "\t" +
                    eval.recall(0) + "\t" +
                    eval.areaUnderROC(0));
        }
    }

    private FilteredClassifier getClassifier() throws Exception {
        Remove remove = new Remove(); // 1 is the first attribute
        remove.setAttributeIndices("1"); // The first index is gene names, which should be removed
        
        // The following options are copied from the GUI
        String[] options = Utils.splitOptions("-R 1.0E-8 -M -1");
        Logistic classifier = new Logistic();
        classifier.setOptions(options);
        
        // Use a FilteredClassifier
        FilteredClassifier fc = new FilteredClassifier();
        fc.setFilter(remove);
        fc.setClassifier(classifier);
        return fc;
    }

    @Test
    public void runClassifierOnMutableNetworkFeature() throws Exception {
        FileUtility.initializeLogging();
        CancerDriverInstancesGenerator generator = getInstancesGenerator();
        
        FICancerDriverPredictor helper = new FICancerDriverPredictor();
        Set<String> driverGenes = helper.loadCancerCensusGenes();
//        driverGenes.addAll(helper.loadNatureGADDriverGenes());
//        Set<String> driverGenes = getDriverGenes("breast");
//        Set<String> driverGenes = getDriverGenes("ovarian");
//        Set<String> driverGenes = getDriverGenes("colorectal");
//        Set<String> driverGenes = new FICancerDriverPredictor().loadNatureDriverGenes();
//        driverGenes = getDriverGenes(null);
        
        logger.info("Loaded driver genes: " + driverGenes.size());
        
        // For test
//        Set<String> fis = fu.loadInteractions(generator.FI_FILE_NAME);
//        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
//        driverGenes = MathUtilities.randomSampling(fiGenes, driverGenes.size());
        
        Instances data = generator.generateInstances(driverGenes);
        
        // For debugging 
        DataSink.write(RESULT_DIR + "OnTheFly.arff", data);
        // The first attribute returned from the above method call is a gene symbol,
        // we want to filter it on the fly
        Remove remove = new Remove(); // 1 is the first attribute
//        remove.setAttributeIndices("1,2,3,4,7"); 
//        remove.setAttributeIndices("1,5,6,8");
//        remove.setAttributeIndices("1,3,5"); // Remove network feature when only mutation is used in the data set.
//        remove.setAttributeIndices("1,6"); // Remove network feature
//        remove.setAttributeIndices("1,2,3,5"); // Network and mutation
//        remove.setAttributeIndices("1,2,3,5,6"); // Mutation only
//        remove.setAttributeIndices("1,2,3,4,5"); // Network feature only
//        remove.setAttributeIndices("1,3"); // Gene names and MutSigScore
//        remove.setAttributeIndices("1,2,3,4,7");
        remove.setAttributeIndices("1,2,3,6"); // The first index is gene names, which should be removed
        
        // The following options are copied from the GUI
        String[] options = Utils.splitOptions("-R 1.0E-8 -M -1");
        Logistic classifier = new Logistic();
        classifier.setOptions(options);
        
        // Use a FilteredClassifier
        FilteredClassifier fc = new FilteredClassifier();
        fc.setFilter(remove);
        fc.setClassifier(classifier);
        
        Classifier copy = AbstractClassifier.makeCopy(fc);
        // For output classifier information. In case the data be changed,
        // we want to use this original dataset first.
        copy.buildClassifier(data);
        
        Evaluation eval = generator.new MutableDataEvalulation(data);
//        Evaluation eval = new Evaluation(data);
        // Here seed 1 seems is used in the GUI of WEKA.
        eval.crossValidateModel(fc, data, 5, new Random(1));
        System.out.println("\n\n=== Classifier model (full training set) ===\n" + copy.toString());
        System.out.println(eval.toSummaryString("=== Stratified cross-validation ===\n=== Summary ===\n", true));
        System.out.println(eval.toClassDetailsString("=== Detailed Accuracy By Class ===\n"));
        System.out.println(eval.toMatrixString("=== Confusion Matrix ===\n"));
//        }
//        if (true)
//            return;
//        // We want to use another set of drivers genes to test
//        Set<String> testDriverGenes = helper.loadScienceReportDriverGenes(false);
//        logger.info("\nDriver genes from SR all: " + testDriverGenes.size());
//        testDriverGenes.removeAll(driverGenes);
//        logger.info("Remove drivers genes in the train data set: " + testDriverGenes.size());
//        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
//        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
//        testDriverGenes.retainAll(fiGenes);
//        fiGenes.removeAll(testDriverGenes);
//        fiGenes.removeAll(driverGenes);
//        logger.info("Test driver genes in the FI network: " + testDriverGenes.size());
//        System.out.println("\n\nResults in the test driver genes:");
//        // Check the score distribution for the testDriverGenes
//        for (int i = 0; i < data.numInstances(); i++) {
//            Instance instance = data.instance(i);
//            String gene = instance.stringValue(0);
//            double[] dist = copy.distributionForInstance(instance);
//            if (testDriverGenes.contains(gene)) {
//                System.out.println("true\t" + dist[0] + "\t" + dist[1]);
//            }
//            else if (fiGenes.contains(gene)) {
//                System.out.println("false\t" + dist[0] + "\t" + dist[1]);
//            }
//        }
    }

    private CancerDriverInstancesGenerator getInstancesGenerator() {
        // The following is used to check mutation and different network features
//        for (NetworkFeature networkFeature : NetworkFeature.values()) {
//            System.out.println("Mutation + Network Feautre (" + networkFeature + ")");
        // The following should be used as the default, which yields best result
        NetworkFeature networkFeature = NetworkFeature.ONE_HOP_SCORE_MAXIMUM;
//        networkFeature = networkFeature.TWO_HOP_SCORE_MINIMUM;
//        networkFeature = NetworkFeature.TWO_HOP_AVERAGE;
        boolean useRandomGenes = false;
        
        // Generating data
        CancerDriverInstancesGenerator generator = new CancerDriverInstancesGenerator();
//        generator.setAaChangeColumnInMAF("amino_acid_change_WU");
        generator.setAaChangeColumnInMAF("amino_acid_change");
//        generator.setFiColumnInMAF("ucsc_cons");
        // ICGC pan cancer
//        generator.setAaChangeColumnInMAF("Protein_Change");
//        generator.setFiColumnInMAF("MA_FI.score");
        generator.setFiColumnInMAF("MetaLR_rankscore");
//        generator.setMutSigResultFile("results/ICGC_PanCancer/mutsigcv_results_052115_sig_genes.txt");
//        generator.setFiPGMScoreFile("datasets/ICGC/2016_04/PGM_FI_Inference_Results_071416_31_19_Annot.txt");
        
        generator.setUseRandomGenesForNegative(useRandomGenes);
        generator.setNetworkFeature(networkFeature);
        // Just for comparison
//        generator.setNoMutationNormalization(true);
        return generator;
    }
    
    /**
     * Compare two sets of drivers genes from two TCGA pan-cancer data analyses.
     */
    @Test
    public void compareTwoSetsOfDriverGenes() throws Exception {
        List<Set<String>> geneLists = new ArrayList<Set<String>>();
        FICancerDriverPredictor helper = new FICancerDriverPredictor();
        Set<String> cancerCensusGenes = helper.loadCancerCensusGenes();
        System.out.println("Total cancer census genes: " + cancerCensusGenes.size());
        geneLists.add(cancerCensusGenes);
//        Set<String> natureGenes = helper.loadNatureDriverGenes();
        Set<String> natureGenes = helper.loadNatureGADDriverGenes();
        System.out.println("Total drivers from Nature: " + natureGenes.size());
        geneLists.add(natureGenes);
        Set<String> shared = InteractionUtilities.getShared(cancerCensusGenes, natureGenes);
        System.out.println("Shared with cancer census genes: " + shared.size());
        if (true)
            return;
        Set<String> srGenesHC = helper.loadScienceReportDriverGenes(true);
        System.out.println("Total HC drivers from SR: " + srGenesHC.size());
        geneLists.add(srGenesHC);
        for (int i = 0; i < geneLists.size(); i++) {
            Set<String> set1 = geneLists.get(i);
            for (int j = i + 1; j < geneLists.size(); j++) {
                Set<String> set2 = geneLists.get(j);
                shared = InteractionUtilities.getShared(set1, set2);
                double pvalue = MathUtilities.calculateHypergeometricPValue(R3Constants.TOTAL_HUMAN_GENES,
                                                                            set1.size(),
                                                                            set2.size(),
                                                                            shared.size());
                System.out.println(i + "\t" + j + "\t" + shared.size() + "\t" + pvalue);
            }
        }
        if (true)
            return;
        shared = InteractionUtilities.getShared(cancerCensusGenes, natureGenes);
        System.out.println("Shared with cancer census genes: " + shared.size());
        Set<String> srGenesAll = helper.loadScienceReportDriverGenes(false);
        System.out.println("Total drivers from SR: " + srGenesAll.size());
        shared = InteractionUtilities.getShared(cancerCensusGenes, natureGenes);
        System.out.println("Shared with cancer census genes: " + shared.size());
        shared = InteractionUtilities.getShared(natureGenes, srGenesHC);
        System.out.println("Shared between Nature and HC SR genes: " + shared.size());
        Set<String> natureCopy = new HashSet<String>(natureGenes);
        natureCopy.removeAll(shared);
        System.out.println("Not shared in Nature drivers: " + natureCopy);
        shared = InteractionUtilities.getShared(natureGenes, srGenesAll);
        System.out.println("Shared between Nature and All SR genes: " + shared.size());
        natureCopy = new HashSet<String>(natureGenes);
        natureCopy.removeAll(shared);
        System.out.println("Not shared in Nature drivers: " + natureCopy);
    }
    
    /**
     * This file is used to generate an ARFF file for WEKA
     * @throws IOException
     */
    @Test
    public void generateARFFFile() throws Exception {
        FICancerDriverPredictor helper = new FICancerDriverPredictor();
//        Set<String> driverGenes = helper.loadNatureDriverGenes();
        Set<String> driverGenes = helper.loadScienceReportDriverGenes(false);
        //      String arffFile = "results/DriverGenes/TCGABrca.arff";
        //      String arffFile = "results/DriverGenes/TCGABrca_120314.arff";
        // In this arff file, all genes, except drivers, are used as non-driver genes.
        // True for using the one-hop neighbor genes; false for using two-hop neighbor genes
        NetworkFeature networkFeature = NetworkFeature.TWO_HOP_AVERAGE;
        boolean useRandomGenes = false;
        
        // Nature driver with one hop neighbors and random genes
//        String arffFile = "TCGABrca_RandomGenes_One_Hop_120414.arff";
        // Nature driver with two hop neighbors and random genes
//        String arffFile = "TCGABrca_RandomGenes_Two_Hop_120414.arff";
        // Nature driver with one hop neighbors and all FI genes
//        String arffFile = "TCGABrca_AllGenes_One_Hop_120414.arff";
        // Nature driver with two hop neighbors and all FI genes
//        String arffFile = "TCGABrca_AllGenes_Two_Hop_120414.arff";
        
        // Science Reports driver with one hop neighbors and all FI genes
//        String arffFile = "TCGABrca_AllGenes_One_Hop_SC_HC_120414.arff";
        // Science Reports driver with two hop neighbors and all FI genes
//        String arffFile = "TCGABrca_AllGenes_Two_Hop_SC_HC_120414.arff";
        // Science Reports driver with two hop neighbors (all) and all FI genes
//        String arffFile = "TCGABrca_AllGenes_Two_Hop_SC_120414.arff";
        String arffFile = "TCGABrca_AllGenes_Two_Hop_SC_072115.arff";
        
        CancerDriverInstancesGenerator generator = new CancerDriverInstancesGenerator();
        generator.setUseRandomGenesForNegative(useRandomGenes);
        generator.setNetworkFeature(networkFeature);
        
        Instances instances = generator.generateInstances(driverGenes);
        for (int i = 0; i < instances.numInstances(); i++) {
            Instance instance = instances.get(i);
            System.out.println(instance.stringValue(0) + "\t" + instance.value(1));
        }
        
//        DataSink.write(RESULT_DIR + arffFile, instances);
        
    }
    
    /**
     * Check how many interactions among driver genes
     * @throws IOException
     */
    @Test
    public void checkInteractions() throws IOException {
//        Set<String> driverGenes = new FICancerDriverPredictor().loadNatureDriverGenes();
//        Set<String> driverGenes = new FICancerDriverPredictor().loadCancerCensusGenes();
        Set<String> driverGenes = new FICancerDriverPredictor().loadScienceReportDriverGenes(false);
        System.out.println("Total driver genes: " + driverGenes.size());
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        driverGenes.retainAll(fiGenes);
        System.out.println("In the FI network: " + driverGenes.size());
        System.out.println("FIs per gene: " + (double)fis.size() / fiGenes.size());
        if (true)
            return;
        
        Set<String> fisInDriver = InteractionUtilities.grepFIsContains(driverGenes, fis);
        System.out.println("Total FIs in driver genes: " + fisInDriver.size());
        
        // Permutation test
        List<Integer> fisNumbers = new ArrayList<Integer>();
        for (int i = 0; i < 1000; i++) {
            Set<String> randomGenes = MathUtilities.randomSampling(fiGenes, driverGenes.size());
            Set<String> randomFIs = InteractionUtilities.grepFIsContains(randomGenes, fis);
            fisNumbers.add(randomFIs.size());
        }
        Collections.sort(fisNumbers, new Comparator<Integer>() {
            public int compare(Integer n1, Integer n2) {
                return n2.compareTo(n1);
            }
        });
        for (Integer number : fisNumbers)
            System.out.println(number);
    }
    
    @Test
    public void checkNetworkFeatures() throws IOException {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Map<String, Set<String>> geneToPartners = InteractionUtilities.generateProteinToPartners(fis);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        CancerDriverInstancesGenerator helper = new CancerDriverInstancesGenerator();
        Set<String> driverGenes = new FICancerDriverPredictor().loadNatureDriverGenes();
        driverGenes.retainAll(fiGenes);
        // One hop network feature
        helper.setNetworkFeature(NetworkFeature.ONE_HOP_SCORE);
        Map<String, Double> geneToDriverRatio = helper.calculateGeneToDriverRatio(driverGenes, geneToPartners);
        SummaryStatistics stat = new SummaryStatistics();
        for (String gene : driverGenes) {
            double value = helper.calculateNetworkFeature(gene, 
                                                          driverGenes,
                                                          geneToPartners, 
                                                          geneToDriverRatio);
            stat.addValue(value);
        }
        System.out.println("Summary of one-hop network feature:");
        System.out.println("Mean: " + stat.getMean());
        System.out.println("SD: " + stat.getStandardDeviation());
        System.out.println("Min: " + stat.getMin());
        System.out.println("Max: " + stat.getMax());
        helper.setNetworkFeature(NetworkFeature.TWO_HOP_AVERAGE);;
        stat.clear();
        for (String gene : driverGenes) {
            double value = helper.calculateNetworkFeature(gene, driverGenes, geneToPartners, geneToDriverRatio);
            stat.addValue(value);
        }
        System.out.println("\n\nSummary of two-hop network feature:");
        System.out.println("Mean: " + stat.getMean());
        System.out.println("SD: " + stat.getStandardDeviation());
        System.out.println("Min: " + stat.getMin());
        System.out.println("Max: " + stat.getMax());
    }
}
