/*
 * Created on Jan 19, 2015
 *
 */
package org.reactome.cancer.driver;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.log4j.Logger;
import org.junit.Test;
import org.reactome.cancer.driver.CancerDriverInstancesGenerator.NetworkFeature;
import org.reactome.cancer.driver.CancerDriverInstancesGenerator.NetworkMutableInstances;
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
import weka.filters.unsupervised.attribute.Remove;

/**
 * This class implements Iterative Classification Algorithm (ICA) as a WEKA meta classifier. The wrapped
 * WEKA classifier is used as a local classifier. This class basically implements the ICA algorithm reported
 * by Sen et al in "Collective Classification in Network Data".
 * Note: This is an unfinished class. Especially, the cross-validation part cannot work correctly, generating
 * results different from using WEKA directly!
 * @author gwu
 *
 */
public class FIDriverICAPredictor  {
    private static Logger logger = Logger.getLogger(FIDriverICAPredictor.class);
    private int maxIteration = 100;
    // Used to generate Instances
    private CancerDriverInstancesGenerator instancesGenerator;
    // Wrapped Classifier from WEK
    private Classifier classifier;
    // Save information to be used later on
    private Map<String, double[]> geneToDist;
    
    /**
     * Default constructor.
     */
    public FIDriverICAPredictor() {
        instancesGenerator = new CancerDriverInstancesGenerator();
        instancesGenerator.setNetworkFeature(NetworkFeature.ONE_HOP_SCORE_MAXIMUM);
        instancesGenerator.setUseRandomGenesForNegative(false); // Assume all other genes are not driver
        instancesGenerator.setCacheData(true); // We don't want to load data multiple times
        
        // Some customization of instanceGenerator for the used maf file.
        instancesGenerator.setAaChangeColumnInMAF("amino_acid_change");
        instancesGenerator.setFiColumnInMAF("MetaLR_rankscore");
        
        // Default
        classifier = new Logistic();
        try {
            // The following options are copied from the GUI
            String[] options = Utils.splitOptions("-R 1.0E-8 -M -1");
            ((Logistic)classifier).setOptions(options);
        }
        catch(Exception e) {
            logger.error("FIDriverICAPredictor: " + e.getMessage(), e);
        }
        geneToDist = new HashMap<String, double[]>();
    }

    /**
     * Classify genes in the unlablledGenes set into either driver or non-driver by using
     * the driverGenes as the training set.
     * @param driverGenes
     * @param unlabeledGenes
     * @return the classifying results in gene to boolean (true for driver, false for non-driver).
     */
    public Map<String, Boolean> classify(Set<String> driverGenes,
                                         Set<String> unlabeledGenes,
                                         double threshold) throws Exception {
        // In the bootstrapping step, we assume unlabeled genes are non-driver
        Map<String, Boolean> geneToDriver = classifyOnce(driverGenes, unlabeledGenes, threshold);
        int iteration = 0;
        // Current set of driver genes
        Set<String> currentDrivers = new HashSet<String>();
        while (iteration < maxIteration) {
            logger.info("Iteration " + iteration);
            currentDrivers.clear();
            currentDrivers.addAll(driverGenes);
            for (String gene : geneToDriver.keySet()) {
                if (geneToDriver.get(gene))
                    currentDrivers.add(gene);
            }
            Map<String, Boolean> newGeneToDriver = classifyOnce(currentDrivers,
                                                                unlabeledGenes,
                                                                threshold);
            if (compareResults(geneToDriver, newGeneToDriver))
                break; // The results have stabilized
            geneToDriver = newGeneToDriver; // Keep the current value
            iteration ++;
        }
        return geneToDriver;
    }
    
    private boolean compareResults(Map<String, Boolean> oldGeneToDriver,
                                   Map<String, Boolean> newGeneToDriver) {
        for (String gene : oldGeneToDriver.keySet()) {
            Boolean oldValue = oldGeneToDriver.get(gene);
            Boolean newValue = newGeneToDriver.get(gene);
            if (!oldValue.equals(newValue))
                return false;
        }
        return true;
    }
    
    private Map<String, Boolean> classifyOnce(Set<String> driverGenes,
                                              Set<String> unlabeledGenes,
                                              double threshold) throws Exception {
        Instances instances = instancesGenerator.generateInstances(driverGenes);
        return classify(instances, unlabeledGenes, threshold);
    }
    
    private Map<String, Boolean> classify(Instances instances,
                                          Set<String> unlabeledGenes,
                                          double threshold) throws Exception {
        Classifier classifier = buildClassifier(instances);
        Map<String, Boolean> geneToLabel = new HashMap<String, Boolean>();
        // Check the score distribution for the testDriverGenes
        for (int i = 0; i < instances.numInstances(); i++) {
            Instance instance = instances.instance(i);
            String gene = instance.stringValue(0);
            if (unlabeledGenes.contains(gene)) {
                double[] dist = classifier.distributionForInstance(instance);           
                // The first value should be for the true class since we listed true before the false
                geneToLabel.put(gene, dist[0] >= threshold);
                // To be used later on
                geneToDist.put(gene, dist);
            }
        }
        return geneToLabel;
    }

    private Classifier buildClassifier(Instances instances) throws Exception {
        // The first attribute returned from the above method call is a gene symbol,
        // we want to filter it on the fly
        Remove remove = new Remove();
        remove.setAttributeIndices("1"); // 1 is the first attribute
        
        // Use a FilteredClassifier
        FilteredClassifier fc = new FilteredClassifier();
        fc.setFilter(remove);
        fc.setClassifier(classifier);
        // Make a copy so that the original classifier will not be used
        Classifier copy = AbstractClassifier.makeCopy(fc);
        // For output classifier information. In case the data be changed,
        // we want to use this original dataset first.
        copy.buildClassifier(instances);
        return copy;
    }
    
    /**
     * The test method.
     * @throws Exception
     */
    @Test
    public void testClassify() throws Exception {
        FileUtility.initializeLogging();
        Set<String> fis = new FileUtility().loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        FICancerDriverPredictor helper = new FICancerDriverPredictor();
        // Nature driver is used as the training data set
        Set<String> natureGenes = helper.loadNatureDriverGenes();
        natureGenes = helper.loadCancerCensusGenes();
        natureGenes.retainAll(fiGenes);
        // Science Reporter is used as the test data set
        Set<String> srGenes = helper.loadScienceReportDriverGenes(false);
        srGenes.retainAll(fiGenes);
        srGenes.removeAll(natureGenes);
        // Create a negative set
        fiGenes.removeAll(natureGenes);
        fiGenes.removeAll(srGenes);
        Set<String> negativeGenes = MathUtilities.randomSampling(fiGenes, srGenes.size());
        Set<String> unlabeled = new HashSet<String>();
        System.out.println("Negative Genes: " + negativeGenes.size());
        unlabeled.addAll(negativeGenes);
        System.out.println("Science Genes: " + srGenes.size());
        unlabeled.addAll(srGenes);
        System.out.println("Unlabeled: " + unlabeled.size());
        
        double threshold = 0.03d;
//        maxIteration = 0;
        Map<String, Boolean> geneToLabel = classify(natureGenes, unlabeled, threshold);
        int allTrue = 0;
        int realTrue = 0;
        for (String gene : geneToLabel.keySet()) {
            if (geneToLabel.get(gene)) {
                allTrue ++;
                if (srGenes.contains(gene))
                    realTrue ++;
            }
        }
        System.out.println("All True: " + allTrue);
        System.out.println("Real True: " + realTrue);
        System.out.println("Total: " + geneToLabel.size());
//        // For comparison
//        System.out.println("One step results:");
//        Set<String> predictedGenes = helper.getPredictedDriverGenes(natureGenes, threshold);
//        Set<String> shared = InteractionUtilities.getShared(predictedGenes, unlabeled);
//        System.out.println("All True: " + shared.size());
//        Set<String> predictedInSR = InteractionUtilities.getShared(predictedGenes, srGenes);
//        System.out.println("Real True: " + predictedInSR.size());
    }
    
    public void generateROCCurveFile(Instances data,
                                     Random random,
                                     int numFolds,
                                     String fileName) throws Exception {
        // The following is modified from Evaluation.crossValidate method:
        // Make a copy of the data we can reorder
        data = new Instances(data);
        data.randomize(random);
        if (data.classAttribute().isNominal()) {
          data.stratify(numFolds);
        }

        int totalPoints = 100;
        double step = 1.0 / totalPoints;
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        fu.printLine("Threshold\tFalse_Positive_Rate\tTrue_Positive_Rate");
        for (int k = 0; k < totalPoints + 1; k++) {
            double threshold = 0.0 + step * k;
            logger.info("Threshold: " + threshold);
            int totalNeg = 0;
            int totalPos = 0;
            int truePos = 0;
            int falsePos = 0;
            // Do the folds
            for (int i = 0; i < numFolds; i++) {
                logger.info("Crossvalidation " + i + "...");
                Instances test = data.testCV(numFolds, i);
                Set<String> unlabeled = new HashSet<String>();
                for (int j = 0; j < test.numInstances(); j++) {
                    Instance inst = test.get(j);
                    String gene = inst.stringValue(0);
                    unlabeled.add(gene);
                }
                Instances train = data.trainCV(numFolds, i, random);
                Set<String> drivers = new HashSet<String>();
                for (int j = 0; j < train.numInstances(); j++) {
                    Instance inst = train.get(j);
                    String cls = inst.stringValue(data.classAttribute());
                    if (cls.toLowerCase().equals("true")) {
                        String gene = inst.stringValue(0);
                        drivers.add(gene);
                    }
                }
                Map<String, Boolean> geneToLabel = classify(drivers, 
                                                            unlabeled, 
                                                            threshold); 
                for (int j = 0; j < test.numInstances(); j++) {
                    Instance inst = test.get(j);
                    String gene = inst.stringValue(0);
                    String label0 = inst.stringValue(data.classAttribute());
                    Boolean label = geneToLabel.get(gene);
                    if (label0.equals(Boolean.TRUE.toString())) {
                        totalPos ++;
                        if (label)
                            truePos ++;
                    }
                    else {
                        totalNeg ++;
                        if (label)
                            falsePos ++;
                    }
                }
            }
            fu.printLine(threshold + "\t" + 
                         falsePos / (double) totalNeg + "\t" + 
                         truePos / (double) totalPos);
        }
        fu.close();
    }
    
    /**
     * Run cross-validation based on WEKA Evaluation class.
     * @param evaluation
     */
    public void runCrossValidation(Instances data,
                                   Random random,
                                   int numFolds,
                                   double threshold) throws Exception {
        // Make a copy of the data we can reorder
        if (data instanceof NetworkMutableInstances) {
            Map<String, Set<String>> geneToPartners = ((NetworkMutableInstances)data).getGeneToPartners();
            // Make a new copy
            data = instancesGenerator.new NetworkMutableInstances(data);
            ((NetworkMutableInstances)data).setGeneToPartners(geneToPartners);
        }
        else
            data = new Instances(data);
        data.randomize(random);
        if (data.classAttribute().isNominal()) {
          data.stratify(numFolds);
        }

        Evaluation eval = new Evaluation(data);
        // Do the folds
        for (int i = 0; i < numFolds; i++) {
            logger.info("Crossvalidation " + i + "...");
            Instances test = data.testCV(numFolds, i);
            Set<String> unlabeled = new HashSet<String>();
            for (int j = 0; j < test.numInstances(); j++) {
                Instance inst = test.get(j);
                String gene = inst.stringValue(0);
                unlabeled.add(gene);
            }
            Instances train = data.trainCV(numFolds, i, random);
            Set<String> drivers = new HashSet<String>();
            for (int j = 0; j < train.numInstances(); j++) {
                Instance inst = train.get(j);
                String cls = inst.stringValue(data.classAttribute());
                if (cls.toLowerCase().equals("true")) {
                    String gene = inst.stringValue(0);
                    drivers.add(gene);
                }
            }
            
//            eval.setPriors(train);
//            Classifier classifier = buildClassifier(train);
//            eval.evaluateModel(classifier, test);
            
            classify(drivers, 
                     unlabeled, 
                     threshold); // The threshold value is used only for the method calling
            for (int j = 0; j < test.numInstances(); j++) {
                Instance inst = test.get(j);
                String gene = inst.stringValue(0);
                double[] dist = geneToDist.get(gene);
                eval.evaluationForSingleInstance(dist, inst, true);
                
//                eval.evaluateModelOnceAndRecordPrediction(classifier, inst);
            }
        }
        // Use the embed classifier for some information
        System.out.println("\n\n=== Classifier model (full training set) ===\n" + classifier.toString());
        System.out.println(eval.toSummaryString("=== Stratified cross-validation ===\n=== Summary ===\n", true));
        System.out.println(eval.toClassDetailsString("=== Detailed Accuracy By Class ===\n"));
        System.out.println(eval.toMatrixString("=== Confusion Matrix ===\n"));
    }
    
    /**
     * Check the performance using WEKA Evaluation class.
     * @throws Exception
     */
    @Test
    public void testCrossValidation() throws Exception {
        FileUtility.initializeLogging();
        maxIteration = 0;
        // Nature driver is used as the training data set
        FICancerDriverPredictor helper = new FICancerDriverPredictor();
//        Set<String> natureGenes = helper.loadNatureDriverGenes();
        Set<String> driverGenes = helper.loadCancerCensusGenes();
        Instances data = instancesGenerator.generateInstances(driverGenes);
        Random random = new Random(1);
        runCrossValidation(data, random, 5, 0.1d);
    }
    
    @Test
    public void testGenerateROCCurve() throws Exception {
        FileUtility.initializeLogging();
        // Nature driver is used as the training data set
        FICancerDriverPredictor helper = new FICancerDriverPredictor();
        Set<String> natureGenes = helper.loadNatureDriverGenes();
        Instances data = instancesGenerator.generateInstances(natureGenes);
        Random random = new Random(1);
        String fileName = CancerDriverAnalyzer.RESULT_DIR + "ROC_10_FOLD_TWO_HOP_MAX_012114.txt";
        generateROCCurveFile(data, random, 10, fileName);
    }
}
