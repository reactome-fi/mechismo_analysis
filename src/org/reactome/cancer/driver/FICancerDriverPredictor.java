/*
 * Created on Dec 2, 2014
 *
 */
package org.reactome.cancer.driver;

import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.apache.log4j.Logger;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;

import weka.classifiers.functions.Logistic;
import weka.classifiers.meta.FilteredClassifier;
import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Utils;
import weka.filters.unsupervised.attribute.Remove;

/**
 * This class implements a PGM_FI way for predicting cancer driver genes.
 * @author gwu
 *
 */
public class FICancerDriverPredictor {
    private final static Logger logger = Logger.getLogger(FICancerDriverPredictor.class);
    private final String DIR_NAME = "datasets/DriverGenes/";
    private final String NATURE_DRIVER_GENE_FILE_NAME = DIR_NAME + "Nature_PMID_24132290.txt";
    private final String NATURE_GAD_DRIVER_GENE_FILE_NAME = DIR_NAME + "Nature_PMID_24390350.txt";
    private final String SCIENCE_REPORT_DRIVER_GENE_FILE_NAME = DIR_NAME + "ScientificReports_PMID_24084849.csv";
    private final String RESULT_DIR_NAME = "results/DriverGenes/";
    
    /**
     * Default constructor.
     */
    public FICancerDriverPredictor() {
    }
    
    /**
     * Perform a running of classification and get a map from gene to predicted score.
     * @param driverGenes
     * @return
     * @throws Exception
     */
    public Instances calculateClassifierScore(Set<String> driverGenes,
                                              CancerDriverInstancesGenerator generator,
                                              Remove remove) throws Exception {
        logger.info("Loaded driver genes: " + driverGenes.size());
        Instances data = generator.generateInstances(driverGenes);
        
        // The first attribute returned from the above method call is a gene symbol,
        // we want to filter it on the fly
        if (remove == null) {
            remove = new Remove();
            remove.setAttributeIndices("1"); // 1 is the first attribute
        }
        
        // The following options are copied from the GUI
        String[] options = Utils.splitOptions("-R 1.0E-8 -M -1");
        Logistic classifier = new Logistic();
        classifier.setOptions(options);
        
        // Use a FilteredClassifier
        FilteredClassifier fc = new FilteredClassifier();
        fc.setFilter(remove);
        fc.setClassifier(classifier);
        fc.buildClassifier(data);
        
        // Add a new Attribute to store scores
        Attribute score = new Attribute("Score");
        data.insertAttributeAt(score, data.numAttributes());
        
        // Check the score distribution for the testDriverGenes
        Map<String, Double> geneToScore = new HashMap<String, Double>();
        for (int i = 0; i < data.numInstances(); i++) {
            Instance instance = data.instance(i);
            String gene = instance.stringValue(0);
            double[] dist = fc.distributionForInstance(instance);
            // The first value should be for the true class since we listed true before the false
            instance.setValue(data.numAttributes() - 1, dist[0]);
//            geneToScore.put(gene, dist[0]);
        }
        
        return data;
    }
    
    public List<Set<String>> getPredictedDriverGenesInBatch(Set<String> driverGenes,
                                                            double[] thresholds) throws Exception {
        Instances instances = calculateClassifierScore(driverGenes, null, null);
        
        // Genes to be returned
        List<Set<String>> geneSets = new ArrayList<Set<String>>();
        for (Double threshold : thresholds) {
            Set<String> geneSet = new HashSet<String>();
            for (Instance inst : instances) {
                String gene = inst.stringValue(0);
                Double score = inst.value(instances.numAttributes() - 1);
                if (score >= threshold)
                    geneSet.add(gene);
            }
            geneSets.add(geneSet);
        }
        return geneSets;
    }

    /**
     * Get a list of predicted driver genes using the default setting.
     * @param threshold
     * @return
     */
    public Set<String> getPredictedDriverGenes(Set<String> driverGenes,
                                               double threshold) throws Exception {
        List<Set<String>> geneSets = getPredictedDriverGenesInBatch(driverGenes, new double[]{threshold});
        return geneSets.get(0);
    }
    
    private double[] calculateConditionalProbs(Map<String, Map<String, Float>> sampleToGeneToValue,
                                               Set<String> driverGenes,
                                               Set<String> randomGenes) {
        double[] values = new double[4];
        // Random genes are treated as non-driver genes
        double[] probs = calculateConditionalProbs(sampleToGeneToValue, 
                                                   randomGenes);
        values[0] = probs[0];
        values[1] = probs[1];
        probs = calculateConditionalProbs(sampleToGeneToValue, driverGenes);
        values[2] = probs[0];
        values[3] = probs[1];
        return values;
    }

    private double[] calculateConditionalProbs(Map<String, Map<String, Float>> sampleToGeneToValue,
                                               Set<String> genes) {
        int count0 = 0; 
        int count1 = 0;
        for (String gene : genes) {
            for (String sample : sampleToGeneToValue.keySet()) {
                Map<String, Float> geneToValue = sampleToGeneToValue.get(sample);
                Integer value = geneToValue.get(gene).intValue();
                if (value == null)
                    continue;
                if (value == 0)
                    count0 ++;
                else
                    count1 ++;
            }
        }
        double[] probs = new double[2];
        probs[0] = (double) count0 / (count0 + count1);
        probs[1] = 1.0d - probs[0];
        return probs;
    }
    
    /**
     * Load the set of driver genes reported in Science Report.
     * @param highConfidenceOnly
     * @return
     * @throws IOException
     */
    public Set<String> loadScienceReportDriverGenes(boolean highConfidenceOnly) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(SCIENCE_REPORT_DRIVER_GENE_FILE_NAME);
        Set<String> genes = new HashSet<String>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(",");
            if (!highConfidenceOnly)
                genes.add(tokens[0]);
            else if (tokens[tokens.length - 1].equals("High Confidence Driver"))
                genes.add(tokens[0]);
        }
        fu.close();
        return genes;
    }
    
    /**
     * Load SR-HC genes that are not predicted based on functional intearctions, which are
     * genes are MutSig, or at least have two hits in CGC and four algorithms.
     * @return
     * @throws IOException
     */
    public Set<String> loadSRHCDriverNoFIGenes() throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(SCIENCE_REPORT_DRIVER_GENE_FILE_NAME);
        Map<String, Set<String>> geneToFeatures = new HashMap<String, Set<String>>();
        String line = fu.readLine();
        String[] headers = line.split(",");
        Set<String> rtn = new HashSet<String>();
        int count = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(",");
            if (!tokens[tokens.length - 1].equals("High Confidence Driver"))
                continue; // HC genes only
            if (tokens[tokens.length - 2].equals("Selected")) {
                rtn.add(tokens[0]);
                continue;
            }
            count = 0;
            for (int i = 1; i < 6; i++) {
                if (tokens[i].equals("CGC") || tokens[i].equals("Selected"))
                    count ++;
            }
            if (count > 1)
                rtn.add(tokens[0]);
        }
        fu.close();
        return rtn;
    }
    
    @Test
    public void checkSRDriverGenes() throws IOException {
        Set<String> sr = loadScienceReportDriverGenes(true);
        System.out.println("SR HR: " + sr.size());
        Set<String> srAll = loadScienceReportDriverGenes(false);
        System.out.println("SR All: " + srAll.size());
        Set<String> srNoFIs = loadSRHCDriverNoFIGenes();
        System.out.println("SR no FIs: " + srNoFIs.size());
        FileUtility fu = new FileUtility();
        fu.setInput(SCIENCE_REPORT_DRIVER_GENE_FILE_NAME);
        Set<String> mutSigGenes = new HashSet<String>();
        Set<String> cgcGenes = new HashSet<String>();
        Map<String, Set<String>> geneToFeatures = new HashMap<String, Set<String>>();
        String line = fu.readLine();
        String[] headers = line.split(",");
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(",");
            if (!tokens[tokens.length - 1].equals("High Confidence Driver"))
                continue; // For HC genes only
            if (tokens[tokens.length - 2].equals("Selected"))
                mutSigGenes.add(tokens[0]);
            if (tokens[1].equals("CGC"))
                cgcGenes.add(tokens[0]);
            for (int i = 2; i < 6; i++) {
                if (tokens[i].equals("Selected"))
                    InteractionUtilities.addElementToSet(geneToFeatures, tokens[0], headers[i]);
            }
        }
        fu.close();
        System.out.println("Total genes in geneToFetaures: " + geneToFeatures.size());
        int count = 0;
        for (String gene : geneToFeatures.keySet()) {
            Set<String> features = geneToFeatures.get(gene);
            if (features.size() > 1)
                count ++;
        }
        System.out.println("Have more than one feature: " + count);
        count = 0;
        for (String gene : geneToFeatures.keySet()) {
            Set<String> features = geneToFeatures.get(gene);
            if (features.size() == 1 && cgcGenes.contains(gene))
                count ++;
        }
        System.out.println("One feature and CGC genes: " + count);
        for (String gene : geneToFeatures.keySet()) {
            Set<String> features = geneToFeatures.get(gene);
            if (features.size() == 1 && !cgcGenes.contains(gene))
                count ++;
        }
        System.out.println("One feature and FIs: " + count);
    }
    
    /**
     * Load the set of driver genes reported in Nature.
     * @return
     * @throws IOException
     */
    public Set<String> loadNatureDriverGenes() throws IOException {
        Set<String> genes = new FileUtility().loadInteractions(NATURE_DRIVER_GENE_FILE_NAME);
        return genes;
    }
    
    @Test
    public void testLoadNatureDriverGenes() throws IOException {
        Set<String> genes = loadNatureDriverGenes();
        List<String> geneList = new ArrayList<String>(genes);
        Collections.sort(geneList);
        System.out.println("Total genes: " + geneList.size());
        for (String gene : geneList) 
            System.out.print(gene + ", ");
        System.out.println();
    }
    
    public Set<String> loadNatureDriverGenes(String cancerType) throws IOException {
        if (cancerType == null)
            return loadNatureDriverGenes();
        int index = NATURE_DRIVER_GENE_FILE_NAME.lastIndexOf(".");
        String fileName = NATURE_DRIVER_GENE_FILE_NAME.substring(0, index) + "_" + cancerType + ".txt";
        return new FileUtility().loadInteractions(fileName);
    }
    
    public Set<String> loadNatureGADDriverGenes() throws IOException {
        return new FileUtility().loadInteractions(NATURE_GAD_DRIVER_GENE_FILE_NAME);
    }
    
    public Set<String> loadNatureGADDriverGenes(String cancerType) throws IOException {
        int index = NATURE_GAD_DRIVER_GENE_FILE_NAME.lastIndexOf(".");
        String fileName = NATURE_GAD_DRIVER_GENE_FILE_NAME.substring(0, index) + "_" + cancerType + ".txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        fu.close();
        String[] tokens = line.split(", ");
        Set<String> genes = new HashSet<String>();
        for (String token : tokens) {
            index = token.indexOf("(");
            genes.add(token.substring(0, index).trim());
        }
        return genes;
    }
    
    public Set<String> loadCancerCensusGenes() throws IOException {
        return loadCancerCensusGenes(null);
    }
    
    public Set<String> loadCancerCensusGenes(String cancerType,
                                             Set<String> allowedMutationTypes) throws IOException {
        Set<String> genes = new HashSet<String>();
        String fileName = DIR_NAME + "cancer_gene_census.csv";
        Reader reader = new FileReader(fileName);
        Iterable<CSVRecord> records = CSVFormat.DEFAULT.parse(reader);
        boolean isHeader = true;
        int mutationTypeIndex = 12;
        int geneIndex = 0;
        int tumorTypeIndex = 7;
        for (CSVRecord record : records) {
            if (isHeader) {
                isHeader = false;
                continue;
            }
            if (cancerType != null) {
                String tumorType = record.get(tumorTypeIndex); // Somatic tumor
                if (!tumorType.contains(cancerType)) {
                    tumorType = record.get(tumorTypeIndex + 1); // Germline tumor
                    if (!tumorType.contains(cancerType)) {
                        continue;
                    }
                }
            }
            // Check mutation types
            String[] mutationTypes = record.get(mutationTypeIndex).split(",");
            boolean isAllowed = false;
            if (allowedMutationTypes == null)
                isAllowed = true;
            else {
                for (String type : mutationTypes) {
                    if (allowedMutationTypes.contains(type.trim())) {
                        isAllowed = true;
                        break;
                    }
                }
            }
            if (!isAllowed)
                continue;
//            String typeText = record.get(12);
            String gene = record.get(geneIndex);
            genes.add(gene);
        }
        return genes;
    }
    
    public Set<String> loadCancerCensusGenes(String cancerType) throws IOException {
        Set<String> allowedMutationTypes = getAllowedCGCMutationTypes();
        return loadCancerCensusGenes(cancerType, allowedMutationTypes);
    }
    
    @Test
    public void checkCancerCensusGenes() throws IOException {
        String fileName = DIR_NAME + "cancer_gene_census.csv";
        Reader reader = new FileReader(fileName);
        Iterable<CSVRecord> records = CSVFormat.DEFAULT.parse(reader);
        Set<String> types = new HashSet<String>();
        int mutationTypeIndex = 12;
        boolean isHeader = true;
        for (CSVRecord record : records) {
            if (isHeader) {
                isHeader = false;
                continue;
            }
            String typeText = record.get(mutationTypeIndex);
            System.out.println(typeText);
            String[] types1 = record.get(mutationTypeIndex).split(",");
            for (String type1 : types1) {
                type1 = type1.trim();
                if (type1.length() == 0)
                    continue;
                types.add(type1);
            }
        }
        System.out.println("\n\nTotal type: " + types.size());
        List<String> typeList = new ArrayList<String>(types);
        Collections.sort(typeList);
        for (String type : typeList)
            System.out.println(type);
    }
    
    private Set<String> getAllowedCGCMutationTypes() {
        String[] types = new String[] {
//                "A", // Amplifcation
//                "D", // Large deletion
//                "F", // Frameshift
//                "Gene Conversion", // Not sure
                // Missense mutations list
                "M", // Not sure. This is more like an error for Mis
                "Mis", // Missense
                "Mis. N", // More like an error: Mis, N
//                
//                "N", // nonsense
//                "O", // Other
//                "Promoter Mis",
//                "S", // Splice
//                "T" // translocation
        };
        Set<String> rtn = new HashSet<String>();
        for (String type : types)
            rtn.add(type);
        return rtn;
    }
    
    @Test
    public void testLoadCancerCensusGenes() throws IOException {
        Set<String> allCGCGenes = loadCancerCensusGenes(null, null);
        System.out.println("All CGC genes: " + allCGCGenes.size());
        Set<String> genes = loadCancerCensusGenes();
        System.out.println("Total cancer census genes: " + genes.size());
        Set<String> breastGenes = loadCancerCensusGenes("breast");
        System.out.println("Total breast cancer genes: " + breastGenes.size());
        for (String gene : breastGenes)
            System.out.println(gene);
    }
    
    /**
     * Combine multiple known driver genes together.
     * @throws IOException
     */
    @Test
    public void combineThreeGeneLists() throws IOException {
        Set<String> cgcGenes = loadCancerCensusGenes(null, null);
        System.out.println("CGC genes: " + cgcGenes.size());
        Set<String> natureGenes = loadNatureDriverGenes();
        System.out.println("Nature genes: " + natureGenes.size());
        Set<String> natureGadGenes = loadNatureGADDriverGenes();
        System.out.println("Nature GAD genes: " + natureGadGenes.size());
        Set<String> allGenes = new HashSet<String>(cgcGenes);
        allGenes.addAll(natureGadGenes);
        allGenes.addAll(natureGadGenes);
        System.out.println("Total genes: " + allGenes.size());
        List<String> allGenesList = new ArrayList<String>(allGenes);
        Collections.sort(allGenesList);
        String fileName = DIR_NAME + "CombinedDrivers_071615.txt";
        new FileUtility().saveCollection(allGenesList, fileName);
    }
    
}
