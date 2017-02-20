/*
 * Created on Jun 7, 2010
 *
 */
package org.reactome.cancer;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.r3.fi.SurvivalAnalysisHelper;
import org.reactome.r3.fi.SurvivalAnalysisResult;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * A set of utilities related to cancer analysis.
 * @author wgm
 *
 */
public class CancerAnalysisUtilitites {
    
    public static void printOutAnnotations(List<GeneSetAnnotation> annotations) {
        System.out.println("Pathway\tHitNumber\tPValue\tFDR\tGenes");
        for (GeneSetAnnotation annotation : annotations) {
            System.out.println(annotation.getTopic() + "\t" + 
                    annotation.getHitNumber() + "\t" + 
                    annotation.getPValue() + "\t" + 
                    annotation.getFdr() + "\t" + 
                    annotation.getHitIds());
        }
    }
    
    /**
     * Export a gene to value map as a Cytoscape attribute file.
     * @param geneToValue
     * @param fileName
     * @param attName
     * @throws IOException
     */
    public static <T> void outputMapAsCytoscapeAttribute(Map<String, T> geneToValue,
                                                         String fileName,
                                                         String attName,
                                                         Class<T> cls) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        // Export as attribute file
        fu.printLine(attName + " (" + cls.getName() + ")");
        for (String gene : geneToValue.keySet()) {
            fu.printLine(gene + "=" + geneToValue.get(gene));
        }
        fu.close();
    }
    
    /**
     * This helper method is used to generate groups of samples that are used for cross-validation.
     * @param sampleToGenes
     * @param fold
     * @return
     */
    public static List<Set<String>> randomCrossValidationSplit(Collection<String> samples,
                                                               int fold) {
        // List of genes used for cross-validations
        List<String> sampleList = new ArrayList<String>(samples);
        // Divide samples into 5 fold
        int size = sampleList.size() / fold;
        List<Set<String>> dividedSamples = new ArrayList<Set<String>>();
        RandomData randomizer = new RandomDataImpl();
        for (int i = 0; i < fold; i++) {
            Set<String> randomSample = MathUtilities.randomSampling(sampleList, size, randomizer);
            sampleList.removeAll(randomSample);
            dividedSamples.add(randomSample); // May have some left-over samples. Just ignore them!
        }
        return dividedSamples;
    }

    
    /**
     * Grep all samples into a set from a gene to sample to value map. The type of value
     * can be Double or Integer.
     * @param geneToSampleToValue
     * @return
     */
    public static <T extends Number> Set<String> grepSamplesInGeneToSampleToValue(Map<String, Map<String, T>> geneToSampleToValue) {
        Set<String> samples = new HashSet<String>();
        for (String gene : geneToSampleToValue.keySet()) {
            Map<String, ?> sampleToValue = geneToSampleToValue.get(gene);
            samples.addAll(sampleToValue.keySet());
        }
        return samples;
    }
    
    /**
     * Do a survival analysis for a list of gene modules based on on gene expression format file.
     * @param geneExpFileName
     * @param clinFileName
     * @param model
     * @param moduleIndex
     * @param clusters
     * @return
     * @throws Exception
     */
    public static SurvivalAnalysisResult doGeneExpClusterSurvivalAnalysis(String geneExpFileName,
                                                                          String clinFileName,
                                                                          String model,
                                                                          String moduleIndex,
                                                                          List<String> checkingSamples,
                                                                          List<Set<String>> clusters) throws Exception {
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = helper.loadGeneExpAfterAverage(geneExpFileName);
        String tmpFileName = R3Constants.TEMP_DIR + "TmpGeneToClusterScore.txt";
        helper.generateSampleToGeneExpClusters(geneToSampleToValue, 
                                               checkingSamples,
                                               tmpFileName, 
                                               clusters,
                                               null);
        File file = new File(tmpFileName);
        SurvivalAnalysisHelper survivalHelper = getSurvivalAnalysisHelper();
        SurvivalAnalysisResult output = survivalHelper.doSurvivalAnalysis(file, 
                                                                          new File(clinFileName),
                                                                          model,
                                                                          moduleIndex,
                                                                          "GeneExp",
                                                                          false);
        return output;
    }
    
    /**
     * Calculate FDRs for p-values generated from CoxPH analysis.
     * @param results
     * @return
     */
    public static String attachFDRsToSurvivalResults(String results) throws IOException {
        StringReader reader = new StringReader(results);
        BufferedReader br = new BufferedReader(reader);
        StringBuilder builder = new StringBuilder();
        String line = br.readLine();
        String[] tokens = line.split("\t");
        for (int i = 0; i < 3; i++) {
            builder.append(tokens[i]).append("\t");
        }
        builder.append("FDR");
        builder.append("\n");
        List<String> lines = new ArrayList<String>();
        List<Double> pvalues = new ArrayList<Double>();
        while ((line = br.readLine()) != null) {
            lines.add(line);
            tokens = line.split("\t");
            Double pvalue = new Double(tokens[2]);
            pvalues.add(pvalue);
        }
        br.close();
        Collections.sort(pvalues);
        List<Double> fdrs = MathUtilities.calculateFDRWithBenjaminiHochberg(pvalues);
        Map<Double, Double> pvalueToFDR = new HashMap<Double, Double>();
        for (int i = 0; i < pvalues.size(); i++)
            pvalueToFDR.put(pvalues.get(i), fdrs.get(i));
        for (String line1 : lines) {
            tokens = line1.split("\t");
            for (int i = 0; i < 3; i++)
                builder.append(tokens[i]).append("\t");
            builder.append(pvalueToFDR.get(new Double(tokens[2])));
            builder.append("\n");
        }
        return builder.toString();
    }
    
    /**
     * Initialize a SurvivalAnalysisHelper for doing survival analysis.
     * @return
     */
    public static SurvivalAnalysisHelper getSurvivalAnalysisHelper() {
        SurvivalAnalysisHelper survivalHelper = new SurvivalAnalysisHelper();
        survivalHelper.setrScript(R3Constants.survivalScript);
        survivalHelper.setTempDirName(R3Constants.TEMP_DIR);
        return survivalHelper;
    }
    
    /**
     * Filter sample to genes map based on sample number cutoff.
     * @param sampleToGenes
     * @param sampleNumber
     */
    public static void filterSampleToGenes(Map<String, Set<String>> sampleToGenes, 
                                           int sampleNumber) {
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        Set<String> totalGenes = new HashSet<String>();
        for (String gene : geneToSamples.keySet()) {
            Set<String> samples = geneToSamples.get(gene);
            if (samples.size() >= sampleNumber)
                totalGenes.add(gene);
        }
        for (String sample : sampleToGenes.keySet()) {
            Set<String> set = sampleToGenes.get(sample);
            set.retainAll(totalGenes);
        }
    }
    
    /**
     * This method is used to select genes in samples based on a sample cutoff value.
     * @param sampleNumber
     * @param sampleToGenes
     * @return
     */
    public static Set<String> selectGenesInSamples(int sampleNumber, Map<String, Set<String>> sampleToGenes) {
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        Set<String> rtn = new HashSet<String>();
        for (String gene : geneToSamples.keySet()) {
            Set<String> samples = geneToSamples.get(gene);
            if (samples.size() >= sampleNumber)
                rtn.add(gene);
        }
        return rtn;
    }
    
}
