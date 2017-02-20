/*
 * Created on Jun 22, 2010
 *
 */
package org.reactome.annotate;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.random.JDKRandomGenerator;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;

/**
 * This class is used to annotate gene sets based on pathways.
 * @author wgm
 *
 */
public class PathwayBasedAnnotator {
    
    // Use to control output
    protected Double pvalueThreshold = null;
    protected Double fdrThreshold = null;
    // This list of genes are used for permutation test
    private Collection<String> randomGenes;
    // For some helping function
    private AnnotationHelper annotationHelper;
    // Two ways to calculate FDRs: default is permutation
    private boolean useBenjaminiHochbergForFDR;
    
    public PathwayBasedAnnotator() {
        annotationHelper = new AnnotationHelper();
    }
    
    public boolean isUseBenjaminiHochbergForFDR() {
        return useBenjaminiHochbergForFDR;
    }

    public void setUseBenjaminiHochbergForFDR(boolean useBenjaminiHochbergForFDR) {
        this.useBenjaminiHochbergForFDR = useBenjaminiHochbergForFDR;
    }

    public void setAnnotationHelper(AnnotationHelper helper) {
        this.annotationHelper = helper;
    }
    
    public AnnotationHelper getAnnotationHelper() {
        return this.annotationHelper;
    }
    
    public void setPValueThreshold(double threshold) {
        this.pvalueThreshold = threshold;
    }
    
    public void setFDRThreshold(double threshold) {
        this.fdrThreshold = threshold;
    }
    
    public void setRandomGenes(Collection<String> genes) {
        this.randomGenes = genes;
    }
    
    /**
     * Annotate a set of genes.
     * @param genes
     * @param type one of types: Pathways, BP, MF and CC in GO.
     * @return
     * @throws Exception
     */
    public List<GeneSetAnnotation> annotateGenesWithFDR(Collection<String> genes,
                                                        AnnotationType type) throws Exception {
        Map<String, Set<String>> geneToPathways = annotationHelper.loadProteinNameToTermsMap(type);
        return annotateGeneSet(genes,
                               geneToPathways);
    }
    
    /**
     * Annotate a set of genes with all Reactome pathways, which are organized in a hierarchical way.
     * @param genes
     * @return
     * @throws Exception
     */
    public List<GeneSetAnnotation> annotateGenesWithReactomePathways(Collection<String> genes) throws Exception {
        Map<String, Set<String>> geneToPathways = annotationHelper.loadProteinNameToReactomePathwaysMap();
        return annotateGeneSet(genes, geneToPathways);
    }
    
    /**
     * Annotate a set of ids of reactions extracted from Reactome based on pathways.
     * @param reactionIds
     * @return
     * @throws Exception
     */
    public List<GeneSetAnnotation> annotateReactionsWithReactomePathways(Collection<String> reactionIds) throws Exception {
        Map<String, Set<String>> reactionIdToPathways = annotationHelper.loadReactomeReactionToPathwaysMap();
        return annotateGeneSet(reactionIds, reactionIdToPathways);
    }
    
    private List<GeneSetAnnotation> outputAnnotationWithFDR(Map<String, Integer> pathwayToGeneNumber,
                                                            Map<String, Double> pathwayToRatio,
                                                            List<GeneSetAnnotation> realTopicInfoes,
                                                            List<GeneSetAnnotation> allTopicInfoes,
                                                            int permutationNumber) {
        annotationHelper.sortGeneSetAnnotation(realTopicInfoes);
        annotationHelper.sortGeneSetAnnotation(allTopicInfoes);
        int realTotal = pathwayToRatio.size();
        int randomTotal = pathwayToRatio.size() * permutationNumber;
//        String header = "Pathway\tRatioOfProteinInPathway\tNumberOfProteinInPathway\tProteinFromSample\tP-Value\tFDR\tIDs";
//        System.out.println(header);
        int index = 0;
        List<GeneSetAnnotation> rtn = new ArrayList<GeneSetAnnotation>();
        for (GeneSetAnnotation info : realTopicInfoes) {
            // Need to get count for p-value less than topicInfo.pvalue
            int randomLessThan = countAnnotationsLessThan(allTopicInfoes, 
                                                          info.getPValue());
            int realLessThan = countAnnotationsLessThan(realTopicInfoes, 
                                                        info.getPValue());
            double fdr = calculateFDR(randomLessThan,
                                      realLessThan,
                                      randomTotal,
                                      realTotal);
            if (fdrThreshold != null && fdr > fdrThreshold)
                continue;
            if (pvalueThreshold != null && info.getPValue() > pvalueThreshold)
                continue;
            String fdrText = fdrToString(fdr, randomLessThan, realLessThan, randomTotal, realTotal);
//            String line = info.getTopic() + "\t" + 
//            String.format("%.3e", pathwayToRatio.get(info.getTopic())) + "\t" + 
//            pathwayToGeneNumber.get(info.getTopic()) + "\t" +
//            info.getHitNumber() + "\t" + 
//            String.format("%.3e", info.getPValue()) + "\t" + 
//            fdrText + "\t" + 
//            info.getHitIds();
//            System.out.println(line);
            index ++;
            // Clone a new GeneSetAnnotation to avoid any possibe error.
            GeneSetAnnotation result = new GeneSetAnnotation();
            result.setTopic(info.getTopic());
            result.setRatioOfTopic(pathwayToRatio.get(info.getTopic()));
            result.setNumberInTopic(pathwayToGeneNumber.get(info.getTopic()));
            result.setHitNumber(info.getHitNumber());
            result.setPValue(info.getPValue());
            result.setFdr(fdrText);
            result.setHitIds(new ArrayList<String>(info.getHitIds()));
            rtn.add(result);
        }
        return rtn;
    }
    
    private String fdrToString(double fdr,
                               int randomLessThan,
                               int realLessThan,
                               int randomTotal,
                               int realTotal) {
        String format = "%." + (int) Math.log10(realTotal) + "e";
        if (fdr != 0.0)
            return String.format(format, fdr);
        // Get the smallest number
        double tmp = calculateFDR(1,
                                  realLessThan,
                                  randomTotal,
                                  realTotal);
        return "<" + String.format(format, tmp);                              
    }
    
    private List<GeneSetAnnotation> _annotateGeneSet(Collection<String> ids,
                                                     Map<String, Set<String>> idToTopics,
                                                     Map<String, Double> topicToRatio,
                                                     Map<String, Integer> topicToNumber) {
        // To track how many ids can be annotated
        int touchedId = 0;
        List<String> topics = new ArrayList<String>();
        Map<String, List<String>> touchedTopicToIds = new HashMap<String, List<String>>();
        for (String id : ids) {
            Set<String> tmp = idToTopics.get(id);
            if (tmp == null || tmp.size() == 0)
                continue;
            touchedId ++;
            for (String topic : tmp) {
                // Use List, instead of Set, to keep redundant gene list,
                // from multiple samples
                List<String> tIds = touchedTopicToIds.get(topic);
                if (tIds == null) {
                    tIds = new ArrayList<String>();
                    touchedTopicToIds.put(topic, tIds);
                }
                tIds.add(id);
            }
        }
        List<GeneSetAnnotation> annotations = new ArrayList<GeneSetAnnotation>();
        for (Iterator<String> it = touchedTopicToIds.keySet().iterator(); it.hasNext();) {
            String topic = it.next();
            List<String> tIds = touchedTopicToIds.get(topic);
            int number = tIds.size();
            // No correction with the following single statement
            // This p-value is based on binomial test
            double pValue = MathUtilities.calculateBinomialPValue(topicToRatio.get(topic),
                                                          touchedId,
                                                          number);
            //                pValue *= MathUtilities.calculateDbinom(touchedId, 
            //                                                        ids.size(), 
            //                                                        coverage);
            //                double correction = MathUtilities.calculatePValue(R3Constants.COVERAGE, 
            //                                                                  ids.size(), 
            //                                                                  touchedId);
            //                pValue *= correction;
            //            double pValue = MathUtilities.calculateEnrichment(topicToRatio.get(topic),
            //                                                              ids.size(),
            //                                                              number);
            GeneSetAnnotation annotation = new GeneSetAnnotation();
            annotation.setHitNumber(number);
            annotation.setPValue(pValue);
            annotation.setTopic(topic);
            annotation.setHitIds(tIds);
            annotation.setRatioOfTopic(topicToRatio.get(topic));
            annotation.setNumberInTopic(topicToNumber.get(topic));
            annotations.add(annotation);
        }
        annotationHelper.sortGeneSetAnnotation(annotations);
        return annotations;
    }
    
    /**
     * Calculate Benjamini-Hochberg FDRs and attach to GeneSetAnnotation objects.
     * @param annotations the list should be sorted based on p-value already.
     */
    private void attachBHFDRs(List<GeneSetAnnotation> annotations) {
        List<Double> pvalues = new ArrayList<Double>();
        for (GeneSetAnnotation annotation : annotations)
            pvalues.add(annotation.getPValue());
        List<Double> fdrs = MathUtilities.calculateFDRWithBenjaminiHochberg(pvalues);
        for (int i = 0; i < annotations.size(); i++) {
            GeneSetAnnotation annotation = annotations.get(i);
            annotation.setFdr(fdrs.get(i) + "");
        }
    }
    
    /**
     * The actual method for doing enrichment analysis.
     * @param genes
     * @param geneToPathways
     * @return
     * @throws IOException
     */
    private List<GeneSetAnnotation> annotateGeneSet(Collection<String> genes,
                                                    Map<String, Set<String>> geneToPathways) throws IOException {
        Map<String, Integer> pathwayToGeneNumber = annotationHelper.countProteinsInTopics(geneToPathways);
        //System.out.println("Total pathways: " + pathwayToGeneNumber.size());
        Map<String, Double> pathwayToRatio = annotationHelper.calculateTopicToRatio(geneToPathways.size(),
                                                                                    pathwayToGeneNumber);
        List<GeneSetAnnotation> realTopicInfoes = _annotateGeneSet(genes,
                                                                   geneToPathways, 
                                                                   pathwayToRatio,
                                                                   pathwayToGeneNumber);
        if (useBenjaminiHochbergForFDR) {
            attachBHFDRs(realTopicInfoes);
            return realTopicInfoes;
        }
        RandomData randomizer = new RandomDataImpl(new JDKRandomGenerator());
        Set<String> sampleGenes = null;
        if (randomGenes != null && randomGenes.size() > 0) {
            sampleGenes = new HashSet<String>(randomGenes);
        }
        else {
            sampleGenes = annotationHelper.loadRandomGenes();
        }
        List<GeneSetAnnotation> allTopicInfoes = new ArrayList<GeneSetAnnotation>();
        int permutationNumber = 1000;
        for (int i = 0; i < permutationNumber; i++) {
            Object[] objs = randomizer.nextSample(sampleGenes, 
                                                  genes.size());
            List<String> sample = InteractionUtilities.convertArrayToList(objs);
            List<GeneSetAnnotation> topicInfoes = _annotateGeneSet(sample, 
                                                                   geneToPathways, 
                                                                   pathwayToRatio,
                                                                   pathwayToGeneNumber);
            allTopicInfoes.addAll(topicInfoes);
        }
        return outputAnnotationWithFDR(pathwayToGeneNumber, 
                                       pathwayToRatio,
                                       realTopicInfoes, 
                                       allTopicInfoes,
                                       permutationNumber);
    }
    
    private double calculateFDR(int randomLessThan,
                                int realLessThan,
                                int randomTotal,
                                int realTotal) {
        double denom = (double) realLessThan / realTotal;
        double num = (double) randomLessThan / randomTotal;
        double fdr = num / denom;
        if (fdr > 1.0d)
            fdr = 1.0d;
        return fdr;
    }
    
    private int countAnnotationsLessThan(List<GeneSetAnnotation> all, double cutoff) {
        // All should be sorted already
        for (int i = 0; i < all.size(); i++) {
            GeneSetAnnotation info = all.get(i);
            if (info.getPValue() > cutoff) {
                return i;
            }
        }
        return all.size();
    }
    
}
