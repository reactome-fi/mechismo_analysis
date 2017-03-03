/*
 * Created on Oct 4, 2012
 *
 */
package org.reactome.cancer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to process files downloaded from Broad Firehose.
 * @author gwu
 *
 */
public class FireHoseFileProcessor {
    private FileUtility fu = new FileUtility();
    
    public FireHoseFileProcessor() {
    }
    
    /**
     * Some of cancer stages have sub-types listed. Use this method to merge
     * them to a coarse analysis.
     * @throws IOException
     */
    @Test
    public void mergeDiseaseStages() throws IOException {
        String cancer = "BRCA";
        cancer = "COADREAD";
        
        String dirName = "datasets/TCGA/" + cancer + "/clinical/gdac.broadinstitute.org_" + cancer + ".Clinical_Pick_Tier1.Level_4.2014071500.0.0/";
        String inFileName = dirName + cancer + ".clin.merged.picked.transformed.txt";
        String outFileName = dirName + cancer + ".clin.merged.picked.transformed.stages.txt";
        
        fu.setInput(inFileName);
        fu.setOutput(outFileName);
        
        String line = fu.readLine();
        fu.printLine(line + "\tdiseasestage");
        
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String stage = tokens[6];
            if (stage.equals("NA")) {
                fu.printLine(line + "\tNA");
                continue;
            }
            System.out.println(stage);
            stage = stage.split(" ")[1];
            if (stage.matches("(\\w+)(a|b|c)")) {
                stage = stage.substring(0, stage.length() - 1);
            }
            fu.printLine(line + "\t" + stage);
        }
        fu.close();
    }
    
    @Test
    public void testLoadMAFFile() throws IOException {
//        String fileName = "test_data/tcga_ov/ov.maf.txt";
        String fileName = "test_data/tcga_brca/brca.maf.txt";
        fileName = "test_data/tcga_coadread/coadread.maf.txt";
        MAFFileLoader mafFileLoader = new MAFFileLoader();
        Map<String, Set<String>> sampleToGenes = mafFileLoader.loadSampleToGenes(fileName, false);
        System.out.println("Sample: " + sampleToGenes.size());
        Set<String> genes = new HashSet<String>();
        for (String sample : sampleToGenes.keySet())
            genes.addAll(sampleToGenes.get(sample));
        System.out.println("Genes: " + genes.size());
        String sample = "TCGA-13-1512";
        Set<String> mutatedGenes = sampleToGenes.get(sample);
        if (mutatedGenes == null)
            return;
        System.out.println(sample + ": " + mutatedGenes.size());
        List<String> list = new ArrayList<String>(mutatedGenes);
        Collections.sort(list);
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Map<String, Set<String>> idToParteners = InteractionUtilities.generateProteinToPartners(fis);
        for (String gene : list) {
            Set<String> partners = idToParteners.get(gene);
            if (partners == null)
                continue;
            System.out.println(gene + "\t" + partners.size());
        }
//        System.out.println();
//        // Find a sample having TP53 mutated
//        for (String sample1 : sampleToGenes.keySet()) {
//            genes = sampleToGenes.get(sample1);
//            if (genes.contains("TP53")) {
//                System.out.println(sample1 + ": " + genes.size());
//                list = new ArrayList<String>(genes);
//                Collections.sort(list);
//                for (String gene : list)
//                    System.out.println(gene);
//                break;
//            }
//        }
    }
    
    /**
     * Check if a line contains "NA" only
     * @param tokens
     * @return
     */
    private boolean isNALine(String[] tokens) {
        for (int i = 1; i < tokens.length; i++) {
            if (i % 4 == 1 && !tokens[i].equals("NA"))
                return false;
        }
        return true;
    }
    
    @Test
    public void zscoreTransformRNAExpression() throws Exception {
//        String dirName = TCGAOvarianCancerAnalyzer.OVARIAN_DIR_NAME + 
//                "FireHose/gdac.broadinstitute.org_OV.mRNA_Preprocess_Median.Level_4.2012082500.0.0/";
//        String inFileName = dirName + "OV.medianexp.transformed.txt";
//        String outFileName = dirName + "OV.medianexp.transformed.global.zscore.txt";
//        String dirName = R3Constants.GBM_DIR + 
//                "FireHose/gdac.broadinstitute.org_GBM.mRNA_Preprocess_Median.Level_4.2012102400.0.0/";
//        String dirName = R3Constants.GBM_DIR + 
//                "FireHose/20130222/gdac.broadinstitute.org_GBM-TP.mRNA_Preprocess_Median.Level_4.2013022200.0.0/";
//
//        String inFileName = dirName + "GBM-TP.medianexp.transformed.txt";
//        String outFileName = dirName + "GBM-TP.medianexp.transformed.global.zscore.txt";
        // Process miRNA Expression file
//        String dirName = "/Users/gwu/datasets/TCGA/OvarianCancer/FireHose/miRNA/gdac.broadinstitute.org_OV.miRseq_Mature_Preprocess.Level_3.2014051800.0.0/";
//        String inFileName = dirName + "OV.miRseq_mature_RPM_log2_processed.txt";
//        String outFileName = dirName + "OV.miRseq_mature_RPM_log2_processed_global_z.txt";
        
        // Process TCGA OV mRNA-seq file
//        String dirName = R3Constants.OVARIAN_DIR_NAME + "FireHose/mRNA/gdac.broadinstitute.org_OV.mRNAseq_Preprocess.Level_3.2014061400.0.0/";
//        // Fix one error related to SLC35E2: SLC35E2|728661 should be SLC35E2B|728661
        
        // BRCA miRNA
//        String dirName = "datasets/TCGA/BRCA/miRNA/";
//        String cancer = "BRCA";
        
        // COADREAD mRNASeq
//        String dirName = "datasets/TCGA/COADREAD/mRNA/gdac.broadinstitute.org_COADREAD.mRNAseq_Preprocess.Level_3.2014071500.0.0/";
//        String cancer = "COADREAD";
//        boolean isMiRNA = false;
//        // COADREAD miRNA
//        dirName = "datasets/TCGA/COADREAD/miRNA/gdac.broadinstitute.org_COADREAD.miRseq_Mature_Preprocess.Level_3.2014071500.0.0/";
//        isMiRNA = true;
        
        // HNSC
        String dirName = "datasets/TCGA/HNSC/gdac.broadinstitute.org_HNSC.mRNAseq_Preprocess.Level_3.2014101700.0.0/";
        String cancer = "HNSC";
        boolean isMiRNA = false;
        
        String inFileName = null;
        String outFileName = null;
        
        if (isMiRNA) {
            inFileName = dirName + cancer + ".miRseq_mature_RPM_log2_processed.txt";
            outFileName = dirName + cancer + ".miRseq_mature_RPM_log2_processed_global_z.txt";
        }
        else {
            inFileName = dirName + cancer + ".uncv2.mRNAseq_RSEM_normalized_log2.transformed.txt";
            outFileName = dirName + cancer + ".uncv2.mRNAseq_RSEM_normalized_log2_transformed_global_z.txt";
        }
        
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        helper.ztransformDataGlobally(inFileName, outFileName, false);
    }
    
    /**
     * Use this method to get patient names and remove annotations rows or columns
     * @throws IOException
     */
    @Test
    public void processmRNAExpressionFile() throws IOException {
        //String dirName = TCGAOvarianCancerAnalyzer.OVARIAN_DIR_NAME + 
        //      "FireHose/gdac.broadinstitute.org_OV.mRNA_Preprocess_Median.Level_4.2012082500.0.0/";
//        String inFileName = dirName + "OV.medianexp.txt";
//        String outFileName = dirName + "OV.medianexp.transformed.txt";
        String dirName = R3Constants.OVARIAN_DIR_NAME + "FireHose/052013/gdac.broadinstitute.org_OV-TP.mRNA_Preprocess_Median.Level_4.2013032600.0.0/";
        String inFileName = dirName + "OV-TP.medianexp.txt";
        String outFileName = dirName + "OV-TP.medianexp.transformed.txt";
        //String dirName = R3Constants.GBM_DIR + 
        //        "FireHose/gdac.broadinstitute.org_GBM.mRNA_Preprocess_Median.Level_4.2012102400.0.0/";
//        String dirName = R3Constants.GBM_DIR + 
//                "FireHose/20130222/gdac.broadinstitute.org_GBM-TP.mRNA_Preprocess_Median.Level_4.2013022200.0.0/";
//        String inFileName = dirName + "GBM-TP.medianexp.txt";
//        String outFileName = dirName + "GBM-TP.medianexp.transformed.txt";
        processFile(inFileName, outFileName);
    }
    
    /**
     * This method is used to process copy number variations data file downloaded
     * from the Firehose web site.
     * @throws IOException
     */
    @Test
    public void processCNVFile() throws IOException {
        //String dirName = TCGAOvarianCancerAnalyzer.OVARIAN_DIR_NAME + 
        //        "FireHose/gdac.broadinstitute.org_OV.CopyNumber_Gistic2.Level_4.2012082500.0.0/";
        //String dirName = R3Constants.GBM_DIR + 
        //        "FireHose/20130222/gdac.broadinstitute.org_GBM-TP.CopyNumber_Gistic2.Level_4.2013022200.0.0/";
        // Breast cancer
        //String dirName = "datasets/TCGA/BRCA/06231013/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2013052300.0.0/";
//        String dirName = R3Constants.OVARIAN_DIR_NAME + "FireHose/052013/gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2013032600.0.0/";
//        String dirName = R3Constants.OVARIAN_DIR_NAME + "FireHose/cnv/gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2014041600.0.0/";
//        String inFileName = dirName + "all_thresholded.by_genes.txt";
//        String outFileName = dirName + "all_thresholded.by_genes.transformed.txt";
        
        // Dir for COADREAD
//        String dirName = "datasets/TCGA/COADREAD/cnv/gdac.broadinstitute.org_COADREAD-TP.CopyNumber_Gistic2.Level_4.2014041600.0.0/";
        
        // For HNSC
        String dirName = "datasets/TCGA/HNSC/gdac.broadinstitute.org_HNSC-TP.CopyNumber_Gistic2.Level_4.2014071500.0.0/";
        
        // Correct file should be used for CNVs is all_date_by_genes.txt
        String inFileName = dirName + "all_data_by_genes.txt";
        String outFileName = dirName + "all_data_by_genes_transformed.txt";
        fu.setInput(inFileName);
        String line = fu.readLine();
        Set<String> allSamples = new HashSet<String>();
        Set<String> allPatients = new HashSet<String>();
        String[] tokens = line.split("\t");
        for (int i = 3; i < tokens.length; i++) {
            allSamples.add(tokens[i]);
            allPatients.add(tokens[i].substring(0, 12));
        }
        fu.close();
        System.out.println("Total samples: " + allSamples.size());
        System.out.println("Total patients: " + allPatients.size());
//        if (true)
//            return;
        // Check the focal event genes only
//        String inFileName = dirName + "focal_data_by_genes.txt";
//        String outFileName = dirName + "focal_data_by_genes_transformed.txt";
        processFile(inFileName, outFileName);
    }
    
    @Test
    public void processMethylationFile() throws IOException {
        String dirName = R3Constants.OVARIAN_DIR_NAME + "FireHose/052013/methylation_jhu_usc_.Level_3.2013042100.0.0/";
        String inFileName = dirName + "OV.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt";
        boolean isBigFile = false;
        String cancer = "OV";
        
        dirName = R3Constants.OVARIAN_DIR_NAME + "FireHose/methylation/gdac.broadinstitute.org_OV.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014061400.0.0/";
        inFileName = dirName + "OV.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt";
        
        // TCGA BRCA 450 files
        dirName = "/Users/gwu/datasets/TCGA/BRCA/methylation/gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014071500.0.0/";
        inFileName = dirName + "BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt";
        isBigFile = true;
        cancer = "BRCA";
        
        // TCGA BRCA 27 platform
        dirName = "/Users/gwu/datasets/TCGA/BRCA/methylation/gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014071500.0.0/";
        inFileName = dirName + "BRCA.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt";
        
        // TCGA COADREAD 450 files
        cancer = "COADREAD";
        dirName = "/Users/gwu/datasets/TCGA/COADREAD/methylation/gdac.broadinstitute.org_COADREAD.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014071500.0.0/";
        
//        inFileName = dirName + cancer + ".methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt";
//        isBigFile = true;
        
        // TCGA COADREAD 27 platform
        dirName = "datasets/TCGA/COADREAD/methylation/gdac.broadinstitute.org_COADREAD.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014071500.0.0/";
        inFileName = dirName + cancer + ".methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt";
        
        String outFileName = dirName + cancer + ".methylation.max.transformed.txt";
        processMethylationFile(inFileName, 
                               outFileName,
                               isBigFile);
    }
    
    /**
     * This method is used to merge two methylation data files that were generated from two microarray platforms
     * (450 and 27).
     * @throws IOException
     */
    @Test
    public void mergeTwoMethylationDataFiles() throws IOException {
        // TCGA BRCA
        String cancer = "BRCA";
        // TCGA COADREAD
        cancer = "COADREAD";
        
        String dirName = "datasets/TCGA/" + cancer + "/methylation/gdac.broadinstitute.org_" + cancer + ".Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014071500.0.0/";
        String fileName450 = dirName + cancer + ".methylation.max.transformed.txt";
        dirName = "datasets/TCGA/" + cancer + "/methylation/gdac.broadinstitute.org_" + cancer + ".Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014071500.0.0/";
        String fileName27 = dirName + cancer + ".methylation.max.transformed.txt";
        String output = "datasets/TCGA/" + cancer + "/methylation/" + cancer + ".methylation.max.transformed.merged.txt";
        
        Set<String> samplesIn450 = grepSamplesInHeader(fileName450);
        Set<String> samplesIn27 = grepSamplesInHeader(fileName27);
        System.out.println("Samples in 450: " + samplesIn450.size());
        System.out.println("Samples in 27: " + samplesIn27.size());
        Set<String> shared = InteractionUtilities.getShared(samplesIn450, samplesIn27);
        if (shared.size() > 0) {
            //throw new IllegalStateException("Two platforms have shared samples: " + shared.size() + " samples: " + shared);
            System.err.println("Two platforms have shared samples: " + shared.size() + " samples (" + shared + ")");
            // In case two platforms have shared samples, we will pick values from the 450 platform
            samplesIn27.removeAll(shared);
        }
        Set<String> genesIn450 = grepGenes(fileName450);
        Set<String> genesIn27 = grepGenes(fileName27);
        shared = InteractionUtilities.getShared(genesIn450, genesIn27);
        System.out.println("Genes in 450: " + genesIn450.size());
        System.out.println("Genes in 27: " + genesIn27.size());
        System.out.println("Shared: " + shared.size());
        // Since 450 have much more genes, we want to use 450 as a base
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = helper.loadGeneExp(fileName27, false);
        StringBuilder builder = new StringBuilder();
        fu.setOutput(output);
        fu.setInput(fileName450);
        String line = fu.readLine();
        // Create a combined header
        builder.append(line);
        for (String sample : samplesIn27)
            builder.append("\t").append(sample);
        fu.printLine(builder.toString());
        builder.setLength(0);
        while ((line = fu.readLine()) != null) {
            builder.append(line);
            int index = line.indexOf("\t");
            // Get gene
            String gene = line.substring(0, index);
            Map<String, Double> sampleToValue27 = geneToSampleToValue.get(gene);
            for (String sample : samplesIn27) {
                builder.append("\t");
                Double value = null;
                if (sampleToValue27 != null)
                    value = sampleToValue27.get(sample);
                if (value != null)
                    builder.append(value);
                else
                    builder.append("NA");
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    private Set<String> grepGenes(String fileName) throws IOException {
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> genes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            int index = line.indexOf("\t");
            genes.add(line.substring(0, index));
        }
        fu.close();
        return genes;
    }
    
    /**
     * Process a downloaded FireHose methylation so that it can be used by our
     * software infrastructure.
     * @param inFileName
     * @param outFileName
     * @throws IOException
     */
    private void processMethylationFile(String inFileName, 
                                        String outFileName,
                                        boolean isBigFile) throws IOException {
        fu.setInput(inFileName);
        fu.setOutput(outFileName);
        String line = fu.readLine();
        // The first line contains sample names
        String[] tokens = line.split("\t");
        StringBuilder builder = new StringBuilder();
        builder.append("Gene");
        int totalSample = 0;
        Set<String> samples = new HashSet<String>();
        Set<String> patients = new HashSet<String>();
        List<String> headers = new ArrayList<String>();
        for (int i = 1; i < tokens.length; i++) {
            headers.add(tokens[i].substring(0, 15)); // Containing sample types information.
            // Each sample has four columns
            if (i % 4 == 1 && isPrimarySample(headers, i - 1)) {
                String sample = tokens[i].substring(0, 12);
                builder.append("\t").append(sample);
                totalSample ++;
                samples.add(tokens[i]);
                patients.add(sample);
            }
        }
        fu.printLine(builder.toString());
        builder.setLength(0);
        //        System.out.println("Total samples: " + totalSample);
        System.out.println("Total samples: " + samples.size());
        System.out.println("Total patients: " + patients.size());
        if (samples.size() != patients.size())
            throw new IllegalStateException("Same patient has multiple samples!");
        // The second line is just some annotation, escape it
        line = fu.readLine();
        // Read the actual data
        // If it is a big file, there is not enough memory to record all mapped genes
        List<String> geneList = null;
        if (!isBigFile)
            geneList = new ArrayList<String>();
        int totalMappedGenes = 0;
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            // Get the gene annotation first
            String gene = tokens[2];
            if (gene.contains(";")) {
                // If a CpG id can be mapped to multiple genes, escape it
                continue;
            }
            if (gene.equals("NA"))
                continue; // Cannot find matched gene, escape it
            builder.append(gene);
            boolean hasValue = false;
            for (int i = 1; i < tokens.length; i++) {
                if (i % 4 == 1 && isPrimarySample(headers, i - 1)) {
                    // Don't output any row having only NA
                    if (!tokens[i].equals("NA")) {
                        hasValue = true;
                    }
                    // Get the beta_values, which should be used to indicate methylation levels
                    builder.append("\t").append(tokens[i]);
                }
            }
            if (!hasValue) { // If it doesn't has any value, escape it
                builder.setLength(0);
                continue;
            }
            fu.printLine(builder.toString());
            // IN case there is a big file, flush the output to keep memory usage low.
            fu.getPrintWriter().flush();
            builder.setLength(0);
            if (geneList != null)
                geneList.add(gene);
            totalMappedGenes ++;
            //            System.out.println("Size of gene list: " + geneList.size());
        }
        fu.close();
        System.out.println("Total genes: " + totalMappedGenes);
        List<String> duplicatedGenes = null;
        if (!isBigFile) {
            Set<String> geneSet = new HashSet<String>(geneList);
            boolean hasDuplicatedGenes = false;
            if (geneList.size() > geneSet.size()) {
                System.out.println("\tHas duplicated genes: " + geneSet.size());
                hasDuplicatedGenes = true;
            }
            else {
                System.out.println("\tNo gene duplicated!");
            }
            if (!hasDuplicatedGenes)
                return;
            Map<String, Integer> geneToCount = InteractionUtilities.countTermUsageInList(geneList);
            duplicatedGenes = new ArrayList<String>();
            for (String gene : geneToCount.keySet()) {
                Integer count = geneToCount.get(gene);
                if (count > 1) {
                    System.out.println(gene + "\t" + count);
                    duplicatedGenes.add(gene);
                }
            }
        }
        extractMaximumBetaValue(outFileName, 
                                duplicatedGenes);
    }


    private void extractMaximumBetaValue(String outFileName,
                                         List<String> duplicatedGenes) throws IOException {
        // as methylation values for those duplicated genes
        // Change file name
        File file = new File(outFileName);
        File destFile = new File(outFileName + ".tmp");
        file.renameTo(destFile);
        fu.setInput(destFile.getAbsolutePath());
        fu.setOutput(outFileName);
        String line = fu.readLine();
        List<String> sampleList = new ArrayList<String>();
        String[] tokens = line.split("\t");
        fu.printLine(line);
        Map<String, List<DescriptiveStatistics>> geneToStats = new HashMap<String, List<DescriptiveStatistics>>();
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            // If there is no list of duplicated gene is provided, treat all genes are duplicated.
            if (duplicatedGenes == null || duplicatedGenes.contains(tokens[0])) {
                List<DescriptiveStatistics> list = geneToStats.get(tokens[0]);
                if (list == null) {
                    list = new ArrayList<DescriptiveStatistics>();
                    for (int i = 1; i < tokens.length; i++) {
                        DescriptiveStatistics stat = new DescriptiveStatistics();
                        list.add(stat);
                    }
                    geneToStats.put(tokens[0], list);
                }
                for (int i = 1; i < tokens.length; i++) {
                    DescriptiveStatistics stat = list.get(i - 1);
                    if (!tokens[i].equals("NA"))
                        stat.addValue(new Double(tokens[i]));
                }
            }
            else {
                fu.printLine(line);
            }
        }
        StringBuilder builder = new StringBuilder();
        // Need to get the maximum value
        for (String gene : geneToStats.keySet()) {
            List<DescriptiveStatistics> list = geneToStats.get(gene);
            builder.setLength(0);
            builder.append(gene);
            for (DescriptiveStatistics stat : list) {
                if (stat.getN() == 0)
                    builder.append("\tNA");
                else
                    builder.append("\t").append(stat.getMax()); // Use the maximum value
                    //builder.append("\t").append(stat.getPercentile(50.0d)); // Median should be 50.0 percental
            }
            fu.printLine(builder.toString());
        }
        fu.close();
        System.out.println("Finished extraction.");
    }
    
    private boolean isPrimarySample(List<String> headers, int col) {
        String sampleId = headers.get(col);
        return sampleId.endsWith("-01");
    }
    
    @Test
    public void processmRNASeqFile() throws IOException {
        //String dirName = R3Constants.OVARIAN_DIR_NAME + "FireHose/052013/gdac.broadinstitute.org_OV.mRNAseq_Preprocess.Level_4.2013042100.0.0/";
//        String dirName = R3Constants.OVARIAN_DIR_NAME + "FireHose/mRNA/gdac.broadinstitute.org_OV.mRNAseq_Preprocess.Level_3.2014061400.0.0/";
        // Fix one error related to SLC35E2: SLC35E2|728661 should be SLC35E2B|728661
//        String inFileName = dirName + "OV.uncv2.mRNAseq_RSEM_normalized_log2_fix.txt";
//        String outFileName = dirName + "OV.uncv2.mRNAseq_RSEM_normalized_log2.transformed.txt";

        // For COADREAD
//        String dirName = "datasets/TCGA/COADREAD/mRNA/gdac.broadinstitute.org_COADREAD.mRNAseq_Preprocess.Level_3.2014071500.0.0/";
//        String cancer = "COADREAD";
        
        // For HNSC
//        String dirName = "datasets/TCGA/HNSC/gdac.broadinstitute.org_HNSC.mRNAseq_Preprocess.Level_3.2014101700.0.0/";
//        String cancer = "HNSC";
//        
//        boolean needValidatePrimary = false;
//        
//        String inFileName = dirName + cancer + ".uncv2.mRNAseq_RSEM_normalized_log2_fix.txt";
//        String outFileName = dirName + cancer + ".uncv2.mRNAseq_RSEM_normalized_log2.transformed.txt";
        
        // For BRCA December release
        String dirName = "datasets/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2014120600.0.0/";
        boolean needValidatePrimary = false;
        String inFileName = dirName + "BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt";
        String outFileName = dirName + "BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.transformed.txt";
        
        processmRNASeqFile(inFileName,
                           outFileName,
                           needValidatePrimary);
    }
    
    @Test
    public void checkMiRNAFile() throws IOException {
        // Ovarian
        String dirName = "/Users/gwu/datasets/TCGA/OvarianCancer/FireHose/miRNA/gdac.broadinstitute.org_OV.miRseq_Mature_Preprocess.Level_3.2014051800.0.0/";
        String cancer = "OV";
        // BRCA
        dirName = "datasets/TCGA/BRCA/miRNA/";
        cancer = "BRCA";
        // COADREAD
        dirName = "datasets/TCGA/COADREAD/miRNA/gdac.broadinstitute.org_COADREAD.miRseq_Mature_Preprocess.Level_3.2014071500.0.0/";
        cancer = "COADREAD";
        
        String fileName = dirName + cancer + ".miRseq_mature_RPM_log2.txt";
        String outFileName = dirName + cancer + ".miRseq_mature_RPM_log2_processed.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        System.out.println("Total samples: " + (tokens.length - 1));
        Set<String> sampleIdsWithTypes = new HashSet<String>();
        Set<String> sampleIds = new HashSet<String>();
        Set<String> primaryIds = new HashSet<String>();
        for (int i = 1; i < tokens.length; i++) {
            sampleIdsWithTypes.add(tokens[i].substring(0, 15));
            sampleIds.add(tokens[i].substring(0, 12));
            if (tokens[i].substring(0, 15).endsWith("-01"))
                primaryIds.add(tokens[i].substring(0, 12));
        }
        System.out.println("Total sample ids with types: " + sampleIdsWithTypes.size());
        System.out.println("Total sample ids: " + sampleIds.size());
        System.out.println("Total primary cancers: " + primaryIds.size());
        fu.close();
        
        fu.setInput(fileName);
        fu.setOutput(outFileName);
        
        StringBuilder builder = new StringBuilder();
        line = fu.readLine();
        tokens = line.split("\t");
        builder.append(tokens[0]);
        List<String> samples = new ArrayList<String>();
        for (int i = 1; i < tokens.length; i++) {
            String sample = tokens[i].substring(0, 15);
            samples.add(sample);
            if (sample.endsWith("-01"))
                builder.append("\t").append(sample.substring(0, 12));
        }
        fu.printLine(builder.toString());
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            if (shouldFilterOut(tokens, primaryIds.size()))
                continue;
            String miRNA = tokens[0].split("\\|")[0];
            builder.setLength(0);
            builder.append(miRNA);
            for (int i = 1; i < tokens.length; i++) {
                String sample = samples.get(i - 1);
                if (sample.endsWith("-01")) {
                    builder.append("\t").append(tokens[i]);
                }
            }
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    /**
     * Filter out rows that have less than 5% values for samples.
     * @param tokens
     * @param totalSample
     * @return
     */
    private boolean shouldFilterOut(String[] tokens, int totalSample) {
        int counter = 0;
        for (int i = 1; i < tokens.length; i++) {
            if (tokens[i].equals("NA"))
                counter ++;
        }
        return (counter > totalSample * 0.95d);
    }
    
    
    /**
     * The actual method that is used to process a mRNA-seq file downloaded from 
     * the FireHose web site. This method does the following this:
     * 1). Convert samples ids into patient ids, assuming all samples are cancer samples.
     * 2). Remove gene ids, and keep gene symbols
     * 3). Remove any rows that have NA only
     * 4). Check if any gene is duplicated.
     */
    private void processmRNASeqFile(String inFileName,
                                    String outFileName,
                                    boolean validPrimarySamples) throws IOException {
        checkSampleNames(inFileName);
        fu.setInput(inFileName);
        fu.setOutput(outFileName);
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        StringBuilder builder = new StringBuilder();
        builder.append("Gene");
        List<String> sampleIds = new ArrayList<String>();
        int sampleCount = 0;
        for (int i = 1; i < tokens.length; i++) {
            String sample = tokens[i];
            if (sample.length() > 15)
                sample = sample.substring(0, 15);
            sampleIds.add(sample);
            int index = sample.lastIndexOf("-");
            int code = new Integer(sample.substring(index + 1));
            if (code > 1 && validPrimarySamples)
                throw new IllegalStateException("Sample is not a primary cancer: " + sample);
            if (code > 1)
                continue;
            sample = sample.substring(0, 12);
            builder.append("\t").append(sample);
            sampleCount ++;
        }
        fu.printLine(builder.toString());
        List<String> geneList = new ArrayList<String>();
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            String gene = tokens[0];
            int index = gene.indexOf("|");
            gene = gene.substring(0, index);
            if (gene.equals("?"))
                continue; // No gene symbol has been mapped
            boolean hasValue = false;
            builder.setLength(0);
            builder.append(gene);
            for (int i = 1; i < tokens.length; i++) {
                String sample = sampleIds.get(i - 1);
                if (!sample.endsWith("-01"))
                    continue; // We need primary cancers only
                if (!tokens[i].equals("NA")) {
                    hasValue = true;
                }
                builder.append("\t").append(tokens[i]);
            }
            if (!hasValue)
                continue;
            fu.printLine(builder.toString());
            geneList.add(gene);
        }
        fu.close();
        System.out.println("Total samples in the output file: " + sampleCount);
        System.out.println("Total genes in list: " + geneList.size());
        Set<String> geneSet = new HashSet<String>(geneList);
        System.out.println("Total genes in set: "  + geneSet.size());
        if (geneList.size() == geneSet.size()) {
            return;
        }
        Map<String, Integer> geneToCounter = InteractionUtilities.countTermUsageInList(geneList);
        for (String gene : geneToCounter.keySet()) {
            Integer counter = geneToCounter.get(gene);
            if (counter > 1)
                System.out.println(gene + "\t" + counter);
        }
    }
    
    /**
     * Remove annotation columns
     * @param inFileName
     * @param outFileName
     * @throws IOException
     */
    private void processFile(String inFileName,
                             String outFileName) throws IOException {
        boolean isExpression = inFileName.contains("mRNA");
        FileUtility fu = new FileUtility();
        fu.setInput(inFileName);
        fu.setOutput(outFileName);
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        StringBuilder builder = new StringBuilder();
        builder.append(tokens[0]);
        int startIndex = isExpression ? 1 : 3;
        for (int i = startIndex; i < tokens.length; i++)
            builder.append("\t").append(tokens[i].substring(0, 12));
        fu.printLine(builder.toString());
        builder.setLength(0);
        int count = 0;
        while ((line = fu.readLine()) != null) {
            count ++;
            if (isExpression && count == 1)
                continue;
            tokens = line.split("\t");
            builder.append(tokens[0]);
            for (int i = startIndex; i < tokens.length; i++)
                builder.append("\t").append(tokens[i]);
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    @Test
    public void checkSampleNames() throws IOException {
        // This a CNV file for OV
//        String fileName = TCGAOvarianCancerAnalyzer.OVARIAN_DIR_NAME + 
//                "FireHose/gdac.broadinstitute.org_OV.CopyNumber_Gistic2.Level_4.2012082500.0.0/all_thresholded.by_genes.txt";
        // This is a mRNA expression file for OV
        //String fileName = TCGAOvarianCancerAnalyzer.OVARIAN_DIR_NAME + 
        //        "FireHose/gdac.broadinstitute.org_OV.mRNA_Preprocess_Median.Level_4.2012082500.0.0/OV.medianexp.transformed.txt";
//        String fileName = R3Constants.GBM_DIR + "Firehose/" + 
//                "gdac.broadinstitute.org_GBM.mRNA_Preprocess_Median.Level_4.2012102400.0.0/GBM.medianexp.txt"; 
        String fileName = "datasets/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2014120600.0.0/" +
                          "BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.transformed.txt";
        
        checkSampleNames(fileName);
    }

    private void checkSampleNames(String fileName) throws IOException {
        boolean isExpression = fileName.contains("mRNA");
        System.out.println(fileName);
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        System.out.println("Total tokens: " + tokens.length);
        List<String> samples = new ArrayList<String>();
        Set<String> patients = new HashSet<String>();
        //int start = isExpression ? 1 : 3;
        int start = 1;
        for (int i = start; i < tokens.length; i++) {
            String sample = tokens[i];
            String patient = sample.substring(0, 12);
            samples.add(sample);
            patients.add(patient);
        }
        fu.close();
        System.out.println("Total samples: " + samples.size());
        System.out.println("Total patients: " + patients.size());
//        if (samples.size() != patients.size())
//            throw new IllegalStateException("Some patients have duplicated samples in file: " + fileName);
    }
    
    @Test
    public void checkSamplesInMethylaiton() throws IOException {
        // TCGA OV data set
        String fileName450 = R3Constants.OVARIAN_DIR_NAME + "Firehose/methylation/" + 
                "gdac.broadinstitute.org_OV.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014071500.0.0/" +
                "OV.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt";
        String fileName27 = R3Constants.OVARIAN_DIR_NAME + "Firehose/methylation/" + 
                "gdac.broadinstitute.org_OV.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014061400.0.0/" + 
                    "OV.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt";
        // TCGA breast cancer
        fileName450 = "datasets/TCGA/BRCA/methylation/" + 
                "gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014071500.0.0/" +
                "BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt";
        fileName27 = "datasets/TCGA/BRCA/methylation/" + 
                "gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014071500.0.0/" + 
                "BRCA.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt";
        // COADREAD
        fileName450 = "datasets/TCGA/COADREAD/methylation/" + 
                "gdac.broadinstitute.org_COADREAD.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014071500.0.0/" + 
                "COADREAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt";
        fileName27 = "datasets/TCGA/COADREAD/methylation/" + 
                "gdac.broadinstitute.org_COADREAD.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014071500.0.0/" +
                "COADREAD.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt";
        
        Set<String> samplesIn450 = grepSamplesInHeader(fileName450);
        System.out.println("Samples in 450: " + samplesIn450.size());
        for (String sample : samplesIn450)
            System.out.println(sample);
        Set<String> samplesIn27 = grepSamplesInHeader(fileName27);
        System.out.println("Samples in 27: " + samplesIn27.size());
        for (String sample : samplesIn27) {
            System.out.println(sample);
        }
        Set<String> shared = InteractionUtilities.getShared(samplesIn450, samplesIn27);
        System.out.println("Shared: " + shared.size());
    }
    
    private Set<String> grepSamplesInHeader(String fileName) throws IOException {
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> samples = new HashSet<String>();
        String[] tokens = line.split("\t");
        for (int i = 1; i < tokens.length; i++) {
            samples.add(tokens[i]);
        }
        fu.close();
        return samples;
    }
    
}
