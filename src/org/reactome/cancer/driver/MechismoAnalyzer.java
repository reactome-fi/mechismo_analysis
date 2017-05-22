/*
 * Created on Apr 10, 2017
 *
 */
package org.reactome.cancer.driver;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.junit.Test;
import org.reactome.annotate.AnnotationHelper;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.px.util.InteractionUtilities;
import org.reactome.r3.Interactome3dAnalyzer;
import org.reactome.r3.ReactomeAnalyzer;
import org.reactome.r3.util.FileUtility;

/**
 * @author gwu
 *
 */
public class MechismoAnalyzer {
    private String dirName = "datasets/Mechismo/";
    private String outputFileName = dirName + "COSMICv74_somatic_noSNPs_GWS_mechismo_output.tsv";
    private FileUtility fu = new FileUtility();
    
    /**
     * Default constructor.
     */
    public MechismoAnalyzer() {
    }
    
    @Test
    public void performPathwayEnrichmentAnalysisForReactions() throws Exception {
        String reactionsFileName = "results/ReactionsInMechismo_051017.txt";
        
        AnnotationHelper annotationHelper = new AnnotationHelper();
        annotationHelper.setReactionIdToPathwayFile("resources/ReactomeReactionsToPathways_051017.txt");
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        annotator.setUseBenjaminiHochbergForFDR(true);
        annotator.setAnnotationHelper(annotationHelper);
        
        // Choose a threshold to get reactions to be analyzed
        double[] thresholds = new double[]{5.0d, 4.0d, 3.0d};
        for (double threshold : thresholds) {
            List<String> reactionIds = getReactionIds(reactionsFileName, threshold);
            List<GeneSetAnnotation> annotations = annotator.annotateReactionsWithReactomePathways(reactionIds);
            System.out.println("Threshold: " + threshold + " with total selected reactions: " + reactionIds.size());
            System.out.println("Pathway\tNumberInPathway\tRatioOfPathway\tHitNumber\tpValue\tFDR\tHitIds");
            for (GeneSetAnnotation annotation : annotations) {
                System.out.println(annotation.getTopic() + "\t" + 
                                   annotation.getNumberInTopic() + "\t" + 
                                   annotation.getRatioOfTopic() + "\t" + 
                                   annotation.getHitNumber() + "\t" + 
                                   annotation.getPValue() + "\t" + 
                                   annotation.getFdr() + "\t" + 
                                   annotation.getHitIds());
            }
            System.out.println();
        }
    }
    
    private List<String> getReactionIds(String fileName, 
                                        double threshold) throws IOException {
        fu.setInput(fileName);
        String line = fu.readLine();
        List<String> rtn = new ArrayList<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 5)
                continue;
            double score = new Double(tokens[4]);
            if (Math.abs(score) >= threshold)
                rtn.add(tokens[0]);
        }
        fu.close();
        return rtn;
    }
    
    /**
     * Load ppi to maximum interaction effect score using the default output filename.
     * @return
     * @throws IOException
     */
    public Map<String, Double> loadPPIToMaxScore() throws IOException {
        return loadPPIToMaxScore(outputFileName);
    }
    
    @Test
    public void checkAllHumanReactions() throws Exception {
        String output = "results/ReactionsInMechisom_051017.txt";
        CancerDriverReactomeAnalyzer reactomeAnalyzer = new CancerDriverReactomeAnalyzer();
        
        checkAllHumanReactions(reactomeAnalyzer, output);
    }
    
    /**
     * Check all human reactions based on mechismo data.
     * @throws Exception
     */
    public void checkAllHumanReactions(CancerDriverReactomeAnalyzer reactomeAnalyzer,
                                       String outputFileName) throws Exception {
        // Load all the non-disease reactions from Reactome
        List<GKInstance> reactions = reactomeAnalyzer.loadHumanReactions();
        System.out.println("Total human reactions in Reactome: " + reactions.size());

        Map<String, Double> ppiToScore = loadPPIToMaxScore();
        System.out.println("Total PPIs in mechismo: " + ppiToScore.size());
        
        ReactomeAnalyzer reactomeDataAnalyzer = new ReactomeAnalyzer();
        fu.setOutput(outputFileName);
        fu.printLine("DB_ID\tName\tTotal_FIs\tMechismo_FIs\tMax_Score");
        int total = 0;
        int checkedReactions = 0;
        for (GKInstance reaction : reactions) {
            Set<String> fis = reactomeDataAnalyzer.generateTentativePPIsForReaction(reaction,
                                                                                    false);
            if (fis.size() == 0)
                continue;
            // Check if there is a structure available from interactome3d
            int overlap = 0;
            Double score = null;
            for (String fi : fis) {
                if (ppiToScore.containsKey(fi)) {
                    overlap ++;
                   Double currentScore = ppiToScore.get(fi);
                   if (score == null || Math.abs(currentScore) > Math.abs(score))
                       score = currentScore;
                }
            }
            fu.printLine(reaction.getDBID() + "\t" +
                    reaction.getDisplayName() + "\t" +
                    fis.size() + "\t" +
                    overlap + "\t" +
                    (score == null ? "" : score));
            if (overlap > 0)
                total ++;
            checkedReactions ++;
        }
        fu.close();
        System.out.println("Total checked reactions: " + checkedReactions);
        System.out.println("Total reactions having FIs in mechsimo: " + total);
        // Output ran on May 10, 2017
//        Total human reactions in Reactome: 9453
//        Total PPIs in mechismo: 13323
//        Total checked reactions: 4589
//        Total reactions having FIs in mechsimo: 1265
    }
    
    @Test
    public void testLoadPPIToMaxScore() throws IOException {
        Map<String, Double> ppiToScore = loadPPIToMaxScore();
        System.out.println("Total PPIs: " + ppiToScore.size());
        // The above output should be: Total PPIs: 13323
    }
    
    @Test
    public void checkKnownCancerDriverGenes() throws IOException {
        CancerDriverAnalyzer driverAnalyzer = new CancerDriverAnalyzer();
        Set<String> knownCancerGenes = driverAnalyzer.getDriverGenes(null);
        System.out.println("Total known driver genes: " + knownCancerGenes.size());
        
    }
    
    /**
     * Load the largest functional impact scores for mutations listed in cosmic between two proteins
     * interacting each other. The functional impact scores for mutations are reported in column 27
     * ( ie: interaction effect for all contacts of the site with this interactor).
     * @param mechismoOutputFilePath
     * @param mechismoInteractionScorePattern
     * @return
     * @throws IOException
     */
    public Map<String, Double> loadPPIToMaxScore(String mechismoFileName) throws IOException{
        fu.setInput(mechismoFileName);
        String line = null;
        Map<String, Double> ppiToScore = new HashMap<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 28)
                continue;
            // If no second protein is reported, escape it.
            if (tokens[19].trim().length() == 0)
                continue;
            String ppi = InteractionUtilities.generateFIFromGene(tokens[1], 
                                                                 tokens[19]);
            Double score = ppiToScore.get(ppi);
            Double currentScore = new Double(tokens[27]);
            if (score == null || Math.abs(currentScore) > Math.abs(score))
                ppiToScore.put(ppi, currentScore);
        }
        fu.close();
        return ppiToScore;
    }
    
    @Test
    public void testGetMechismoInteractions() throws IOException {
        Set<String> interactions = getInteractions(outputFileName);
        System.out.println("Total interactions: " + interactions.size());
        // Output is: 13323
    }
    
    @Test
    public void checkOverlapWithInteractome3d() throws IOException {
        Set<String> mechismoPPIs = getInteractions(outputFileName);
        System.out.println("Total PPIs in mechismo output: " + mechismoPPIs.size());
        fu.saveInteractions(mechismoPPIs, "results/MechismoPPIs_051017.txt");
        
        Interactome3dAnalyzer interactomeAnalyser = new Interactome3dAnalyzer();
        String interactomeDirName = "datasets/interactome3d/2016_06/prebuilt/representative/";
        Map<String, File> fiToPDB = interactomeAnalyser.loadPPIToPDBFile(interactomeDirName,
                                                                         false);
        System.out.println("Total PPIs in interactome3d: " + fiToPDB.size());
        fu.saveInteractions(fiToPDB.keySet(), "results/Interactome3dPPIs_051017.txt");
        
        Set<String> shared = InteractionUtilities.getShared(mechismoPPIs, fiToPDB.keySet());
        System.out.println("Shared: " + shared.size());
        // Output:
        //        Total PPIs in mechismo output: 13,323
        //        Total PPIs in interactome3d: 9,875
        //        Shared: 4011
        // TODO: The 50% coverage of interactome3d by mechismo may result from no cosmic mutations are found in other
        // PPIs. We should check the original cosmic mutation data and then compare the outputs.
    }
    
    /**
     * Get the interactions reported in the mechismo output file.
     * @param mechismoOutputFilePath
     * @param mechismoInteractionScorePattern
     * @return
     * @throws IOException
     */
    public Set<String> getInteractions(String fileName) throws IOException {
        Set<String> interactions = new HashSet<>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 20)
                continue;
            if (tokens[19].trim().length() == 0)
                continue;
            String ppi = InteractionUtilities.generateFIFromGene(tokens[1], tokens[19]);
            interactions.add(ppi);
        }
        fu.close();
        return interactions;
    }
    
    /**
     * The mechismo scores are calculated for each mutation in a specific coordinate in a specific protein. It is 
     * collected from scores for this specific mutation across all interactions invovled in the protein having 
     * the mutation. So they are unique regard mutation, coordinate, and protein. However, scores used for cancer
     * analysis are not: the interacting partners need to be considered. In other words,
     * they are specific to mutations, coordinates, and two proteins involved in the interaction.
     * @throws IOException
     */
    @Test
    public void checkMechismoScores() throws IOException {
        String fileName = dirName + "SOS1_RAS.txt";
        fileName = outputFileName;
        
        fu.setInput(fileName);
        String line = null;
        Map<String, String> mutationToScore = new HashMap<>();
        Set<String> scores = new HashSet<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            
            //            for (int i = 0; i < tokens.length; i++) {
            //                System.out.println(i + "\t" + tokens[i]);
            //            }
            //            if (true)
            //                break;
            
            if (tokens.length < 20) {
                //                System.out.println(line);
                continue;
            }
            // Make sure there are AA mutations
            if (tokens[6].length() > 0 && tokens[19].length() > 0) {
                //                mutationToScore.put(tokens[6] + "\t" + tokens[19], tokens[14]);
                mutationToScore.put(tokens[6], tokens[14]);
                scores.add(tokens[14]);
            }
        }
        System.out.println("Total mutations: " + mutationToScore.size()); // Should be 20 * 19 = 380 
        System.out.println("Total scores: " + new HashSet<>(mutationToScore.values()).size());
        System.out.println("Actual scores: " + scores.size());
        for (String mutation : mutationToScore.keySet())
            System.out.println(mutation + "\t" + mutationToScore.get(mutation));
        System.out.println("\nActual extra scores: ");
        scores.removeAll(mutationToScore.values());
        for (String score : scores)
            System.out.println(score);
    }
    
    @Test
    public void pickRows() throws IOException {
        String[] targetGenes = new String[] {
                "HRAS", "NRAS", "KRAS", "SOS1"
        };
        Set<String> geneSet = new HashSet<>();
        for (String gene : targetGenes)
            geneSet.add(gene);
        
        String fileName = outputFileName;
        //        fileName = dirName + "HRAS_mechismo_op.tsv";
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 19)
                continue; // Cannot get interactors
            //            System.out.println(tokens[0] + "\t" + tokens[18]);
            if (geneSet.contains(tokens[0]) && geneSet.contains(tokens[18]))
                System.out.println(line);
        }
        fu.close();
    }
    
}
