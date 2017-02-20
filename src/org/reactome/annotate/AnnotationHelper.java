/*
 * Created on Jun 22, 2010
 *
 */
package org.reactome.annotate;

import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to help annotation done in other classes in this package.
 * @author wgm
 *
 */
public class AnnotationHelper {
    //private String proteinNameToPathwayFile = R3Constants.RESULT_DIR + "ProteinNameToTopics051109.txt";
    private String proteinNameToPathwayFile = R3Constants.RESULT_DIR + "ProteinNameToTopics072512.txt";
    private String geneFIFile = R3Constants.GENE_FI_FILE_NAME;
    private FileUtility fu;
    private GOTermLoader goTermLoader;
    // For Reactome pathways only
    private String proteinNameToReactomePathwayFile;
    private String reactionIdToPathwayFile;
    
    public AnnotationHelper() {
        fu = new FileUtility();
    }
    
    public String getReactionIdToPathwayFile() {
        return reactionIdToPathwayFile;
    }

    public void setReactionIdToPathwayFile(String reactionIdToPathwayFile) {
        this.reactionIdToPathwayFile = reactionIdToPathwayFile;
    }

    public GOTermLoader getGoTermLoader() {
        return goTermLoader;
    }

    public void setGoTermLoader(GOTermLoader goTermLoader) {
        this.goTermLoader = goTermLoader;
    }

    public String getGeneFIFile() {
        return geneFIFile;
    }
    
    public void setGeneFIFile(String geneFIFile) {
        this.geneFIFile = geneFIFile;
    }
    
    public String getProteinNameToReactomePathwayFile() {
        return proteinNameToReactomePathwayFile;
    }

    public void setProteinNameToReactomePathwayFile(String proteinNameToReactomePathwayFile) {
        this.proteinNameToReactomePathwayFile = proteinNameToReactomePathwayFile;
    }

    public void setProteinNameToPathwayFile(String fileName) {
        this.proteinNameToPathwayFile = fileName;
    }
    
    public String getProteinNameToPathwayFile() {
        return this.proteinNameToPathwayFile;
    }
    
    public Map<String, Set<String>> loadProteinNameToPathwaysMap() throws IOException {
        return fu.loadSetMap(proteinNameToPathwayFile);
    }
    
    public Map<String, Set<String>> loadProteinNameToReactomePathwaysMap() throws IOException {
        if (proteinNameToReactomePathwayFile == null)
            return null;
        return fu.loadSetMap(proteinNameToReactomePathwayFile);
    }
    
    public Map<String, Set<String>> loadReactomeReactionToPathwaysMap() throws IOException {
        if (reactionIdToPathwayFile == null)
            return null;
        return fu.loadSetMap(reactionIdToPathwayFile);
    }
    
    public Map<String, Set<String>> loadProteinNameToTermsMap(AnnotationType type) throws IOException {
        Map<String, Set<String>> proteinToTerms = null;
        boolean isGOTerm = true;
        switch (type) {
            case CC :
                proteinToTerms = goTermLoader.loadProteinToGOCCTerms();
                break;
            case BP :
                proteinToTerms = goTermLoader.loadProteinToGOBPTerms();
                break;
            case MF :
                proteinToTerms = goTermLoader.loadProteinToGOMFTerms();
                break;
            default :
                proteinToTerms = loadProteinNameToPathwaysMap(); // Default using pathways
                isGOTerm = false;
        }
        if (isGOTerm) {
            // Need to convert protein in UniProt to gene names
            proteinToTerms = goTermLoader.convertProteinIdToNameForGO(proteinToTerms);
            // Convert GO Term ids to names
            proteinToTerms = goTermLoader.convertGOIdsToTerms(proteinToTerms);
        }
        return proteinToTerms;
    }
    
    public Set<String> loadRandomGenes() throws IOException {
        FileUtility fu = new FileUtility();
        Set<String> fis = fu.loadInteractions(geneFIFile); 
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> randomGenes = new HashSet<String>(allGenes);
        return randomGenes;
    }
    
    public Map<String, Integer> countProteinsInTopics(Map<String, Set<String>> protein2Topic) {
        Map<String, Integer> topicToProteinNumber = new HashMap<String, Integer>();
        for (Iterator<String> it = protein2Topic.keySet().iterator(); it.hasNext();) {
            String id = it.next();
            Set<String> topics = protein2Topic.get(id);
            for (String t : topics) {
                Integer c = topicToProteinNumber.get(t);
                if (c == null)
                    topicToProteinNumber.put(t, 1);
                else {
                    topicToProteinNumber.put(t, ++c);
                }
            }
        }
        return topicToProteinNumber;
    }
    
    public Map<String, Double> calculateTopicToRatio(int total, 
                                                     Map<String, Integer> topicToIdNumber) {
        Map<String, Double> topicToRatio = new HashMap<String, Double>();
        for (Iterator<String> it = topicToIdNumber.keySet().iterator(); it.hasNext();) {
            String topic = it.next();
            Integer number = topicToIdNumber.get(topic);
            Double ratio = (double) number / total;
            topicToRatio.put(topic, ratio);
        }
        return topicToRatio;
    }   
    
    public void sortGeneSetAnnotation(List<GeneSetAnnotation> annotations) {
        Collections.sort(annotations, new Comparator<GeneSetAnnotation>() {
            public int compare(GeneSetAnnotation info1,
                               GeneSetAnnotation info2) {
                double diff = info1.getPValue() - info2.getPValue();
                if (diff > 0)
                    return 1;
                if (diff < 0)
                    return -1;
                return 0;
            }
        });
    } 
    
//    @Test
//    public void testLoadProteinNameToPathwayMap() throws IOException {
//        Map<String, Set<String>> geneToPathways = loadProteinNameToPathwaysMap();
//        Map<String, Set<String>> pathwayToGenes = InteractionUtilities.switchKeyValues(geneToPathways);
//        System.out.println("Pathway\tProtein_Number\tProteins");
//        List<String> pathwayList = new ArrayList<String>(pathwayToGenes.keySet());
//        Collections.sort(pathwayList);
//        for (String pathway : pathwayList) {
//            if (pathway.endsWith("(N)")) {
//                Set<String> genes = pathwayToGenes.get(pathway);
//                List<String> geneList = new ArrayList<String>(genes);
//                Collections.sort(geneList);
//                String geneText = InteractionUtilities.joinStringElements(",", geneList);
//                System.out.println(pathway + "\t" + 
//                                   geneList.size() + "\t" +
//                                   geneText);
//            }
//        }
//    }
    
}
