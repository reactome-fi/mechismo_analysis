package org.reactome.cancer.driver;

import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.annotate.AnnotationHelper;
import org.reactome.annotate.AnnotationType;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.r3.CosmicAnalyzer;
import org.reactome.r3.ReactomeAnalyzer;
import org.reactome.structure.model.ProteinMutation;
import org.reactome.structure.model.ReactionMutationProfile;

/**
 * Check Reactome pathways for known cancer driver genes.
 * @author wug
 *
 */
public class ReactomePathwayAnalyzer extends CancerDriverReactomeAnalyzer {
    private final static Logger logger = Logger.getLogger(ReactomePathwayAnalyzer.class);
    
    public ReactomePathwayAnalyzer() {
    }
    
    /**
     * Analyze reactions in a specific pathway.
     * @throws Exception
     */
    @Test
    public void analyzePathway() throws Exception {
        // Pick up a pathway for detailed analysis
        Long dbId = 8848021L; // Signaling by PTK6: 66 genes in total, 22 are cancer drivers
        dbId = 1433557L; // SCF-KIT signaling: 290 genes in total, 44 cancer drivers
        
        MySQLAdaptor dba = getDBA();
        GKInstance pathway = dba.fetchInstance(dbId);
        logger.info(pathway);
        
        Set<GKInstance> reactions = InstanceUtilities.grepPathwayEventComponents(pathway);
        logger.info("Total reactions: " + reactions.size());
        
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        Set<String> pathwayGenes = reactomeAnalyzer.grepGenesFromPathway(pathway);
        logger.info("Total genes: " + pathwayGenes.size());
        
        CosmicAnalyzer cosmicAnalyzer = new CosmicAnalyzer();
        Map<String, List<ProteinMutation>> geneToMutations = cosmicAnalyzer.loadMutations(pathwayGenes);
        logger.info("Total mutation size: " + geneToMutations.size());
        
        ReactomeReactionAnalyzer reactionAnalyzer = new ReactomeReactionAnalyzer();
        reactionAnalyzer.setGeneToMutations(geneToMutations);
        
        for (GKInstance reaction : reactions) {
            System.out.println();
            System.out.println(reaction);
            ReactionMutationProfile profile = reactionAnalyzer.analyzeReaction(reaction);
            if (profile == null)
                continue;
            String output = profile.exportResults();
            if (output == null)
                System.out.println("No information!");
            else
                System.out.println(output);
        }
    }
    
    @Test
    public void performEnrichmentAnalysis() throws Exception {
        CancerDriverAnalyzer driverAnalyzer = new CancerDriverAnalyzer();
        Set<String> knownDrivers = driverAnalyzer.getDriverGenes(null);
        
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        AnnotationHelper helper = new AnnotationHelper();
        annotator.setAnnotationHelper(helper);
        String geneToPathwayFile = "resources/ProteinNameToTopics022717.txt";
        helper.setProteinNameToPathwayFile(geneToPathwayFile);
        
        List<GeneSetAnnotation> annotations = annotator.annotateGenesWithFDR(knownDrivers, AnnotationType.Pathway);
        System.out.println("Pathway\tNumberInPathway\tHitNumber\tRatioOfPathway\tpValue\tFDR\tHitIds");
        annotations.stream()
                   .filter(annotation -> annotation.getTopic().endsWith("(R)"))
                   .forEach(System.out::println);
    }

}
