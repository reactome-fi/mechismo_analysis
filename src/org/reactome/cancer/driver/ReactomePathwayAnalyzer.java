package org.reactome.cancer.driver;

import java.util.List;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.annotate.AnnotationHelper;
import org.reactome.annotate.AnnotationType;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;

/**
 * Check Reactome pathways for known cancer driver genes.
 * @author wug
 *
 */
public class ReactomePathwayAnalyzer extends CancerDriverReactomeAnalyzer {
    
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
        MySQLAdaptor dba = getDBA();
        GKInstance pathway = dba.fetchInstance(dbId);
        Set<GKInstance> reactions = InstanceUtilities.grepPathwayEventComponents(pathway);
        System.out.println(pathway + ": " + reactions.size());
        
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
