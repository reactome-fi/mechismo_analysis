package org.reactome.cancer.driver;

import java.io.File;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.reactome.r3.ReactomeAnalyzer;
import org.reactome.structure.model.MutationObservation;
import org.reactome.structure.model.ReactionMutationProfile;

/**
 * Perform cancer driving analysis for Reactome reactions.
 * @author wug
 *
 */
public class ReactomeReactionAnalyzer {
    private ReactomeAnalyzer reactomeDataAnalyzer;
    private Map<String, Set<File>> ppiToPDBs;
    private Map<String, String> uniprotIdToGene;
    private Map<String, Set<MutationObservation>> geneToMutations;
    
    public ReactomeReactionAnalyzer() {
    }
    
    private boolean validatePreCondition() {
        if (reactomeDataAnalyzer == null)
            return false;
        if (ppiToPDBs == null || ppiToPDBs.size() == 0)
            return false;
        if (uniprotIdToGene == null || uniprotIdToGene.size() == 0)
            return false;
        if (geneToMutations == null || geneToMutations.size() == 0)
            return false;
        return true;
    }
    
    public ReactomeAnalyzer getReactomeDataAnalyzer() {
        return reactomeDataAnalyzer;
    }

    public void setReactomeDataAnalyzer(ReactomeAnalyzer reactomeDataAnalyzer) {
        this.reactomeDataAnalyzer = reactomeDataAnalyzer;
    }

    public Map<String, Set<File>> getPpiToPDBs() {
        return ppiToPDBs;
    }

    public void setPpiToPDBs(Map<String, Set<File>> ppiToPDBs) {
        this.ppiToPDBs = ppiToPDBs;
    }
    
    public Set<String> extractFIs(GKInstance reaction) throws Exception {
        if (reactomeDataAnalyzer == null)
            throw new IllegalStateException("Set ReactomeAnalyzer for the object!");
        Set<String> fis = reactomeDataAnalyzer.generateTentativePPIsForReaction(reaction, false);
        return fis;
    }

    /**
     * Perform mutation interface analysis for a single reaction.
     * @param reaction
     * @return
     * @throws Exception
     */
    public ReactionMutationProfile analyzeReaction(GKInstance reaction,
                                                   Set<String> fis) throws Exception {
        if (!validatePreCondition())
            throw new IllegalStateException("Make sure needed properties have been set for this object!");;
        if (fis == null || fis.size() == 0)
            return null;
        ReactionMutationProfile mutationProfile = new ReactionMutationProfile();
        mutationProfile.setExtractFIs(fis);
        fis.forEach(fi -> {
            Set<File> pdbs = ppiToPDBs.get(fi);
            
        });
//        Set<String> sigGenes = checkInterfaces(accessions,
//                uniprotIdToGene,
//                geneToMutations,
//                pdbFiles);
//        System.out.println(reaction.getDBID() + "\t" +
//                reaction.getDisplayName() + "\t" +
//                sigGenes.size() + "\t" +
//                sigGenes);
        return null;
    }

}
