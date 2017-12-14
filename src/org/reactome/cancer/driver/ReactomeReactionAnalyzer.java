package org.reactome.cancer.driver;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.cancer.InteractionInterfaceAnalyzer;
import org.reactome.r3.CosmicAnalyzer;
import org.reactome.r3.Interactome3dAnalyzer;
import org.reactome.r3.ReactomeAnalyzer;
import org.reactome.structure.model.InteractionMutationProfile;
import org.reactome.structure.model.ProteinMutation;
import org.reactome.structure.model.ReactionMutationProfile;

/**
 * Perform cancer driving analysis for Reactome reactions.
 * @author wug
 *
 */
public class ReactomeReactionAnalyzer {
    private final static Logger logger = Logger.getLogger(ReactomeReactionAnalyzer.class);
    private ReactomeAnalyzer reactomeDataAnalyzer;
    private Map<String, List<File>> ppiToPDBs;
    private Map<String, List<ProteinMutation>> geneToMutations;
    
    public ReactomeReactionAnalyzer() {
    }
    
    private boolean validatePreCondition() throws IOException {
        if (reactomeDataAnalyzer == null)
            reactomeDataAnalyzer = new ReactomeAnalyzer();
        if (ppiToPDBs == null || ppiToPDBs.size() == 0) {
            // We will use the default interactome3d data in case this is not there
            loadInteractome3d();
//            throw new IllegalStateException("PPIToPDBs map has not been loaded!");
        }
        if (geneToMutations == null || geneToMutations.size() == 0)
            throw new IllegalStateException("geneToMutations map has not been assidnged!");
        return true;
    }
    
    private void loadInteractome3d() throws IOException {
        Interactome3dAnalyzer analyzer = new Interactome3dAnalyzer();
        String dirName = "datasets/interactome3d/2017_01/complete/";
        ppiToPDBs = analyzer.loadPPIToCompletePDBFiles(dirName, 
                                                       true);
    }
    
    public ReactomeAnalyzer getReactomeDataAnalyzer() {
        return reactomeDataAnalyzer;
    }

    public void setReactomeDataAnalyzer(ReactomeAnalyzer reactomeDataAnalyzer) {
        this.reactomeDataAnalyzer = reactomeDataAnalyzer;
    }

    public Map<String, List<File>> getPpiToPDBs() {
        return ppiToPDBs;
    }

    public void setPpiToPDBs(Map<String, List<File>> ppiToPDBs) {
        this.ppiToPDBs = ppiToPDBs;
    }
    
    public void setGeneToMutations(Map<String, List<ProteinMutation>> geneToMutations) {
        this.geneToMutations = geneToMutations;
    }
    
    public Set<String> extractFIs(GKInstance reaction) throws Exception {
        if (reactomeDataAnalyzer == null)
            throw new IllegalStateException("Set ReactomeAnalyzer for the object!");
        Set<String> fis = reactomeDataAnalyzer.generateTentativePPIsForReaction(reaction, true);
        return fis;
    }
    
    @Test
    public void testExtractFIs() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                "reactome_59_plus_i",
                "root",
                "macmysql01");
        reactomeDataAnalyzer = new ReactomeAnalyzer();
        reactomeDataAnalyzer.setMySQLAdaptor(dba);
        
        Long dbId = 5672969L;
        dbId = 171026L;
        
        GKInstance reaction = dba.fetchInstance(dbId);
        Set<String> fis = extractFIs(reaction);
        System.out.println("Total FIs: " + fis.size());
        fis.forEach(System.out::println);
    }
    
    @Test
    public void testAnalyzeReaction() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                             "reactome_59_plus_i",
                                             "root",
                                             "macmysql01");
        
        // [Reaction:8849032] PTK6 phosphorylates KHDRBS2
        Long dbId = 8849032L;
        //5658435 RAS GAPs bind RAS:GTP
        dbId = 5658435L; // We should not see null from this reaction
        // There is nothing displayed for this reaction
//        dbId = 8848751L; // [Reaction:8848751] PTK6 binds AKT1
        dbId = 5675431L; // [Reaction:5675431] PP2A dephosphorylates RAF1
        
        GKInstance reaction = dba.fetchInstance(dbId);
        
        System.out.println("Analyzing " + reaction + "...");
        reactomeDataAnalyzer = new ReactomeAnalyzer();
        Set<String> reactionGenes = reactomeDataAnalyzer.grepGenesFromReaction(reaction);
        System.out.println("Total genes in reaction: " + reactionGenes.size());
        // Initialize conditions
        CosmicAnalyzer cosmicAnalyzer = new CosmicAnalyzer();
        Map<String, List<ProteinMutation>> geneToCosmicEntries = cosmicAnalyzer.loadMutations(reactionGenes);
        System.out.println("Total mutations loaded: " + geneToCosmicEntries.size());
        setGeneToMutations(geneToCosmicEntries);
        
        ReactionMutationProfile profile = analyzeReaction(reaction);
        System.out.println(profile.exportResults());
    }

    /**
     * Perform mutation interface analysis for a single reaction.
     * @param reaction
     * @return
     * @throws Exception
     */
    public ReactionMutationProfile analyzeReaction(GKInstance reaction) throws Exception {
        if (!validatePreCondition())
            throw new IllegalStateException("Make sure needed properties have been set for this object!");;
        Set<String> fis = reactomeDataAnalyzer.generateTentativePPIsForReaction(reaction, false);
        if (fis == null || fis.size() == 0) {
            logger.info(reaction + " cannot generate any FIs!");
            return null;
        }
        logger.info("Extracted FIs: " + fis.size());
        InteractionInterfaceAnalyzer interfaceAnalyzer = new InteractionInterfaceAnalyzer();
        interfaceAnalyzer.setGeneToMutations(geneToMutations);
        ReactionMutationProfile mutationProfile = new ReactionMutationProfile();
        mutationProfile.setExtractFIs(fis);
        fis.forEach(fi -> {
            List<File> pdbs = ppiToPDBs.get(fi);
            if (pdbs == null || pdbs.size() == 0) {
                logger.info(fi + ": No PDB file!");
                return;
            }
            // Want to pick up the result with minimum p-value
            InteractionMutationProfile intMutProfile = null;
            for (File pdbFile : pdbs) {
                try {
                    InteractionMutationProfile tmp = interfaceAnalyzer.analyze(pdbFile);
                    if (intMutProfile == null || intMutProfile.getP_value() > tmp.getP_value())
                        intMutProfile = tmp;
                }
                catch(Exception e) {
                    logger.error(e.getMessage(), e);
                }
            }
            mutationProfile.addPpiToProfile(fi, intMutProfile);
        });
        return mutationProfile;
    }

}
