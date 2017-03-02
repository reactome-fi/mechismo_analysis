/*
 * Created on May 3, 2006
 *
 */
package org.reactome.r3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.PersistenceAdaptor;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.InvalidAttributeException;
import org.gk.schema.SchemaAttribute;
import org.gk.schema.SchemaClass;
import org.junit.Test;
import org.reactome.r3.util.InteractionUtilities;

@SuppressWarnings("unchecked")
public class ReactomeAnalyzer {
    private boolean excludeComplex = false;
    protected PersistenceAdaptor dba;
    // A helper class
    private ReactomeAnalyzerTopicHelper topicHelper;
    
    public ReactomeAnalyzer() {
        topicHelper = new ReactomeAnalyzerTopicHelper();
    }
    
    public void setExcludeComplex(boolean value) {
        this.excludeComplex = value;
    }
    
    public void setMySQLAdaptor(MySQLAdaptor dba) {
        this.dba = dba;
    }
    
    public PersistenceAdaptor getMySQLAdaptor() throws Exception {
        return this.dba;
    }
    
    @SuppressWarnings("rawtypes")
    public void extractInteractorsFromReaction(GKInstance rxn, 
                                               Set<GKInstance> interactors) throws Exception {
        List input = rxn.getAttributeValuesList(ReactomeJavaConstants.input);
//        // Something special for gene regulatory reaction annotated in BlackBoxEvent
//        boolean isGeneRegulatory = false;
//        if (rxn.getSchemClass().isa(ReactomeJavaConstants.BlackBoxEvent) && input != null && input.size() == 1) {
//            GKInstance input1 = (GKInstance) input.get(0);
//            if (input1.getSchemClass().isValidAttribute(ReactomeJavaConstants.referenceEntity)) {
//                GKInstance refEntity = (GKInstance) input1.getAttributeValue(ReactomeJavaConstants.referenceEntity);
//                if (refEntity != null && refEntity.getSchemClass().isa(ReactomeJavaConstants.ReferenceDNASequence)) {
//                    isGeneRegulatory = true;
//                    GKInstance output = (GKInstance) rxn.getAttributeValue(ReactomeJavaConstants.output);
//                    if (output != null)
//                        interactors.add(output);
//                }
//            }
//        }
//        if (!isGeneRegulatory && input != null) {
//            interactors.addAll(input);
//        }
        if (input != null)
            interactors.addAll(input);
        // Get catalyst
        List cas = rxn.getAttributeValuesList(ReactomeJavaConstants.catalystActivity);
        if (cas != null) {
            for (Iterator it1 = cas.iterator(); it1.hasNext();) {
                GKInstance ca = (GKInstance) it1.next();
                List catalysts = ca.getAttributeValuesList(ReactomeJavaConstants.physicalEntity);
                if (catalysts != null)
                    interactors.addAll(catalysts);
            }
        }
        // Check regulators
        Collection regulations = rxn.getReferers(ReactomeJavaConstants.regulatedEntity);
        if (regulations != null) {
            for (Iterator it1 = regulations.iterator(); it1.hasNext();) {
                GKInstance regulation = (GKInstance) it1.next();
                List regulators = regulation.getAttributeValuesList(ReactomeJavaConstants.regulator);
                for (Iterator it2 = regulators.iterator(); it2.hasNext();) {
                    GKInstance regulator = (GKInstance) it2.next();
                    if (regulator.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity))
                        interactors.add(regulator);
                }
            }
        }
    }
    
    protected void generateInteractionsWithComplexAndSet(Set<GKInstance> interactors,
                                                         Set<String> interactions,
                                                         GKInstance reaction) throws Exception {
        List<GKInstance> list = new ArrayList<GKInstance>(interactors);
        int size = list.size();
        for (int i = 0; i < size - 1; i++) {
            GKInstance interactor1 = list.get(i);
            for (int j = i + 1; j < size; j++) {
                GKInstance interactor2 = list.get(j);
                generateInteractionsWithComplexAndSet(interactor1, 
                                                      interactor2, 
                                                      interactions);
            }
        }
    }
    
    private void generateInteractionsWithComplexAndSet(GKInstance interactor1,
                                                       GKInstance interactor2,
                                                       Set<String> interactions) throws Exception {
        Set<GKInstance> refPepSeqs1 = grepRefPepSeqs(interactor1);
        if (refPepSeqs1.size() == 0)
            return;
        Set<GKInstance> refPepSeqs2 = grepRefPepSeqs(interactor2);
        if (refPepSeqs2.size() == 0)
            return;
        String int1 = convertRefPepSeqToString(refPepSeqs1, interactor1);
        if (int1.length() == 0)
            return;
        String int2 = convertRefPepSeqToString(refPepSeqs2, interactor2);
        if (int2.length() == 0)
            return;
        int compare = int1.compareTo(int2);
        if (compare < 0)
            interactions.add(int1 + " " + int2);
        else
            interactions.add(int2 + " " + int1);
    }
              
    private String convertRefPepSeqToString(Set<GKInstance> refPepSeqs,
                                            GKInstance interactor) throws Exception {
        List<String> ids = new ArrayList<String>();
        for (Iterator<GKInstance> it = refPepSeqs.iterator(); it.hasNext();) {
            GKInstance refPepSeq = it.next();
            String identifier = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
            if (identifier == null)
                continue; // maybe
            ids.add(identifier);
        }
        Collections.sort(ids);
        String delimit = "?"; // default: should not be used
        if (interactor.getSchemClass().isa(ReactomeJavaConstants.Complex))
            delimit = ",";
        else if (interactor.getSchemClass().isa(ReactomeJavaConstants.EntitySet))
            delimit = "|";
        StringBuilder builder = new StringBuilder();
        for (Iterator<String> it = ids.iterator(); it.hasNext();) {
            builder.append(it.next());
            if (it.hasNext())
                builder.append(delimit);
        }
        return builder.toString();
    }
    
    private Set<GKInstance> grepComplexes(Set<GKInstance> interactors) throws Exception {
        Set<GKInstance> complexes = new HashSet<GKInstance>();
        for (GKInstance interactor : interactors) {
            if (interactor.getSchemClass().isa(ReactomeJavaConstants.Complex))
                complexes.add(interactor);
            Set<GKInstance> temp = InstanceUtilities.getContainedInstances(interactor,
                                                                           ReactomeJavaConstants.hasMember,
                                                                           ReactomeJavaConstants.hasComponent);
            for (GKInstance inst : temp) {
                if (inst.getSchemClass().isa(ReactomeJavaConstants.Complex))
                    complexes.add(inst);
            }
        }
        return complexes;
    }
    
    private Set<String> grepGenesFromComplex(GKInstance complex) throws Exception {
        Set<String> genes = new HashSet<String>();
        grepGenesFromEntity(complex, genes);
        return genes;
    }
    
    private void grepGenesFromEntity(GKInstance pe,
                                     Set<String> genes) throws Exception {
        Set<GKInstance> refEntities = InstanceUtilities.grepReferenceEntitiesForPE(pe);
        for (GKInstance refEntity : refEntities) {
            if (refEntity.getSchemClass().isValidAttribute(ReactomeJavaConstants.geneName)) {
                String geneName = (String) refEntity.getAttributeValue(ReactomeJavaConstants.geneName);
                if (geneName != null)
                    genes.add(geneName);
            }
        }
    }

    private Set<String> grepGenesFromReaction(GKInstance rxn) throws Exception {
        Set<GKInstance> participants = InstanceUtilities.getReactionParticipants(rxn);
        Set<String> genes = new HashSet<String>();
        for (GKInstance participant : participants) {
            grepGenesFromEntity(participant, genes);
        }
        return genes;
    }
    
    @Test
    public void generateFIsForOneReaction() throws Exception {
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        Long dbId = 5672965L; // RAS GEFs promote RAS nucleotide exchange
//        dbId = 5672972L; // MAP2Ks and MAPKs bind to the activated RAF complex
//        dbId = 69213L; // Formation of Cyclin D:Cdk4/6 complexes 
//        dbId = 5617896L; // Retinoic acid activates HOXD4 chromatin
        String dirName = "/Users/gwu/Documents/EclipseWorkspace/caBigR3/results/DriverGenes/Drivers_0816/";
        String fileName = dirName + "FIsInReaction" + dbId + ".txt";
        generateFIsForOneReaction(dba,
                                  dbId,
                                  fileName);
    }
    
    private void generateFIsForOneReaction(MySQLAdaptor dba,
                                           Long dbId,
                                           String fileName) throws Exception, IOException {
        List<Long> dbIds = new ArrayList<Long>();
        dbIds.add(dbId);
        generateFIsForReactionsWithFeatures(dba, dbIds, fileName);
    }
    
    /**
     * Extract a set of FIs in gene names from the passed reaction.
     * @param reaction
     * @return
     * @throws Exception
     */
    public Set<String> generateAttentativePPIsForReaction(GKInstance rxn,
                                                          boolean useGeneName) throws Exception {
        Set<String> rxtFIs = new HashSet<String>();
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        // We will have two steps to generate FIs that most likely interact physically.
        // The first step is for FIs among inputs and catalysts
        List input = rxn.getAttributeValuesList(ReactomeJavaConstants.input);
        if (input != null)
            interactors.addAll(input);
        // Get catalyst
        List cas = rxn.getAttributeValuesList(ReactomeJavaConstants.catalystActivity);
        if (cas != null) {
            for (Iterator it1 = cas.iterator(); it1.hasNext();) {
                GKInstance ca = (GKInstance) it1.next();
                List catalysts = ca.getAttributeValuesList(ReactomeJavaConstants.physicalEntity);
                if (catalysts != null)
                    interactors.addAll(catalysts);
            }
        }
        if (interactors.size() == 0)
            return rxtFIs; // Nothing needs to be done
        if (interactors.size() > 1) {
            generateInteractions(interactors,
                                 rxtFIs, 
                                 rxn,
                                 useGeneName);
        }
        // The second step is for FIs between inputs and catalysts and regulators. This is different from
        // FIs generation, where FIs have also been extracted between two regulators working on the same
        // reaction.
        Collection regulations = rxn.getReferers(ReactomeJavaConstants.regulatedEntity);
        if (regulations == null || regulations.size() == 0)
            return rxtFIs;
        for (Iterator it1 = regulations.iterator(); it1.hasNext();) {
            GKInstance regulation = (GKInstance) it1.next();
            List regulators = regulation.getAttributeValuesList(ReactomeJavaConstants.regulator);
            for (Iterator it2 = regulators.iterator(); it2.hasNext();) {
                GKInstance regulator = (GKInstance) it2.next();
                if (regulator.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity)) {
                    for (GKInstance interactor : interactors) {
                        if (regulator == interactor)
                            continue;
                        generateInteractions(regulator, 
                                             interactor,
                                             rxtFIs,
                                             rxn,
                                             useGeneName);
                    }
                }
            }
        }
        return rxtFIs;
    }
    
    @Test
    public void testGenerateFIsForReaction() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_59_plus_i",
                                            "root",
                                            "macmysql01");
//        Long dbId = 5672965L;
        Long dbId = 5617454L;
        dbId = 2730888L;
        GKInstance reaction = dba.fetchInstance(dbId);
        Set<String> fis = generateAttentativePPIsForReaction(reaction, false);
        System.out.println("Total fis in reaction " + reaction + ": " + fis.size());
        for (String fi : fis)
            System.out.println(fi);
    }
    
    public void generateFIsForReactionsWithFeatures(MySQLAdaptor dba,
                                                    Collection<Long> dbIds,
                                                    String fileName) throws Exception, IOException {
        Set<String> interactions = new HashSet<String>();
        Set<String> rxtFIs = new HashSet<String>();
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        Map<String, Set<Long>> fiToRxtIDs = new HashMap<String, Set<Long>>();
        for (Long dbId : dbIds) {
            GKInstance reaction = dba.fetchInstance(dbId);
            interactors.clear();
            rxtFIs.clear();
            generateFIsForSingleReaction(interactors, rxtFIs, reaction, false, false);
            if (rxtFIs.size() == 0)
                continue;
            interactions.addAll(rxtFIs);
            for (String fi : rxtFIs) {
                InteractionUtilities.addElementToSet(fiToRxtIDs,
                                                     fi, 
                                                     dbId);
            }
        }
        
//        // Feature handler
//        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
//        List<String> featureList = featureHandler.getFeatureList();
//        Map<String, Value> fiToValue = featureHandler.convertPairsToValues(interactions, true);
//        
//        Map<String, String> idToGene = getUniProtToGeneMap(dba);
//        
////        for (String fi : interactions)
////            System.out.println(fi.replace('\t', ' '));
//        FileUtility fu = new FileUtility();
//        fu.setOutput(fileName);
//        fu.printLine("UniProt1\tUniProt2\tGene1\tGene2\tPositiveFeature\tHumanPPI\tMousePPI\tFlyPPI\tWormPPI\tYeastPPI\tDomainInt\tReactions");
//        for (String fi : interactions) {
//            Value value = fiToValue.get(fi);
//            Boolean posFeature = value.humanInteraction |
//                                 value.mousePPI |
//                                 value.dmePPI |
//                                 value.celPPI |
//                                 value.scePPI |
//                                 value.pfamDomainInt;
//            int index = fi.indexOf("\t");
//            String gene1 = idToGene.get(fi.substring(0, index));
//            String gene2 = idToGene.get(fi.substring(index + 1));
//            String rxtIds = fiToRxtIDs.get(fi).toString();
//            rxtIds = rxtIds.substring(1, rxtIds.length() - 1);
//            fu.printLine(fi + "\t" + 
//                               gene1 + "\t" + 
//                               gene2 + "\t" +
//                               posFeature + "\t" + 
//                               value.humanInteraction + "\t" + 
//                               value.mousePPI + "\t" + 
//                               value.dmePPI + "\t" + 
//                               value.celPPI + "\t" + 
//                               value.scePPI + "\t" +
//                               value.pfamDomainInt + "\t" + 
//                               rxtIds);
//        }
//        fu.close();
    }
    
    private void generateFIsForSingleReaction(Set<GKInstance> interactors,
                                              Set<String> interactions,
                                              GKInstance rxn,
                                              boolean includeFIsInComplex,
                                              boolean useGeneName) throws Exception {
        extractInteractorsFromReaction(rxn, interactors);
        generateInteractions(interactors, interactions, rxn, useGeneName);
        if (!includeFIsInComplex)
            return;
        // Collect FIs from complexes involved in reactions
        Set<GKInstance> rxnComplexes = grepComplexes(interactors);
        for (GKInstance complex : rxnComplexes) {
            Set<GKInstance> complexInteractors = new HashSet<GKInstance>();
            grepComplexComponents(complex, complexInteractors);
            generateInteractions(complexInteractors, 
                                 interactions,
                                 complex,
                                 useGeneName);
        }
    }
    
    public Set<String> extractInteractionSet() throws Exception {
        Collection reactions = prepareReactions();
        Collection complexes = prepareComplexes();
        return extractInteractionSet(reactions, 
                                     complexes);
    }

    private Set<String> extractInteractionSet(Collection reactions,
                                              Collection complexes) throws Exception {
        Set<String> interactions = new HashSet<String>();
        GKInstance rxn = null;
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        long time1 = System.currentTimeMillis();
        for (Iterator it = reactions.iterator(); it.hasNext();) {
            rxn = (GKInstance) it.next();
            //System.out.println("Reaction: " + c++);
            extractInteractorsFromReaction(rxn, interactors);
            generateInteractions(interactors, interactions, rxn);
            interactors.clear();
        }
        System.out.println("Total interactions from reactions: " + interactions.size());
        if (!excludeComplex) {
            GKInstance complex = null;
            for (Iterator it = complexes.iterator(); it.hasNext();) {
                complex = (GKInstance) it.next();
                //System.out.println("Complex: " + c++ + " " + complex.getDBID());
                interactors.clear();
                grepComplexComponents(complex, interactors);
                // No need
                //if (interactors.size() > 10)
                //    continue; // cutoff set manually
                generateInteractions(interactors, interactions, complex);
            }
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Time for looping: " + (time2 - time1));
        System.out.println("Total interactions from Reactome: " + interactions.size());
        return interactions;
    }
    
    @SuppressWarnings("rawtypes")
    public void grepComplexComponents(GKInstance complex, Set<GKInstance> interactors) throws Exception {
        Set<GKInstance> current = new HashSet<GKInstance>();
        current.add(complex);
        Set<GKInstance> next = new HashSet<GKInstance>();
        while (current.size() > 0) {
            for (Iterator it = current.iterator(); it.hasNext();) {
                GKInstance tmp = (GKInstance) it.next();
                List components = tmp.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
                if (components == null || components.size() == 0)
                    continue;
                for (Iterator it1 = components.iterator(); it1.hasNext();) {
                    GKInstance tmp1 = (GKInstance) it1.next();
                    if (tmp1.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence))
                        interactors.add(tmp1);
                    else if (tmp1.getSchemClass().isa(ReactomeJavaConstants.EntitySet))
                        interactors.add(tmp1);
                    else if (tmp1.getSchemClass().isa(ReactomeJavaConstants.Complex))
                        next.add(tmp1);
                }
            }
            current.clear();
            current.addAll(next);
            next.clear();
        }
    }
    
    protected void generateInteractions(Set<GKInstance> interactors, 
                                        Set<String> interactions,
                                        GKInstance source) throws Exception {
        generateInteractions(interactors, interactions, source, false);
    }
    
    private void generateInteractions(Set<GKInstance> interactors, 
                                      Set<String> interactions,
                                      GKInstance source,
                                      boolean useGeneNames) throws Exception {
        List<GKInstance> list = new ArrayList<GKInstance>(interactors);
        int size = list.size();
        for (int i = 0; i < size - 1; i++) {
            GKInstance interactor1 = list.get(i);
            for (int j = i + 1; j < size; j++) {
                GKInstance interactor2 = list.get(j);
                generateInteractions(interactor1, 
                                     interactor2,
                                     interactions,
                                     source,
                                     useGeneNames);
            }
        }
    }
    
    /**
     * Generate interactions from the passed interactions into interactions as strings.
     * @param interactors
     * @param interactions
     * @param source
     * @throws Exception
     */
    public void generateInteractionsWithDBNames(Set<GKInstance> interactors,
                                                Set<String> interactions,
                                                GKInstance source) throws Exception {
        List<GKInstance> list = new ArrayList<GKInstance>(interactors);
        int size = list.size();
        for (int i = 0; i < size - 1; i++) {
            GKInstance interactor1 = list.get(i);
            for (int j = i + 1; j < size; j++) {
                GKInstance interactor2 = list.get(j);
                generateInteractionsWithDBNames(interactor1, interactor2, interactions, source);
            }
        }        
    }
    
    private void generateInteractionsWithDBNames(GKInstance interactor1,
                                                 GKInstance interactor2,
                                                 Set<String> interactions,
                                                 GKInstance source) throws Exception {
        Set<GKInstance> refPepSeqs1 = grepRefPepSeqs(interactor1);
        if (refPepSeqs1.size() == 0)
            return;
        Set<GKInstance> refPepSeqs2 = grepRefPepSeqs(interactor2);
        if (refPepSeqs1.size() == 0)
            return;
        // Permutate members in these two sets
        int comp = 0;
        String pair = null;
        for (GKInstance ref1 : refPepSeqs1) {
            String uid1 = (String) ref1.getAttributeValue(ReactomeJavaConstants.identifier);
            if (uid1 == null)
                continue;
            String dbName1 = null;
            GKInstance db1 = (GKInstance) ref1.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
            if (db1 != null)
                dbName1 = db1.getDisplayName();
            else
                dbName1 = "unknown";
            for (GKInstance ref2 : refPepSeqs2) {
                String uid2 = (String) ref2.getAttributeValue(ReactomeJavaConstants.identifier);
                if (uid2 == null)
                    continue;
                String dbName2 = null;
                GKInstance db2 = (GKInstance) ref2.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
                if (db2 != null)
                    dbName2 = db2.getDisplayName();
                else
                    dbName2 = "unknown";
                comp = uid1.compareTo(uid2);
                if (comp < 0)
                    pair = dbName1 + ":" + uid1 + "\t" + dbName2 + ":" + uid2;
                else if (comp > 0)
                    pair = dbName2 + ":" + uid2 + "\t" + dbName1 + ":" + uid1;
                if (pair != null) {
                    interactions.add(pair); //exclude self interaction
                }
            }
        }
    }
    
    private void generateInteractions(GKInstance interactor1, 
                                      GKInstance interactor2,
                                      Set<String> interactions,
                                      GKInstance source,
                                      boolean useGeneNames) throws Exception {
        if (excludeComplex) {
            if (interactor1.getSchemClass().isa(ReactomeJavaConstants.Complex) ||
                interactor2.getSchemClass().isa(ReactomeJavaConstants.Complex))
                return;
        }
        Set<GKInstance> refPepSeqs1 = grepRefPepSeqs(interactor1);
        Set<GKInstance> refPepSeqs2 = grepRefPepSeqs(interactor2);
        generateFIs(refPepSeqs1, 
                    refPepSeqs2,
                    interactions,
                    useGeneNames);
    }

    protected void generateFIs(Set<GKInstance> refPepSeqs1, 
                               Set<GKInstance> refPepSeqs2, 
                               Set<String> interactions,
                               boolean useGeneNames) throws InvalidAttributeException, Exception {
        if (refPepSeqs1.size() == 0 || refPepSeqs2.size() == 0)
            return;
        // Permutate members in these two sets
        int comp = 0;
        String pair = null;
        for (GKInstance ref1 : refPepSeqs1) {
            String uid1 = null;
            if (useGeneNames && ref1.getSchemClass().isValidAttribute(ReactomeJavaConstants.geneName))
                uid1 = (String) ref1.getAttributeValue(ReactomeJavaConstants.geneName);
            else
                uid1 = (String) ref1.getAttributeValue(ReactomeJavaConstants.identifier);
            if (uid1 == null)
                continue;
            for (GKInstance ref2 : refPepSeqs2) {
                String uid2 = null;
                if (useGeneNames && ref2.getSchemClass().isValidAttribute(ReactomeJavaConstants.geneName))
                    uid2 = (String) ref2.getAttributeValue(ReactomeJavaConstants.geneName);
                else
                    uid2 = (String) ref2.getAttributeValue(ReactomeJavaConstants.identifier);
                if (uid2 == null)
                    continue;
                comp = uid1.compareTo(uid2);
                if (comp < 0)
                    pair = uid1 + "\t" + uid2;
                else if (comp > 0)
                    pair = uid2 + "\t" + uid1;
                if (pair != null) {
                    interactions.add(pair); //exclude self interaction
                    // Used for debugging
                    //if (pair.equals("O95405 P01270"))
                    //    System.out.println(pair + " < " + source.getDBID());
                }
            }
        }
    }
    
    protected Set<GKInstance> grepRefPepSeqs(GKInstance interactor) throws Exception {
        return topicHelper.grepRefPepSeqs(interactor);
    }
    
    protected Collection prepareComplexes() throws Exception {
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        GKInstance homosapiens = dba.fetchInstance(48887L);
        Collection complexes = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Complex,
                                                            ReactomeJavaConstants.species,
                                                            "=",
                                                            homosapiens);
        // Change made on January 31, 2017. Previously all human complexes are used,
        // which should be regarded as a bug. Do the following for filtering.
        for (Iterator it = complexes.iterator(); it.hasNext();) {
            GKInstance complex = (GKInstance) it.next();
            GKInstance dataSource = (GKInstance) complex.getAttributeValue(ReactomeJavaConstants.dataSource);
            if (dataSource != null)
                it.remove();
        }
        SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.Complex);
        SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.hasComponent);
        dba.loadInstanceAttributeValues(complexes, att);
        return complexes;
    }
    
    protected Collection prepareReactions() throws Exception {
        // Load necessary attributes
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        // Load all reactions for analyzed
        GKInstance homosapiens = dba.fetchInstance(48887L);
        Collection reactions = null;
        SchemaClass reactionCls = null;
        // Adjust for new schema
        if (dba.getSchema().isValidClass(ReactomeJavaConstants.ReactionlikeEvent)) {
            reactions = dba.fetchInstanceByAttribute(ReactomeJavaConstants.ReactionlikeEvent,
                                                    ReactomeJavaConstants.species,
                                                    "=",
                                                    homosapiens);
            reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.ReactionlikeEvent);
        }
        else {
            reactions = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Reaction,
                                                     ReactomeJavaConstants.species,
                                                     "=",
                                                     homosapiens);
            reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.Reaction);
        }
        // Need a little bit filtering for Reactome reactions only
        for (Iterator it = reactions.iterator(); it.hasNext();) {
            GKInstance rxt = (GKInstance) it.next();
            GKInstance dataSource = (GKInstance) rxt.getAttributeValue(ReactomeJavaConstants.dataSource);
            if (dataSource != null)
                it.remove();
        }
        Collection cas = dba.fetchInstancesByClass(ReactomeJavaConstants.CatalystActivity);
        Collection regulations = dba.fetchInstancesByClass(ReactomeJavaConstants.Regulation);
        Collection entities = dba.fetchInstanceByAttribute(ReactomeJavaConstants.EntityWithAccessionedSequence,
                                                           ReactomeJavaConstants.species,
                                                           "=",
                                                           homosapiens);
        SchemaAttribute att = reactionCls.getAttribute(ReactomeJavaConstants.input);
        dba.loadInstanceAttributeValues(reactions, att);
        att = reactionCls.getAttribute(ReactomeJavaConstants.output);
        dba.loadInstanceAttributeValues(reactions, att);
        att = reactionCls.getAttribute(ReactomeJavaConstants.catalystActivity);
        dba.loadInstanceAttributeValues(reactions, att);
        reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.CatalystActivity);
        att = reactionCls.getAttribute(ReactomeJavaConstants.physicalEntity);
        dba.loadInstanceAttributeValues(cas, att);
        reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.Regulation);
        att = reactionCls.getAttribute(ReactomeJavaConstants.regulatedEntity);
        dba.loadInstanceAttributeValues(regulations, att);
        att = reactionCls.getAttribute(ReactomeJavaConstants.regulator);
        dba.loadInstanceAttributeValues(regulations, att);
        reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.EntityWithAccessionedSequence);
        att = reactionCls.getAttribute(ReactomeJavaConstants.referenceEntity);
        dba.loadInstanceAttributeValues(entities, att);
        return reactions;
    }
    
    private Map<String, String> getUniProtToGeneMap(MySQLAdaptor dba) throws Exception {
        Collection<GKInstance> refSeqs = dba.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceGeneProduct,
                                                                      ReactomeJavaConstants.species,
                                                                      "=",
                                                                      48887L);
        dba.loadInstanceAttributeValues(refSeqs, new String[]{
                ReactomeJavaConstants.identifier,
                ReactomeJavaConstants.geneName
        });
        Map<String, String> idToGene = new HashMap<String, String>();
        for (GKInstance inst : refSeqs) {
            String id = (String) inst.getAttributeValue(ReactomeJavaConstants.identifier);
            String gene = (String) inst.getAttributeValue(ReactomeJavaConstants.geneName);
            idToGene.put(id, gene);
        }
        return idToGene;
    }
}
