package org.reactome.cancer;

import org.reactome.cancer.driver.CancerDriverReactomeAnalyzer;
import org.reactome.r3.ReactomeAnalyzer;

import java.util.*;

public class ReactomeMechismoDataMap {
    private Set<FI> fis;
    private Set<Reaction> reactions;
    private Map<Long, Reaction> reactionIDToReaction;
    private Map<String, FI> alphaUniprotsToFI;
    private Map<Patient, Set<FI>> patientsToFIs;
    private Map<FI, Set<Patient>> fisToPatients;
    private Map<Patient, Map<FI, Set<Mutation>>> patientToFIsToMutations;
    private Map<Patient, Set<Reaction>> patientToReactions;
    private Map<Reaction, Set<FI>> reactionToFIs;
    private Map<Reaction, Set<Patient>> reactionToPatients;
    private CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer;
    private ReactomeReactionGraphLoader reactomeReactionGraphLoader;

    public ReactomeMechismoDataMap(CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer,
                                   MechismoOutputLoader mechismoOutputLoader,
                                   ReactomeReactionGraphLoader reactomeReactionGraphLoader) throws Exception {
        this.cancerDriverReactomeAnalyzer = cancerDriverReactomeAnalyzer;
        this.reactomeReactionGraphLoader = reactomeReactionGraphLoader;
        this.patientsToFIs = mechismoOutputLoader.ExtractPatientsToFIs();
        this.fisToPatients = mechismoOutputLoader.ExtractFIs2Samples();
        this.fis = mechismoOutputLoader.ExtractMechismoFIs();
        BuildReactions();
        BuildSortedUniprotsToFIs();
        BuildReactionToFIs();
        BuildReactionToPatients();
        this.patientToFIsToMutations = mechismoOutputLoader.ExtractSamples2FIs2Muts();
    }

    private Set<Mutation> getReactionPatientMutations(Reaction reaction, Patient patient){
       Set<Mutation> mutations = new HashSet<>();
       if(this.reactionToFIs.keySet().contains(reaction)) {
           Set<FI> reactionFIs = new HashSet<>(this.reactionToFIs.get(reaction));
           reactionFIs.retainAll(this.patientsToFIs.get(patient));
           for (FI fi : reactionFIs) {
               if (this.patientToFIsToMutations.containsKey(patient)) {
                   mutations.addAll(this.patientToFIsToMutations.get(patient).get(fi));
               }
           }
       }
       return mutations;
    }

    public Set<Mutation> getReactionPatientsMutations(Reaction reaction, Set<Patient> patients){
        Set<Mutation> mutations = new HashSet<>();
        for(Patient patient : patients){
            mutations.addAll(getReactionPatientMutations(reaction,patient));
        }
        return mutations;
    }

    public Set<FI> getReactionPatientsFIs(Reaction reaction, Set<Patient> patients) {
        Set<FI> reactionPatientsFIs = new HashSet<>();
        for (Patient patient : patients) {
            Set<FI> patientFIs = new HashSet<>(this.patientsToFIs.get(patient));
            patientFIs.retainAll(this.reactionToFIs.get(reaction));
            reactionPatientsFIs.addAll(patientFIs);
        }
        return reactionPatientsFIs;
    }

    private void BuildReactions() throws Exception {
        this.reactions = new HashSet<>();
        this.reactionIDToReaction = new HashMap<>();
        Map<Long,String> reactionIDToName = cancerDriverReactomeAnalyzer.loadReactionLongDBIDToName();
        for (Long reactionID : this.reactomeReactionGraphLoader.getReactionSet()) {
            if(reactionIDToName.containsKey(reactionID)) {
                Reaction reaction = new Reaction(
                        reactionID,
                        reactionIDToName.get(reactionID));
                this.reactions.add(reaction);
                this.reactionIDToReaction.put(reactionID, reaction);
            }
        }
    }

    private void BuildSortedUniprotsToFIs() {
        this.alphaUniprotsToFI = new HashMap<>();
        for (FI fi : this.fis) {
            Gene[] genes = fi.getGenes();
            String uniprot1 = genes[0].getUniprotID();
            String uniprot2 = genes[1].getUniprotID();
            String[] uniprots = new String[]{uniprot1, uniprot2};
            Arrays.sort(uniprots);
            String key = String.format("%s:%s", uniprots[0], uniprots[1]);
            alphaUniprotsToFI.put(key, fi);
        }
    }

    public Set<Patient> getPatients() {
        return this.patientsToFIs.keySet();
    }

    public Reaction getReaction(Long reactionID) {
        return this.reactionIDToReaction.get(reactionID);
    }

    public Set<FI> getFIs(Reaction reaction) {
        return this.reactionToFIs.get(reaction);
    }

    public Set<Patient> getPatients(Reaction reaction) {
        Set<Patient> rxnPatients = this.reactionToPatients.get(reaction);
        return rxnPatients == null
                ? new HashSet<>()
                : rxnPatients;
    }

    public Set<FI> getFIs() {
        return this.fis;
    }

    public Integer getNumPatients() {
        return this.patientsToFIs.keySet().size();
    }

    private void BuildReactionToFIs() throws Exception {
        this.reactionToFIs = new HashMap<>();
        Map<FI, Set<Reaction>> fisToReactions = new HashMap<>();
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        for(Reaction reaction : this.reactions) {
            List<Long> reactionIDList = new ArrayList<>();
            reactionIDList.add(reaction.getReactionID());
            Set<String> reactionFIStrings = new HashSet<>();
            reactomeAnalyzer.generateFIsForReactionsWithFeatures(
                    this.cancerDriverReactomeAnalyzer.getDBA(),
                    reactionIDList,
                    null,
                    reactionFIStrings);

            Set<FI> reactionFIs = new HashSet<>();
            for (String fiString : reactionFIStrings) {
                String[] uniprots = fiString.split("\t");
                Arrays.sort(uniprots);
                String key = String.format("%s:%s", uniprots[0], uniprots[1]);
                reactionFIs.add(this.alphaUniprotsToFI.get(key));
            }

            //reaction FIs now in interactions set
            Set<FI> intersection = new HashSet<>(reactionFIs);
            intersection.retainAll(fis);
            if (intersection.size() > 0) {
                //if Mechismo FIs in reaction interactions
                this.reactionToFIs.put(reaction, intersection);
                for (FI interaction :
                        intersection) {
                    Set<Reaction> interactionRxns;
                    if (fisToReactions.containsKey(interaction)) {
                        interactionRxns = fisToReactions.get(interaction);
                    } else {
                        interactionRxns = new HashSet<>();
                    }
                    interactionRxns.add(reaction);
                    fisToReactions.put(interaction, interactionRxns);
                }
            }
        }
    }

    public Set<Reaction> getSupportedReactions() {
        return this.reactionToFIs.keySet();
    }

    public Set<Reaction> getReactions(Patient patient){
        return this.patientToReactions.get(patient);
    }

    private void BuildReactionToPatients() {
        this.reactionToPatients = new HashMap<>();
        this.patientToReactions = new HashMap<>();
        for (Reaction reaction : this.reactionToFIs.keySet()) {
            Set<FI> rxnFis = this.reactionToFIs.get(reaction);
            for (FI rxnFi : rxnFis) {
                Set<Patient> rxnPatients;
                if (reactionToPatients.containsKey(reaction)) {
                    rxnPatients = new HashSet<>(reactionToPatients.get(reaction));
                } else {
                    rxnPatients = new HashSet<>();
                }
                rxnPatients.addAll(fisToPatients.get(rxnFi));
                reactionToPatients.put(reaction, rxnPatients);
                for(Patient patient : rxnPatients) {
                    Set<Reaction> patientRxns;
                    if (patientToReactions.containsKey(patient)){
                        patientRxns = new HashSet<>(patientToReactions.get(patient));
                    }else{
                        patientRxns = new HashSet<>();
                    }
                    patientRxns.add(reaction);
                    patientToReactions.put(patient,patientRxns);
                }
            }
        }
    }
}
