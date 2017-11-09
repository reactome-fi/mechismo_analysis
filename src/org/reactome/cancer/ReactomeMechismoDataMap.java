package org.reactome.cancer;

import org.reactome.cancer.driver.CancerDriverReactomeAnalyzer;
import org.reactome.r3.ReactomeAnalyzer;

import java.util.*;

public class ReactomeMechismoDataMap {
    private Set<FI> fis;
    private Map<String,FI> alphaUniprotsToFI;
    private Set<Reaction> reactions;
    private Map<Long,Reaction> reactionIDToReaction;
    private Map<Patient,Set<FI>> samples2FIs;
    private Map<FI,Set<Patient>> fis2Samples;
    private Map<FI,Set<Reaction>> fisToReactions;
    private Map<Reaction,Set<FI>> reactionToFIs;
    private Map<Reaction,Set<Patient>> reactionToPatientMap;
    private Map<Patient, Map<FI, Set<Mutation>>> reactionToFIToMutationMap;
    private CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer;
    private Map<Mutation, Integer> mut2SampleCount;
    private Map<Gene, Integer> gene2SampleCount;
    private MechismoOutputLoader mechismoOutputLoader;
    private ReactomeReactionGraphLoader reactomeReactionGraphLoader;
    public ReactomeMechismoDataMap(CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer,
                                   MechismoOutputLoader mechismoOutputLoader,
                                   ReactomeReactionGraphLoader reactomeReactionGraphLoader) throws Exception {
        this.cancerDriverReactomeAnalyzer = cancerDriverReactomeAnalyzer;
        this.mechismoOutputLoader = mechismoOutputLoader;
        this.reactomeReactionGraphLoader = reactomeReactionGraphLoader;
        this.samples2FIs = this.mechismoOutputLoader.ExtractSamples2FIs();
        this.fis2Samples = this.mechismoOutputLoader.ExtractFIs2Samples();
        this.fis = mechismoOutputLoader.ExtractMechismoFIs();
        this.mut2SampleCount = mechismoOutputLoader.getMut2SampleCount();
        this.gene2SampleCount = mechismoOutputLoader.getGene2SampleCount();
        BuildReactions();
        BuildSortedUniprotsToFIs();
        BuildReactionToFIs();
        BuildReactionToSamples();
        this.reactionToFIToMutationMap = mechismoOutputLoader.ExtractSamples2FIs2Muts();
    }

    private void BuildReactions() throws Exception {
        this.reactions = new HashSet<>();
        this.reactionIDToReaction = new HashMap<>();
        Map<Long,String> reactionLongDBIDToName = cancerDriverReactomeAnalyzer.loadReactionLongDBIDToName();
        for(Long reactionID : reactionLongDBIDToName.keySet()){
            Reaction reaction = new Reaction(reactionID,reactionLongDBIDToName.get(reactionID));
            this.reactions.add(reaction);
            this.reactionIDToReaction.put(reactionID,reaction);
        }
    }

    private Integer getSampleCount(Mutation mutation){
        return this.mut2SampleCount.get(mutation);
    }

    private Integer getSampleCount(Gene gene){
        return this.gene2SampleCount.get(gene);
    }

    private void BuildSortedUniprotsToFIs(){
        this.alphaUniprotsToFI = new HashMap<>();
        for(FI fi : this.fis){
            Gene[] genes = fi.getGenes();
            String uniprot1 = genes[0].getUniprotID();
            String uniprot2 = genes[1].getUniprotID();
            String[] uniprots = new String[]{uniprot1,uniprot2};
            Arrays.sort(uniprots);
            String key = String.format("%s:%s",uniprots[0],uniprots[1]);
            alphaUniprotsToFI.put(key,fi);
        }
    }

    public Reaction getReaction(Long reactionID){
        return this.reactionIDToReaction.get(reactionID);
    }

    public Set<FI> getFIs(Patient patient){
        return this.samples2FIs.get(patient);
    }

    public Set<FI> getFIs(Reaction reaction){
        return this.reactionToFIs.get(reaction);
    }

    public Set<Patient> getPatients(FI fi){
        return this.fis2Samples.get(fi);
    }

    public Set<Reaction> getReactions(FI fi){
        return this.fisToReactions.get(fi);
    }

    public Set<Patient> getPatients(Reaction reaction){
        return this.reactionToPatientMap.get(reaction);
    }

    public Set<Mutation> getMutations(Patient patient,FI fi){
        return this.reactionToFIToMutationMap.get(patient).get(fi);
    }

    public Set<FI> getFIs(){
        return this.fis;
    }

    public Integer getNumPatients(){
        return this.samples2FIs.keySet().size();
    }

    private void BuildReactionToFIs() throws Exception {
        this.reactionToFIs = new HashMap<>();
        this.fisToReactions = new HashMap<>();
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        Iterator<Reaction> rxnItr = this.reactions.iterator();
        while (rxnItr.hasNext()) {
            Reaction reaction = rxnItr.next();
            List<Long> reactionIDList = new ArrayList<>();
            reactionIDList.add(reaction.getReactionID());
            Set<String> reactionFIStrings = new HashSet<>();
            reactomeAnalyzer.generateFIsForReactionsWithFeatures(
                    this.cancerDriverReactomeAnalyzer.getDBA(),
                    reactionIDList,
                    null,
                    reactionFIStrings);

            Set<FI> reactionFIs = new HashSet<>();
            for(String fiString : reactionFIStrings){
                String[] uniprots = fiString.split("\t");
                Arrays.sort(uniprots);
                String key = String.format("%s:%s",uniprots[0],uniprots[1]);
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

    public Set<Reaction> getSupportedReactions(){
        return this.reactionToFIs.keySet();
    }

    private void BuildReactionToSamples() {
        this.reactionToPatientMap = new HashMap<>();
        for (Reaction reaction : this.reactionToFIs.keySet()) {
            Set<FI> rxnFis = this.reactionToFIs.get(reaction);
            for (FI rxnFi : rxnFis) {
                Set<Patient> rxnSamples;
                if (reactionToPatientMap.containsKey(reaction)) {
                    rxnSamples = new HashSet<>(reactionToPatientMap.get(reaction));
                } else {
                    rxnSamples = new HashSet<>();
                }
                rxnSamples.addAll(fis2Samples.get(rxnFi));
                reactionToPatientMap.put(reaction, rxnSamples);
            }
        }
    }
}
