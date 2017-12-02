package org.reactome.cancer;

import org.gk.util.StringUtils;
import org.reactome.cancer.driver.CancerDriverReactomeAnalyzer;
import org.reactome.r3.ReactomeAnalyzer;
import org.reactome.r3.util.FileUtility;

import java.io.IOException;
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

    public Set<FI> getFIs(Patient patient){
        return this.patientToFIsToMutations.get(patient).keySet();
    }

    private Set<Mutation> getReactionPatientMutations(Reaction reaction, Patient patient) {
        Set<Mutation> mutations = new HashSet<>();
        if (this.reactionToFIs.keySet().contains(reaction)) {
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

    public Set<Mutation> getReactionPatientsMutations(Reaction reaction, Set<Patient> patients) {
        Set<Mutation> mutations = new HashSet<>();
        for (Patient patient : patients) {
            mutations.addAll(getReactionPatientMutations(reaction, patient));
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
        Map<Long, String> reactionIDToName = cancerDriverReactomeAnalyzer.loadReactionLongDBIDToName();
        for (Long reactionID : this.reactomeReactionGraphLoader.getReactionSet()) {
            if (reactionIDToName.containsKey(reactionID)) {
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
        for (Reaction reaction : this.reactions) {
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

    public Set<Reaction> getReactions(Patient patient) {
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
                for (Patient patient : rxnPatients) {
                    Set<Reaction> patientRxns;
                    if (patientToReactions.containsKey(patient)) {
                        patientRxns = new HashSet<>(patientToReactions.get(patient));
                    } else {
                        patientRxns = new HashSet<>();
                    }
                    patientRxns.add(reaction);
                    patientToReactions.put(patient, patientRxns);
                }
            }
        }
    }

    private void FindIntersectionUnionClusters(Set<Set<String>> fiIntersectingSetUnionClusters,
                                               Set<Set<String>> allFIs) {
        Iterator<Set<String>> allFIsItr = allFIs.iterator();
        while (allFIsItr.hasNext()) {
            Set<String> fiSet = allFIsItr.next();
            Set<String> union = new HashSet<>(fiSet);
            Iterator<Set<String>> fiClusterItr = fiIntersectingSetUnionClusters.iterator();
            while (fiClusterItr.hasNext()) {
                Set<String> fiCluster = fiClusterItr.next();
                Set<String> intersection = new HashSet<>(fiSet);
                intersection.retainAll(fiCluster);
                if (!intersection.isEmpty()) {
                    union.addAll(fiCluster);
                    fiClusterItr.remove(); //remove out-dated cluster
                }
            }
            fiIntersectingSetUnionClusters.add(union);
        }
    }

    private void CalculateClusterStats(Set<Set<String>> fiIntersectingSetUnionClusters) {
        Iterator<Set<String>> fiClusterItr = fiIntersectingSetUnionClusters.iterator();
        Double min = (double) fiClusterItr.next().size();
        Double max = min;
        Double sum = max;
        while (fiClusterItr.hasNext()) {
            Double sz = (double) fiClusterItr.next().size();
            min = sz < min ? sz : min;
            max = sz > max ? sz : max;
            sum += sz;
        }
        Double mean = sum / (double) fiIntersectingSetUnionClusters.size();

        System.out.println(String.format("Clusters min = %f, max = %f, mean = %f, sum = %f",
                min,
                max,
                mean,
                sum));
    }

    private void WriteClustersToFile(String outputDir,
                                     String outputFilePrefix,
                                     Set<Set<String>> fiIntersectingSetUnionClusters) {
        //perform pathway enrichment analysis among FI clusters
        FileUtility fileUtility0 = new FileUtility();
        String outFilePath0 = outputDir + outputFilePrefix + "fiClusters.csv";
        try {
            fileUtility0.setOutput(outFilePath0);
            for (Set<String> fis0 : fiIntersectingSetUnionClusters) {
                fileUtility0.printLine(
                        fis0.size() > 1
                                ? StringUtils.join(" ",
                                new ArrayList(fis0)).replace("\t", " ")
                                : Collections.singletonList(fis0).get(0).toString()
                                .replace("[", "")
                                .replace("]", "")
                                .replace("\t", " "));
            }
            fileUtility0.close();
        } catch (IOException ioe) {
            System.out.println(String.format("Couldn't use %s, %s: %s",
                    outFilePath0,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    private void WriteReaction2FIsToFile(String outputDir,
                                         String outputFilePrefix,
                                         Map<Long, Set<String>> reaction2FiSet) {
        //write reaction2FIs to file
        FileUtility fileUtility1 = new FileUtility();
        String outFilePath1 = outputDir + outputFilePrefix + "reactions2FIs.csv";
        try {
            fileUtility1.setOutput(outFilePath1);
            fileUtility1.printLine("RxnId,Num FIs,Mapped FIs");
            for (Map.Entry<Long, Set<String>> pair : reaction2FiSet.entrySet()) {
                Long rxnId = pair.getKey();
                Set<String> mappedFis = pair.getValue();
                fileUtility1.printLine(String.format("%s,%d,%s",
                        rxnId.toString(),
                        mappedFis.size(),
                        mappedFis.size() > 1
                                ? StringUtils.join(" ",
                                new ArrayList(mappedFis)).replace("\t", " ")
                                : Collections.singletonList(mappedFis).get(0).toString()
                                .replace("[", "")
                                .replace("]", "")
                                .replace("\t", " ")));
            }
            fileUtility1.close();
        } catch (IOException ioe) {
            System.out.println(String.format("Couldn't use %s, %s: %s",
                    outFilePath1,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    private void WriteSamples2FIsToFile(String outputDir,
                                        String outputFilePrefix,
                                        Map<String, Set<String>> samples2FIs) {
        //write samples2FIs to file
        FileUtility fileUtility2 = new FileUtility();
        String outFilePath2 = outputDir + outputFilePrefix + "samples2FIs.csv";
        try {
            fileUtility2.setOutput(outFilePath2);
            fileUtility2.printLine("Sample Barcode,Num FIs,Mapped FIs");
            for (Map.Entry<String, Set<String>> pair : samples2FIs.entrySet()) {
                String sampleBarcode = pair.getKey();
                Set<String> mappedFis = pair.getValue();
                fileUtility2.printLine(String.format("%s,%d,%s",
                        sampleBarcode,
                        mappedFis.size(),
                        mappedFis.size() > 1
                                ? StringUtils.join(" ",
                                new ArrayList(mappedFis)).replace("\t", " ")
                                : Collections.singletonList(mappedFis).get(0).toString()
                                .replace("[", "")
                                .replace("]", "")
                                .replace("\t", " ")));
            }
            fileUtility2.close();
        } catch (IOException ioe) {
            System.out.println(String.format("Couldn't use %s, %s: %s",
                    outFilePath2,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    private void WriteFIs2SamplesToFile(String outputDir,
                                        String outputFilePrefix,
                                        Map<String, Set<String>> fis2Samples,
                                        Map<String, Set<String>> samples2FIs) {
        //write fis2Samples to file
        FileUtility fileUtility3 = new FileUtility();
        String outFilePath3 = outputDir + outputFilePrefix + "fis2Samples.csv";
        try {
            fileUtility3.setOutput(outFilePath3);
            fileUtility3.printLine("FI,Num Samples,FI Frequency,Mapped Samples");
            for (Map.Entry<String, Set<String>> pair : fis2Samples.entrySet()) {
                String fi = pair.getKey();
                Set<String> mappedSamples = pair.getValue();
                fileUtility3.printLine(String.format("%s,%d,%f,%s",
                        fi,
                        mappedSamples.size(),
                        (double) mappedSamples.size() /
                                (double) samples2FIs.keySet().size(),
                        mappedSamples.size() > 1
                                ? StringUtils.join(" ",
                                new ArrayList(mappedSamples)).replace(
                                "\t", " ")
                                : Arrays.asList(mappedSamples).get(0).toString()
                                .replace("[", "")
                                .replace("]", "")
                                .replace("\t", " ")));
            }
            fileUtility3.close();
        } catch (IOException ioe) {
            System.out.println(String.format("Couldn't use %s, %s: %s",
                    outFilePath3,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    private void WriteTargetReactionSummariesToFile(String outputDir,
                                                    String outputFilePrefix,
                                                    Set<TargetReactionCandidate> targetReactionSummaries) {
        //write targetReactionSummaries to file
        FileUtility fileUtility = new FileUtility();
        String outFilePath = outputDir + outputFilePrefix + "dnUpReactions.csv";
        try {
            fileUtility.setOutput(outFilePath);
            fileUtility.printLine(TargetReactionCandidate.getHeaderLine());
            for (TargetReactionCandidate targetReactionCandidate : targetReactionSummaries) {
                fileUtility.printLine(targetReactionCandidate.toString());
            }
            fileUtility.close();
        } catch (IOException ioe) {
            System.out.println(String.format("Couldn't use %s, %s: %s",
                    outFilePath,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }

        System.out.println(String.format("Output written to %s",
                outFilePath));
    }

    public void WriteRxn2SamplesToFile(String outputDir,
                                       String outputFilePrefix) {
        //write rxn2Samples to file
        FileUtility fileUtility = new FileUtility();
        String outFilePath = outputDir + outputFilePrefix + "rxn2Samples.csv";
        try {
            fileUtility.setOutput(outFilePath);
            fileUtility.printLine("Rxn,Num Samples,Rxn Frequency,Mapped Samples");
            for (Reaction reaction : getSupportedReactions()) {
                Set<Patient> mappedSamples = getPatients(reaction);
                fileUtility.printLine(String.format("%s,%d,%f,%s",
                        reaction.toString(),
                        mappedSamples.size(),
                        (double) mappedSamples.size() /
                                (double) getSupportedReactions().size(),
                        mappedSamples.size() > 1
                                ? org.gk.util.StringUtils.join(" ",
                                new ArrayList<>(mappedSamples)).replace("\t", " ")
                                : Collections.singletonList(mappedSamples).get(0).toString()
                                .replace("[", "")
                                .replace("]", "")
                                .replace("\t", " ")));
            }
            fileUtility.close();
        } catch (IOException ioe) {
            System.out.println(String.format("Couldn't use %s, %s: %s",
                    outFilePath,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }
}
