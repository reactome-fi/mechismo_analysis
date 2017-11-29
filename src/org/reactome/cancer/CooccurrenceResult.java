package org.reactome.cancer;

import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.MathUtilities;

import java.io.IOException;
import java.util.*;

public class CooccurrenceResult {
    private final Double MAX_CLUSTER_EMPIRICAL_P_VALUE = 0.05;
    private List<Reaction> targetRxns;
    private List<Set<Reaction>> cooccurringUpstreamRxns;
    private List<Set<FI>> cooccurringUpstreamRxnFIs;
    private List<Set<Patient>> samplesWTargetRxnMutations;
    private List<Integer> numSamplesW0MutatedUpstreamRxns;
    private List<Set<Patient>> samplesW1MutatedUpstreamRxn;
    private List<Set<Patient>> samplesW2MutatedUpstreamRxns;
    private List<Set<Patient>> samplesW3plusMutatedUpstreamRxns;
    private List<Set<Mutation>> superIndirectMutations;
    private List<Set<Mutation>> indirectMutations;
    private List<Set<Mutation>> superDirectMutations;
    private List<Set<Mutation>> directMutations;
    private List<Double> pValues;
    private List<Set<Gene>> superIndirectMutatedGenes;
    private List<Set<Gene>> indirectMutatedGenes;
    private List<Set<Gene>> superDirectMutatedGenes;
    private List<Set<Gene>> directMutatedGenes;
    private Map<Double, Double> pValue2BHAdjustedPValueMap;
    private Map<Double, Double> pValue2EmpiricalPValueMap;
    private Map<Integer, List<List<Patient>>> patientsWith123TMutations;
    private Map<Integer, Set<Patient>> hashCodeToAgPatients;

    public CooccurrenceResult(
            List<Reaction> targetRxns,
            List<Set<Reaction>> cooccurringUpstreamRxns,
            List<Set<FI>> cooccurringUpstreamRxnFIs,
            List<Set<Patient>> samplesWTargetRxnMutations,
            List<Integer> numSamplesW0MutatedUpstreamRxns,
            List<Set<Patient>> samplesW1MutatedUpstreamRxn,
            List<Set<Patient>> samplesW2MutatedUpstreamRxns,
            List<Set<Patient>> samplesW3plusMutatedUpstreamRxns,
            List<Set<Mutation>> superIndirectMutations,
            List<Set<Mutation>> indirectMutations,
            List<Set<Mutation>> superDirectMutations,
            List<Set<Mutation>> directMutations,
            List<Double> pValues
    ) {
        this.targetRxns = targetRxns;
        this.cooccurringUpstreamRxns = cooccurringUpstreamRxns;
        this.cooccurringUpstreamRxnFIs = cooccurringUpstreamRxnFIs;
        this.samplesWTargetRxnMutations = samplesWTargetRxnMutations;
        this.numSamplesW0MutatedUpstreamRxns = numSamplesW0MutatedUpstreamRxns;
        this.samplesW1MutatedUpstreamRxn = samplesW1MutatedUpstreamRxn;
        this.samplesW2MutatedUpstreamRxns = samplesW2MutatedUpstreamRxns;
        this.samplesW3plusMutatedUpstreamRxns = samplesW3plusMutatedUpstreamRxns;
        this.superIndirectMutations = superIndirectMutations;
        this.indirectMutations = indirectMutations;
        this.superDirectMutations = superDirectMutations;
        this.directMutations = directMutations;
        this.pValues = pValues;

        this.hashCodeToAgPatients = null;
        this.superIndirectMutatedGenes = null;
        this.indirectMutatedGenes = null;
        this.superDirectMutatedGenes = null;
        this.directMutatedGenes = null;
        this.pValue2BHAdjustedPValueMap = null;
        this.pValue2EmpiricalPValueMap = null;
        this.patientsWith123TMutations = null;
    }

    private void FillMutatedGeneSets() {

        this.superIndirectMutatedGenes = new ArrayList<>();
        this.indirectMutatedGenes = new ArrayList<>();
        this.superDirectMutatedGenes = new ArrayList<>();
        this.directMutatedGenes = new ArrayList<>();

        for (int i = 0; i < this.pValues.size(); i++) {
            Set<Gene> superIndirectTargetMutatatedGenes = new HashSet<>();
            Set<Gene> indirectTargetMutatedGenes = new HashSet<>();
            Set<Gene> superDirectTargetMutatedGenes = new HashSet<>();
            Set<Gene> directTargetMutatedGenes = new HashSet<>();
            for (Mutation mutation : this.superIndirectMutations.get(i)) {
                superIndirectTargetMutatatedGenes.add(mutation.getGene());
            }
            for (Mutation mutation : this.indirectMutations.get(i)) {
                indirectTargetMutatedGenes.add(mutation.getGene());
            }
            for (Mutation mutation : this.superDirectMutations.get(i)) {
                superDirectTargetMutatedGenes.add(mutation.getGene());
            }
            for (Mutation mutation : this.directMutations.get(i)) {
                directTargetMutatedGenes.add(mutation.getGene());
            }
            this.superIndirectMutatedGenes.add(superIndirectTargetMutatatedGenes);
            this.indirectMutatedGenes.add(indirectTargetMutatedGenes);
            this.superDirectMutatedGenes.add(superDirectTargetMutatedGenes);
            this.directMutatedGenes.add(directTargetMutatedGenes);
        }
    }

    private List<Set<Gene>> getSuperIndirectMutatedGenes() {
        if (this.superIndirectMutatedGenes == null) {
            FillMutatedGeneSets();
        }
        return this.superIndirectMutatedGenes;
    }

    private List<Set<Gene>> getIndirectMutatedGenes() {
        if (this.indirectMutatedGenes == null) {
            FillMutatedGeneSets();
        }
        return this.indirectMutatedGenes;
    }

    private List<Set<Gene>> getSuperDirectMutatedGenes() {
        if (this.superDirectMutatedGenes == null) {
            FillMutatedGeneSets();
        }
        return this.superDirectMutatedGenes;
    }

    private List<Set<Gene>> getDirectMutatedGenes() {
        if (this.directMutatedGenes == null) {
            FillMutatedGeneSets();
        }
        return this.directMutatedGenes;
    }

    public void CalculateBHAdjustedPValues() {
        //TODO: check map is not null
        //The FDR calculation has the side effect of sorting the passed list...
        List<Double> pValuesSorted = new ArrayList<>(this.pValues);

        List<Double> BHAdjustedPValues =
                MathUtilities.calculateFDRWithBenjaminiHochberg(pValuesSorted);

        this.pValue2BHAdjustedPValueMap = new HashMap<>();
        for (int i = 0; i < pValues.size(); i++) {
            this.pValue2BHAdjustedPValueMap.put(pValuesSorted.get(i), BHAdjustedPValues.get(i));
        }
    }

    //GSEA FDR Calculation
    public void CalculateEmpiricalPValues(List<CooccurrenceResult> rewiredNetworkResults) {
        //TODO: check map is not null

        List<Double> rewiredNetworkPvalues = new ArrayList<>();

        double p5 = this.pValues.size() * 1.5d;
        double m5 = this.pValues.size() * 0.5d;
        //loop over rewirings
        for (CooccurrenceResult rewiredNetworkResult : rewiredNetworkResults) {

            if (rewiredNetworkResult.pValues.size() > p5 ||
                    rewiredNetworkResult.pValues.size() < m5) {
                System.out.println(String.format(
                        "INFO: permuted network had %d target reactions, real network had %d",
                        rewiredNetworkResult.pValues.size(),
                        this.pValues.size()));
            }

            rewiredNetworkPvalues.addAll(rewiredNetworkResult.pValues);
        }

        this.pValue2EmpiricalPValueMap = new HashMap<>();

        //loop over pValues
        for (Double pValue : this.pValues) {
            this.pValue2EmpiricalPValueMap.put(pValue,
                    CalculateEmpiricalPValue(
                            CountValuesLTEThresh(rewiredNetworkPvalues, pValue),
                            CountValuesLTEThresh(this.pValues, pValue),
                            rewiredNetworkPvalues.size(),
                            this.pValues.size()));
        }
    }

    private double CalculateEmpiricalPValue(int randomLessThan,
                                            int realLessThan,
                                            int randomTotal,
                                            int realTotal) {

        double fdr = ((double) randomLessThan /
                (double) randomTotal) /
                ((double) realLessThan /
                        (double) realTotal);
        return MathUtilities.boundDouble01(fdr);
    }

    private int CountValuesLTEThresh(List<Double> vals, Double threshold) {
        int count = 0;
        for (Double val : vals) {
            count = val <= threshold ?
                    count + 1 :
                    count;
        }
        return count;
    }

    public void MagicallyShrinkMemoryFootprint() {
        //leave only this.pValues
        if (targetRxns != null) targetRxns.clear();
        if (cooccurringUpstreamRxns != null) cooccurringUpstreamRxns.clear();
        if (cooccurringUpstreamRxnFIs != null) cooccurringUpstreamRxnFIs.clear();
        if (samplesWTargetRxnMutations != null) samplesWTargetRxnMutations.clear();
        if (numSamplesW0MutatedUpstreamRxns != null) numSamplesW0MutatedUpstreamRxns.clear();
        if (samplesW1MutatedUpstreamRxn != null) samplesW1MutatedUpstreamRxn.clear();
        if (samplesW2MutatedUpstreamRxns != null) samplesW2MutatedUpstreamRxns.clear();
        if (samplesW3plusMutatedUpstreamRxns != null) samplesW3plusMutatedUpstreamRxns.clear();
        if (superIndirectMutations != null) superIndirectMutations.clear();
        if (indirectMutations != null) indirectMutations.clear();
        if (superDirectMutations != null) superDirectMutations.clear();
        if (directMutations != null) directMutations.clear();
        if (superIndirectMutatedGenes != null) superIndirectMutatedGenes.clear();
        if (indirectMutatedGenes != null) indirectMutatedGenes.clear();
        if (superDirectMutatedGenes != null) superDirectMutatedGenes.clear();
        if (directMutatedGenes != null) directMutatedGenes.clear();
        if (pValue2BHAdjustedPValueMap != null) pValue2BHAdjustedPValueMap.clear();
        if (pValue2EmpiricalPValueMap != null) pValue2EmpiricalPValueMap.clear();
        System.gc();
    }

    private Integer calculateReactionDistanceBetweenPatients(Patient patient1,
                                                             Patient patient2,
                                                             ReactomeMechismoDataMap reactomeMechismoDataMap) {

        Set<Reaction> patient1Rxns = reactomeMechismoDataMap.getReactions(patient1);
        Set<Reaction> patient2Rxns = reactomeMechismoDataMap.getReactions(patient2);

        if (patient1Rxns != null && patient2Rxns != null) {
            Set<Reaction> patient1RxnsCpy = new HashSet<>(patient1Rxns);
            Set<Reaction> patient2RxnsCpy = new HashSet<>(patient2Rxns);
            patient1RxnsCpy.removeAll(patient2Rxns);
            Integer patient1OnlyRxns = patient1RxnsCpy.size();
            patient2RxnsCpy.removeAll(patient1Rxns);
            Integer patient2OnlyRxns = patient2RxnsCpy.size();
            return patient1OnlyRxns + patient2OnlyRxns;
        }
        return reactomeMechismoDataMap.getSupportedReactions().size();
    }

    private Double calculateMeanReactionDistanceToReferencePatients(Patient queryPatient,
                                                                    Collection<Patient> referencePatients,
                                                                    ReactomeMechismoDataMap reactomeMechismoDataMap) {
        Integer totalDistance = 0;
        for (Patient referencePatient : referencePatients) {
            totalDistance += calculateReactionDistanceBetweenPatients(queryPatient, referencePatient, reactomeMechismoDataMap);
        }
        return (double) totalDistance / (double) referencePatients.size();

    }

    private Set<Patient> aggregateAllPatients(List<Set<Patient>> listOfSets) {
        if (this.hashCodeToAgPatients == null) {
            this.hashCodeToAgPatients = new HashMap<>();
        }
        if (!this.hashCodeToAgPatients.containsKey(listOfSets.hashCode())) {
            Set<Patient> ag = new HashSet<>();
            for (int i = 0; i < targetRxns.size(); i++) {
                ag.addAll(listOfSets.get(i));
            }
            this.hashCodeToAgPatients.put(listOfSets.hashCode(), ag);
        }
        return this.hashCodeToAgPatients.get(listOfSets.hashCode());
    }

    private Double getMeanDistanceToAllClusteredPatients(Patient queryPatient,
                                                         ReactomeMechismoDataMap reactomeMechismoDataMap) {

        Set<Patient> ag1 = aggregateAllPatients(samplesW1MutatedUpstreamRxn);
        Set<Patient> ag2 = aggregateAllPatients(samplesW2MutatedUpstreamRxns);
        Set<Patient> ag3 = aggregateAllPatients(samplesW3plusMutatedUpstreamRxns);
        Set<Patient> agT = aggregateAllPatients(samplesWTargetRxnMutations);

        Integer totalSize = ag1.size() +
                ag2.size() +
                ag3.size() +
                agT.size();

        Double ag1d = calculateMeanReactionDistanceToReferencePatients(
                queryPatient,
                ag1,
                reactomeMechismoDataMap) * (double) ag1.size();
        Double ag2d = calculateMeanReactionDistanceToReferencePatients(
                queryPatient,
                ag2,
                reactomeMechismoDataMap) * (double) ag2.size();
        Double ag3d = calculateMeanReactionDistanceToReferencePatients(
                queryPatient,
                ag3,
                reactomeMechismoDataMap) * (double) ag3.size();
        Double agTd = calculateMeanReactionDistanceToReferencePatients(
                queryPatient,
                agT,
                reactomeMechismoDataMap) * (double) agT.size();

        ag1d = ag1d.isNaN()
                ? 0.0d
                : ag1d;
        ag2d = ag2d.isNaN()
                ? 0.0d
                : ag2d;
        ag3d = ag3d.isNaN()
                ? 0.0d
                : ag3d;
        agTd = agTd.isNaN()
                ? 0.0d
                : agTd;

        Double ret = (ag1d + ag2d + ag3d + agTd) / (double) totalSize;

        if(ret.isNaN()){
            int debug = 1;
        }

        return ret;
    }

    private Double getMeanDistanceToTargetReactionGroupAllClusteredPatients(Patient queryPatient,
                                                                            ReactomeMechismoDataMap reactomeMechismoDataMap,
                                                                            List<List<Patient>> referencePatients) {

        Collection<Patient> ag1 = referencePatients.get(0);
        Collection<Patient> ag2 = referencePatients.get(1);
        Collection<Patient> ag3 = referencePatients.get(2);
        Collection<Patient> agT = referencePatients.get(3);

        Integer totalSize = ag1.size() +
                ag2.size() +
                ag3.size() +
                agT.size();

        Double ag1d = calculateMeanReactionDistanceToReferencePatients(
                queryPatient,
                ag1,
                reactomeMechismoDataMap) * (double) ag1.size();
        Double ag2d = calculateMeanReactionDistanceToReferencePatients(
                queryPatient,
                ag2,
                reactomeMechismoDataMap) * (double) ag2.size();
        Double ag3d = calculateMeanReactionDistanceToReferencePatients(
                queryPatient,
                ag3,
                reactomeMechismoDataMap) * (double) ag3.size();
        Double agTd = calculateMeanReactionDistanceToReferencePatients(
                queryPatient,
                agT,
                reactomeMechismoDataMap) * (double) agT.size();

        ag1d = ag1d.isNaN()
                ? 0.0d
                : ag1d;
        ag2d = ag2d.isNaN()
                ? 0.0d
                : ag2d;
        ag3d = ag3d.isNaN()
                ? 0.0d
                : ag3d;
        agTd = agTd.isNaN()
                ? 0.0d
                : agTd;

        Double ret = (ag1d + ag2d + ag3d + agTd) / (double) totalSize;

        if(ret.isNaN()){
            int debug = 1;
        }

        return ret;
    }

    public void writePatientDistancesToFile(String outputDir,
                                            String outputFilePrefix,
                                            ReactomeMechismoDataMap reactomeMechismoDataMap) {
        Map<Patient, List<Double>> patientToDistanceList = new HashMap<>();
        Map<Patient, List<Double>> patientToMembershipList = new HashMap<>();
        Map<Integer, List<List<Patient>>> targetReactionGroupClusters = getPatientsWith123TMutations();
        List<Integer> orderedClusterIDs = new ArrayList<>(targetReactionGroupClusters.keySet());
        Collections.sort(orderedClusterIDs);

        Integer patientCounter = 1;
        for (Patient patient : reactomeMechismoDataMap.getPatients()) {
            List<Double> distanceList = new ArrayList<>();
            List<Double> membershipList = new ArrayList<>();

            //all Target Reaction Groups

            //TARGET REACTION GROUP DISTANCES
            distanceList.add(calculateMeanReactionDistanceToReferencePatients(patient,
                    aggregateAllPatients(samplesW1MutatedUpstreamRxn),
                    reactomeMechismoDataMap));
            distanceList.add(calculateMeanReactionDistanceToReferencePatients(patient,
                    aggregateAllPatients(samplesW2MutatedUpstreamRxns),
                    reactomeMechismoDataMap));
            distanceList.add(calculateMeanReactionDistanceToReferencePatients(patient,
                    aggregateAllPatients(samplesW3plusMutatedUpstreamRxns),
                    reactomeMechismoDataMap));
            distanceList.add(calculateMeanReactionDistanceToReferencePatients(patient,
                    aggregateAllPatients(samplesWTargetRxnMutations),
                    reactomeMechismoDataMap));
            distanceList.add(getMeanDistanceToAllClusteredPatients(patient, reactomeMechismoDataMap));

            //TARGET REACTION GROUP MEMBERSHIP
            membershipList.add(aggregateAllPatients(samplesW1MutatedUpstreamRxn).contains(patient)
                    ? 1.0d
                    : 0.0d);
            membershipList.add(aggregateAllPatients(samplesW2MutatedUpstreamRxns).contains(patient)
                    ? 1.0d
                    : 0.0d);
            membershipList.add(aggregateAllPatients(samplesW3plusMutatedUpstreamRxns).contains(patient)
                    ? 1.0d
                    : 0.0d);
            membershipList.add(aggregateAllPatients(samplesWTargetRxnMutations).contains(patient)
                    ? 1.0d
                    : 0.0d);

            //each Target Reaction Group
            for (Integer clusterID : orderedClusterIDs) {
                //TARGET REACTION GROUP DISTANCES
                distanceList.add(calculateMeanReactionDistanceToReferencePatients(patient,
                        targetReactionGroupClusters.get(clusterID).get(0),
                        reactomeMechismoDataMap));
                distanceList.add(calculateMeanReactionDistanceToReferencePatients(patient,
                        targetReactionGroupClusters.get(clusterID).get(1),
                        reactomeMechismoDataMap));
                distanceList.add(calculateMeanReactionDistanceToReferencePatients(patient,
                        targetReactionGroupClusters.get(clusterID).get(2),
                        reactomeMechismoDataMap));
                distanceList.add(calculateMeanReactionDistanceToReferencePatients(patient,
                        targetReactionGroupClusters.get(clusterID).get(3),
                        reactomeMechismoDataMap));
                distanceList.add(getMeanDistanceToTargetReactionGroupAllClusteredPatients(patient,
                        reactomeMechismoDataMap,
                        targetReactionGroupClusters.get(clusterID)));

                //TARGET REACTION GROUP MEMBERSHIP
                membershipList.add(targetReactionGroupClusters.get(clusterID).get(0).contains(patient)
                        ? 1.0d
                        : 0.0d);
                membershipList.add(targetReactionGroupClusters.get(clusterID).get(1).contains(patient)
                        ? 1.0d
                        : 0.0d);
                membershipList.add(targetReactionGroupClusters.get(clusterID).get(2).contains(patient)
                        ? 1.0d
                        : 0.0d);
                membershipList.add(targetReactionGroupClusters.get(clusterID).get(3).contains(patient)
                        ? 1.0d
                        : 0.0d);
            }

            patientToDistanceList.put(patient, distanceList);
            patientToMembershipList.put(patient, membershipList);


            patientCounter++;
            if (patientCounter % 500 == 0) {
                System.out.println(String.format("Calculated %d of %d patient distances...",
                        patientCounter,
                        reactomeMechismoDataMap.getPatients().size()));
            }
        }

        String outFilePath = outputDir + outputFilePrefix + "TargetReactionGroupPatientDistances.csv";
        FileUtility fileUtility = new FileUtility();
        try {
            fileUtility.setOutput(outFilePath);
            StringBuilder headerLine = new StringBuilder("Patient Barcode,Cancer Type,");
            headerLine.append("Distance to Samples with 1 Mutated Upstream Reaction (all TRGs)," +
                    "Distance to Samples with 2 Mutated Upstream Reactions (all TRGs)," +
                    "Distance to Samples with 3+ Mutated Upstream Reactions (all TRGs)," +
                    "Distance to Samples with Mutated Target Reaction (all TRGs)," +
                    "Distance to Samples with 123T Mutated Reactions (all TRGs),");

            for(Integer clusterID : orderedClusterIDs){
             headerLine.append(String.format("Distance to Samples with 1 Mutated Upstream Reaction (TRG %1$d)," +
                             "Distance to Samples with 2 Mutated Upstream Reactions (TRG %1$d)," +
                             "Distance to Samples with 3+ Mutated Upstream Reactions (TRG %1$d)," +
                             "Distance to Samples with Mutated Target Reaction (TRG %1$d)," +
                             "Distance to Samples with 123T Mutated Reactions (TRG %1$d),",
                     clusterID));
            }

            headerLine.append("Has 1 Mutated Upstream Reaction (all TRGs)," +
                    "Has 2 Mutated Upstream Reactions (all TRGs)," +
                    "Has 3+ Mutated Upstream Reactions (all TRGs)," +
                    "Has Mutated Target Reaction (all TRGs),");

            for(Integer clusterID : orderedClusterIDs){

                headerLine.append(String.format("Has 1 Mutated Upstream Reaction (TRG %1$d)," +
                        "Has 2 Mutated Upstream Reactions (TRG %1$d)," +
                        "Has 3+ Mutated Upstream Reactions (TRG %1$d)," +
                        "Has Mutated Target Reaction (TRG %1$d),",
                        clusterID));
            }

            headerLine = new StringBuilder(headerLine.substring(0, headerLine.length() - 1));//remove trailing ','

            fileUtility.printLine(headerLine.toString());

            for (Patient patient : patientToDistanceList.keySet()) {
                StringBuilder patientData = new StringBuilder();
                for (int i = 0; i < patientToDistanceList.get(patient).size(); i++) {
                    patientData.append(String.format("%.100e,",
                            patientToDistanceList.get(patient).get(i)));
                }
                for (int i = 0; i < patientToMembershipList.get(patient).size(); i++) {
                    patientData.append(String.format("%.100e,",
                            patientToMembershipList.get(patient).get(i)));
                }
                patientData = new StringBuilder(patientData.substring(0, patientData.length() - 1));//remove trailing ','

                fileUtility.printLine(String.format("%s,%s,%s",
                        patient.getTcgaBarcode(),
                        patient.getCancerType(),
                        patientData.toString()));
            }
            fileUtility.close();
        } catch (IOException ioe) {
            System.out.println(String.format("Couldn't use %s, %s: %s",
                    outFilePath,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    private Map<Integer, List<List<Patient>>> getPatientsWith123TMutations() {
        if (this.patientsWith123TMutations == null) {
            Map<Integer, List<List<Patient>>> clusterIDToPatients = new HashMap<>();
            for (int i = 0; i < targetRxns.size(); i++) {
                if (pValue2EmpiricalPValueMap.get(pValues.get(i)) <= MAX_CLUSTER_EMPIRICAL_P_VALUE) {
                    Integer clusterID = Math.abs(cooccurringUpstreamRxns.get(i).hashCode());
                    List<List<Patient>> patientSetList;
                    if (!clusterIDToPatients.containsKey(clusterID)) {
                        patientSetList = new ArrayList<>();
                        patientSetList.add(new ArrayList<>(samplesW1MutatedUpstreamRxn.get(i)));
                        patientSetList.add(new ArrayList<>(samplesW2MutatedUpstreamRxns.get(i)));
                        patientSetList.add(new ArrayList<>(samplesW3plusMutatedUpstreamRxns.get(i)));
                        patientSetList.add(new ArrayList<>(samplesWTargetRxnMutations.get(i)));
                    }else{
                        patientSetList = clusterIDToPatients.get(clusterID);
                        Set<Patient> pw1 = new HashSet<>(patientSetList.get(0));
                        pw1.addAll(samplesW1MutatedUpstreamRxn.get(i));
                        Set<Patient> pw2 = new HashSet<>(patientSetList.get(1));
                        pw2.addAll(samplesW2MutatedUpstreamRxns.get(i));
                        Set<Patient> pw3 = new HashSet<>(patientSetList.get(2));
                        pw3.addAll(samplesW3plusMutatedUpstreamRxns.get(i));
                        Set<Patient> pwT = new HashSet<>(patientSetList.get(3));
                        pwT.addAll(samplesWTargetRxnMutations.get(i));
                        patientSetList.add(0,new ArrayList<>(pw1));
                        patientSetList.add(1,new ArrayList<>(pw2));
                        patientSetList.add(2,new ArrayList<>(pw3));
                        patientSetList.add(3,new ArrayList<>(pwT));
                    }
                    clusterIDToPatients.put(clusterID, patientSetList);
                }
            }
            this.patientsWith123TMutations = clusterIDToPatients;
        }
        return this.patientsWith123TMutations;
    }

    public void writePatientGroupingsToFile(String outputDir, String outputFilePrefix) {
        Map<Integer, List<List<Patient>>> clusterIDToPatients = getPatientsWith123TMutations();

        for (Integer clusterID : clusterIDToPatients.keySet()) {
            String outFilePath = outputDir + outputFilePrefix + "TargetReactionGroup" + clusterID + ".csv";
            FileUtility fileUtility = new FileUtility();
            List<Patient> patientGroup1 = clusterIDToPatients.get(clusterID).get(0);
            List<Patient> patientGroup2 = clusterIDToPatients.get(clusterID).get(1);
            List<Patient> patientGroup3 = clusterIDToPatients.get(clusterID).get(2);
            List<Patient> patientGroup4 = clusterIDToPatients.get(clusterID).get(3);
            int numLines = Math.max(
                    Math.max(patientGroup1.size(),
                            patientGroup2.size()),
                    patientGroup3.size());
            try {
                fileUtility.setOutput(outFilePath);
                fileUtility.printLine("1 Upstream Reaction Mutated," +
                        "2 Upstream Reactions Mutated," +
                        "3+ Upstream Reactions Mutated," +
                        "Target Reaction Mutated");
                for (int i = 0; i < numLines; i++) {
                    fileUtility.printLine(String.format("%s,%s,%s,%s",
                            patientGroup1.size() > i
                                    ? patientGroup1.get(i).getTcgaBarcode()
                                    : null,
                            patientGroup2.size() > i
                                    ? patientGroup2.get(i).getTcgaBarcode()
                                    : null,
                            patientGroup3.size() > i
                                    ? patientGroup3.get(i).getTcgaBarcode()
                                    : null,
                            patientGroup4.size() > i
                                    ? patientGroup4.get(i).getTcgaBarcode()
                                    : null));
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

    public void writeToFile(String outputDir, String outputFilePrefix) {
        String outFilePath = outputDir + outputFilePrefix + "rxnCooccurrence.csv";
        FileUtility fileUtility = new FileUtility();

        try {
            fileUtility.setOutput(outFilePath);
            fileUtility.printLine(
                    "Target Reaction," +
                            "#Co-occurring Upstream Reactions," +
                            "#Co-occurring Upstream Reaction (Mutated) FIs," +
                            "#Samples With Mutated Target Reaction," +
                            "#Samples With 0 Mutated Upstream Reactions," +
                            "#Samples With 1 Mutated Upstream Reaction," +
                            "#Samples With 2 Mutated Upstream Reactions," +
                            "#Samples With 3+ Mutated Upstream Reactions," +
                            "#Super-Indirect (Mutated) Genes," +
                            "#Indirect (Mutated) Genes," +
                            "#Super-Direct (Mutated) Genes," +
                            "#Direct (Mutated) Genes," +
                            "#Unique Super-Indirect Mutations," +
                            "#Unique Indirect Mutations," +
                            "#Unique Super-Direct Mutations," +
                            "#Unique Direct Mutations," +
                            "Co-occurring Upstream Reaction Cluster ID," +
                            "Co-occurring Upstream Reactions," +
                            "Co-occurring Upstream (Mutated) FIs," +
                            "Samples with Mutated Target Reaction," +
                            "Samples With 1 Mutated Upstream Reaction," +
                            "Samples With 2 Mutated Upstream Reactions," +
                            "Samples With 3+ Mutated Upstream Reactions," +
                            "Super-Indirect (Mutated) Genes," +
                            "Indirect (Mutated) Genes," +
                            "Super-Direct (Mutated) Genes," +
                            "Direct (Mutated) Genes," +
                            "Super-Indirect Mutations," +
                            "Indirect Mutations," +
                            "Super-Direct Mutations," +
                            "Direct Mutations," +
                            "Fishers Method Combined P-value," +
                            "BH Adjusted P-value," +
                            "Permutation-Based Empirical P-value");

            for (int i = 0; i < targetRxns.size(); i++) {
                WriteLineToFile(
                        fileUtility,
                        i);
            }
            fileUtility.close();
        } catch (IOException ioe) {
            System.out.println(String.format("Couldn't use %s, %s: %s",
                    outFilePath,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    private void WriteLineToFile(FileUtility fileUtility,
                                 int i) throws IOException {

        Double bhAdjustedP = (pValue2BHAdjustedPValueMap == null)
                ? 1.0d
                : pValue2BHAdjustedPValueMap.get(pValues.get(i));

        Double empiricalP = (pValue2EmpiricalPValueMap == null)
                ? 1.0d
                : pValue2EmpiricalPValueMap.get(pValues.get(i));

        //TODO: add sample support to entity classes

        fileUtility.printLine(String.format(
                "%s," + //Target Reaction
                        "%d," + //#Co-occurring Upstream Reactions
                        "%d," + //#Co-occurring Upstream Reaction FIs
                        "%d," + //#Samples With Mutated Target Reaction
                        "%d," + //#Samples With 0 Mutated Upstream Reactions
                        "%d," + //#Samples With 1 Mutated Upstream Reaction
                        "%d," + //#Samples With 2 Mutated Upstream Reactions
                        "%d," + //#Samples With 3+ Mutated Upstream Reactions
                        "%d," + //#Super-Indirect Genes
                        "%d," + //#Indirect Genes
                        "%d," + //#Super-Direct Genes
                        "%d," + //#Direct Genes
                        "%d," + //#Unique Super-Indirect Mutations
                        "%d," + //#Unique Indirect Mutations
                        "%d," + //#Unique Super-Direct Mutations
                        "%d," + //#Unique Direct Mutations
                        "%d," + //Co-occurring Upstream Reaction Cluster ID
                        "%s," + //Co-occurring Upstream Reactions (#Samples)
                        "%s," + //FIs (#Samples)
                        "%s," + //Samples With Mutated Target Reaction
                        "%s," + //Samples With 1 Mutated Upstream Reaction
                        "%s," + //Samples With 2 Mutated Upstream Reactions
                        "%s," + //Samples With 3+ Mutated Upstream Reactions
                        "%s," + //Super-Indirect Mutated Genes (#Samples)
                        "%s," + //Indirect Mutated Genes (#Samples)
                        "%s," + //Super-Direct Mutated Genes (#Samples)
                        "%s," + //Direct Mutated Genes (#Samples)
                        "%s," + //Super-Indirect Mutations (#Samples)
                        "%s," + //Indirect Mutations (#Samples)
                        "%s," + //Super-Direct Mutations (#Samples)
                        "%s," + //Direct Mutations (#Samples)
                        "%.100e," + //Fishers Method Combined P-value
                        "%.100e," + //BH Adjusted P-value
                        "%.100e", //Permutation-Based Empirical P-value
                targetRxns.get(i),
                cooccurringUpstreamRxns.get(i).size(),
                cooccurringUpstreamRxnFIs.get(i).size(),
                samplesWTargetRxnMutations.get(i).size(),
                numSamplesW0MutatedUpstreamRxns.get(i),
                samplesW1MutatedUpstreamRxn.get(i).size(),
                samplesW2MutatedUpstreamRxns.get(i).size(),
                samplesW3plusMutatedUpstreamRxns.get(i).size(),
                getSuperIndirectMutatedGenes().get(i).size(),
                getIndirectMutatedGenes().get(i).size(),
                getSuperDirectMutatedGenes().get(i).size(),
                getDirectMutatedGenes().get(i).size(),
                superIndirectMutations.get(i).size(),
                indirectMutations.get(i).size(),
                superDirectMutations.get(i).size(),
                directMutations.get(i).size(),
                Math.abs(cooccurringUpstreamRxns.get(i).hashCode()),
                cooccurringUpstreamRxns.get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(cooccurringUpstreamRxns.get(i)))
                        : Collections.singletonList(cooccurringUpstreamRxns.get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                cooccurringUpstreamRxnFIs.get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(cooccurringUpstreamRxnFIs.get(i)))
                        : Collections.singletonList(cooccurringUpstreamRxnFIs.get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                samplesWTargetRxnMutations.get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(samplesWTargetRxnMutations.get(i)))
                        : Collections.singletonList(samplesWTargetRxnMutations.get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                samplesW1MutatedUpstreamRxn.get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(samplesW1MutatedUpstreamRxn.get(i)))
                        : Collections.singletonList(samplesW1MutatedUpstreamRxn.get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                samplesW2MutatedUpstreamRxns.get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(samplesW2MutatedUpstreamRxns.get(i)))
                        : Collections.singletonList(samplesW2MutatedUpstreamRxns.get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                samplesW3plusMutatedUpstreamRxns.get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(samplesW3plusMutatedUpstreamRxns.get(i)))
                        : Collections.singletonList(samplesW3plusMutatedUpstreamRxns.get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                getSuperIndirectMutatedGenes().get(i).size() > 1
                        ? org.gk.util.StringUtils.join(" ", //space for copy-pasting
                        new ArrayList<>(getSuperIndirectMutatedGenes().get(i)))
                        : Collections.singletonList(getSuperIndirectMutatedGenes().get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                getIndirectMutatedGenes().get(i).size() > 1
                        ? org.gk.util.StringUtils.join(" ", //space for copy-pasting
                        new ArrayList<>(getIndirectMutatedGenes().get(i)))
                        : Collections.singletonList(getIndirectMutatedGenes().get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                getSuperDirectMutatedGenes().get(i).size() > 1
                        ? org.gk.util.StringUtils.join(" ", //space for copy-pasting
                        new ArrayList<>(getSuperDirectMutatedGenes().get(i)))
                        : Collections.singletonList(getSuperDirectMutatedGenes().get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                getDirectMutatedGenes().get(i).size() > 1
                        ? org.gk.util.StringUtils.join(" ", //space for copy-pasting
                        new ArrayList<>(getDirectMutatedGenes().get(i)))
                        : Collections.singletonList(getDirectMutatedGenes().get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                superIndirectMutations.get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(superIndirectMutations.get(i)))
                        : Collections.singletonList(superIndirectMutations.get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                indirectMutations.get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(indirectMutations.get(i)))
                        : Collections.singletonList(indirectMutations.get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                superDirectMutations.get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(superDirectMutations.get(i)))
                        : Collections.singletonList(superDirectMutations.get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                directMutations.get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(directMutations.get(i)))
                        : Collections.singletonList(directMutations.get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                pValues.get(i),
                bhAdjustedP,
                empiricalP));
    }
}
