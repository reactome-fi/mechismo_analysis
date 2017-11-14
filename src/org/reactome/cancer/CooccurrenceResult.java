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
    private Set<Reaction> patientCluster0Union;

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

        this.superIndirectMutatedGenes = null;
        this.indirectMutatedGenes = null;
        this.superDirectMutatedGenes = null;
        this.directMutatedGenes = null;
        this.pValue2BHAdjustedPValueMap = null;
        this.pValue2EmpiricalPValueMap = null;
        this.patientsWith123TMutations = null;
        this.patientCluster0Union = null;
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

    private List<Double> getpValues() {
        return pValues;
    }

    private List<Set<FI>> getCooccurringUpstreamRxnFIs() {
        return cooccurringUpstreamRxnFIs;
    }

    private List<Reaction> getTargetRxns() {
        return targetRxns;
    }

    private List<Set<Reaction>> getCooccurringUpstreamRxns() {
        return cooccurringUpstreamRxns;
    }

    private Map<Double, Double> getpValue2BHAdjustedPValueMap() {
        return pValue2BHAdjustedPValueMap;
    }

    private Map<Double, Double> getpValue2EmpiricalPValueMap() {
        return pValue2EmpiricalPValueMap;
    }

    private List<Integer> getNumSamplesW0MutatedUpstreamRxns() {
        return numSamplesW0MutatedUpstreamRxns;
    }

    private List<Set<Patient>> getSamplesW1MutatedUpstreamRxn() {
        return samplesW1MutatedUpstreamRxn;
    }

    private List<Set<Patient>> getSamplesW3plusMutatedUpstreamRxns() {
        return samplesW3plusMutatedUpstreamRxns;
    }

    private List<Set<Mutation>> getSuperIndirectMutations() {
        return superIndirectMutations;
    }

    private List<Set<Mutation>> getIndirectMutations() {
        return indirectMutations;
    }

    private List<Set<Mutation>> getSuperDirectMutations() {
        return superDirectMutations;
    }

    private List<Set<Mutation>> getDirectMutations() {
        return directMutations;
    }

    private List<Set<Patient>> getSamplesW2MutatedUpstreamRxns() {
        return samplesW2MutatedUpstreamRxns;
    }

    private Integer getRxnDistanceToPatientCluster0Union(Patient patient,
                                                         ReactomeMechismoDataMap reactomeMechismoDataMap) {
        if (this.patientCluster0Union == null) {
            getPatientCluster0Union(reactomeMechismoDataMap);
        }
        Set<Reaction> patientReactions = reactomeMechismoDataMap.getReactions(patient);

        if (patientReactions != null) {
            Set<Reaction> diff1 = new HashSet<>(patientReactions);
            diff1.removeAll(this.patientCluster0Union);

            Set<Reaction> diff2 = new HashSet<>(this.patientCluster0Union);
            diff2.removeAll(patientReactions);

            return diff1.size() + diff2.size();
        } else {
            return reactomeMechismoDataMap.getSupportedReactions().size();
        }
    }

    private Set<Reaction> getPatientCluster0Union(ReactomeMechismoDataMap reactomeMechismoDataMap) {
        if (this.patientCluster0Union == null) {
            if (this.patientsWith123TMutations == null) {
                getPatientsWith123TMutations();
            }
            List<List<Patient>> group0 = this.getPatientsWith123TMutations().get(0);
            Set<Reaction> patientCluster0Union = new HashSet<>();
            for (List<Patient> patientsInGroup : group0) {
                for (Patient patient : patientsInGroup) {
                    patientCluster0Union.addAll(reactomeMechismoDataMap.getReactions(patient));
                }
            }
            this.patientCluster0Union = patientCluster0Union;
        }
        return this.patientCluster0Union;
    }

    public void writePatientCluster0UnionDistancesToFile(String outputDir,
                                                         String outputFilePrefix,
                                                         ReactomeMechismoDataMap reactomeMechismoDataMap) {
        Map<Patient, Integer> patientToCluster0UnionDistance = new HashMap<>();
        Set<Patient> patients = reactomeMechismoDataMap.getPatients();
        Integer sumDistance = 0;
        Integer numPatients = patients.size();
        for (Patient patient : patients) {
            Integer distance = getRxnDistanceToPatientCluster0Union(patient, reactomeMechismoDataMap);
            sumDistance += distance;
            patientToCluster0UnionDistance.put(patient,
                    distance);
        }
        Integer meanDistance = sumDistance / numPatients;

        String outFilePath = outputDir + outputFilePrefix + "cluster0UnionDistances.csv";
        FileUtility fileUtility = new FileUtility();
        try {
            fileUtility.setOutput(outFilePath);
            fileUtility.printLine("Patient,Distance");
            for (Patient patient : patientToCluster0UnionDistance.keySet()) {
                if (patientToCluster0UnionDistance.get(patient) < meanDistance) {
                    fileUtility.printLine(String.format("%s,%d",
                            patient,
                            patientToCluster0UnionDistance.get(patient)
                    ));
                }
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
            Set<Patient> affectedBy1Patients = new HashSet<>();
            Set<Patient> affectedBy2Patients = new HashSet<>();
            Set<Patient> affectedBy3Patients = new HashSet<>();
            Set<Patient> affectedByTPatients = new HashSet<>();
            Map<Integer, List<List<Patient>>> clusterIDToPatients = new HashMap<>();
            for (int i = 0; i < getTargetRxns().size(); i++) {
                if (getpValue2EmpiricalPValueMap().get(getpValues().get(i)) <= MAX_CLUSTER_EMPIRICAL_P_VALUE) {
                    affectedBy1Patients.addAll(getSamplesW1MutatedUpstreamRxn().get(i));
                    affectedBy2Patients.addAll(getSamplesW2MutatedUpstreamRxns().get(i));
                    affectedBy3Patients.addAll(getSamplesW3plusMutatedUpstreamRxns().get(i));
                    affectedByTPatients.addAll(getSamplesWTargetRxnMutations().get(i));
                    Integer clusterID = Math.abs(getCooccurringUpstreamRxns().get(i).hashCode());
                    if (!clusterIDToPatients.containsKey(clusterID)) {
                        List<List<Patient>> patientSetList = new ArrayList<>();
                        patientSetList.add(new ArrayList<>(getSamplesW1MutatedUpstreamRxn().get(i)));
                        patientSetList.add(new ArrayList<>(getSamplesW2MutatedUpstreamRxns().get(i)));
                        patientSetList.add(new ArrayList<>(getSamplesW3plusMutatedUpstreamRxns().get(i)));
                        patientSetList.add(new ArrayList<>(getSamplesWTargetRxnMutations().get(i)));
                        clusterIDToPatients.put(clusterID, patientSetList);
                    }
                }
            }
            affectedBy1Patients.removeAll(affectedBy2Patients);
            affectedBy1Patients.removeAll(affectedBy3Patients);
            affectedBy1Patients.removeAll(affectedByTPatients);

            affectedBy2Patients.removeAll(affectedBy3Patients);
            affectedBy2Patients.removeAll(affectedByTPatients);

            affectedBy3Patients.removeAll(affectedByTPatients);

            List<List<Patient>> patientSetList = new ArrayList<>();
            patientSetList.add(new ArrayList<>(affectedBy1Patients));
            patientSetList.add(new ArrayList<>(affectedBy2Patients));
            patientSetList.add(new ArrayList<>(affectedBy3Patients));
            patientSetList.add(new ArrayList<>(affectedByTPatients));

            clusterIDToPatients.put(0, patientSetList);
            this.patientsWith123TMutations = clusterIDToPatients;
        }
        return this.patientsWith123TMutations;
    }

    public void writePatientGroupingsToFile(String outputDir, String outputFilePrefix) {
        //TODO: use a 'cancertype' argument as file prefix
        Map<Integer, List<List<Patient>>> clusterIDToPatients = getPatientsWith123TMutations();

        for (Integer clusterID : clusterIDToPatients.keySet()) {
            String outFilePath = outputDir + outputFilePrefix + "group" + clusterID + ".csv";
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
                fileUtility.printLine("1,2,3,4");
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

            for (int i = 0; i < getTargetRxns().size(); i++) {
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

        Double bhAdjustedP = (getpValue2BHAdjustedPValueMap() == null)
                ? 1.0d
                : getpValue2BHAdjustedPValueMap().get(getpValues().get(i));

        Double empiricalP = (getpValue2EmpiricalPValueMap() == null)
                ? 1.0d
                : getpValue2EmpiricalPValueMap().get(getpValues().get(i));

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
                getTargetRxns().get(i),
                getCooccurringUpstreamRxns().get(i).size(),
                getCooccurringUpstreamRxnFIs().get(i).size(),
                getSamplesWTargetRxnMutations().get(i).size(),
                getNumSamplesW0MutatedUpstreamRxns().get(i),
                getSamplesW1MutatedUpstreamRxn().get(i).size(),
                getSamplesW2MutatedUpstreamRxns().get(i).size(),
                getSamplesW3plusMutatedUpstreamRxns().get(i).size(),
                getSuperIndirectMutatedGenes().get(i).size(),
                getIndirectMutatedGenes().get(i).size(),
                getSuperDirectMutatedGenes().get(i).size(),
                getDirectMutatedGenes().get(i).size(),
                getSuperIndirectMutations().get(i).size(),
                getIndirectMutations().get(i).size(),
                getSuperDirectMutations().get(i).size(),
                getDirectMutations().get(i).size(),
                Math.abs(getCooccurringUpstreamRxns().get(i).hashCode()),
                getCooccurringUpstreamRxns().get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(getCooccurringUpstreamRxns().get(i)))
                        : Collections.singletonList(getCooccurringUpstreamRxns().get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                getCooccurringUpstreamRxnFIs().get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(getCooccurringUpstreamRxnFIs().get(i)))
                        : Collections.singletonList(getCooccurringUpstreamRxnFIs().get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                getSamplesWTargetRxnMutations().get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(getSamplesWTargetRxnMutations().get(i)))
                        : Collections.singletonList(getSamplesWTargetRxnMutations().get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                getSamplesW1MutatedUpstreamRxn().get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(getSamplesW1MutatedUpstreamRxn().get(i)))
                        : Collections.singletonList(getSamplesW1MutatedUpstreamRxn().get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                getSamplesW2MutatedUpstreamRxns().get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(getSamplesW2MutatedUpstreamRxns().get(i)))
                        : Collections.singletonList(getSamplesW2MutatedUpstreamRxns().get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                getSamplesW3plusMutatedUpstreamRxns().get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(getSamplesW3plusMutatedUpstreamRxns().get(i)))
                        : Collections.singletonList(getSamplesW3plusMutatedUpstreamRxns().get(i)).get(0).toString()
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
                getSuperIndirectMutations().get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(getSuperIndirectMutations().get(i)))
                        : Collections.singletonList(getSuperIndirectMutations().get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                getIndirectMutations().get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(getIndirectMutations().get(i)))
                        : Collections.singletonList(getIndirectMutations().get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                getSuperDirectMutations().get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(getSuperDirectMutations().get(i)))
                        : Collections.singletonList(getSuperDirectMutations().get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                getDirectMutations().get(i).size() > 1
                        ? org.gk.util.StringUtils.join("|",
                        new ArrayList<>(getDirectMutations().get(i)))
                        : Collections.singletonList(getDirectMutations().get(i)).get(0).toString()
                        .replace("[", "")
                        .replace("]", ""),
                getpValues().get(i),
                bhAdjustedP,
                empiricalP));
    }

    public List<Set<Patient>> getSamplesWTargetRxnMutations() {
        return samplesWTargetRxnMutations;
    }
}
