package org.reactome.cancer;

import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.MathUtilities;

import java.io.IOException;
import java.util.*;

public class CooccurrenceResult {
    private List<Reaction> targetRxns;
    private List<Set<Reaction>> cooccurringUpstreamRxns;
    private List<Set<FI>> cooccurringUpstreamRxnFIs;
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

    public CooccurrenceResult(
            List<Reaction> targetRxns,
            List<Set<Reaction>> cooccurringUpstreamRxns,
            List<Set<FI>> cooccurringUpstreamRxnFIs,
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
        this.targetRxns.clear();
        this.cooccurringUpstreamRxns.clear();
        this.cooccurringUpstreamRxnFIs.clear();
        this.numSamplesW0MutatedUpstreamRxns.clear();
        this.samplesW1MutatedUpstreamRxn.clear();
        this.samplesW3plusMutatedUpstreamRxns.clear();
        this.superIndirectMutations.clear();
        this.indirectMutations.clear();
        this.superDirectMutations.clear();
        this.directMutations.clear();
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

    public void writeToFile(String outputDir,String outputFilePrefix){
        String outFilePath5 = outputDir + outputFilePrefix + "rxnCooccurrence.csv";
        FileUtility fileUtility = new FileUtility();

        try {
            fileUtility.setOutput(outFilePath5);
            fileUtility.printLine(
                    "Target Reaction," +
                            "#Co-occurring Upstream Reactions," +
                            "#Co-occurring Upstream Reaction (Mutated) FIs," +
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
                    outFilePath5,
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
                getCooccurringUpstreamRxns().get(i).hashCode(),
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
}
