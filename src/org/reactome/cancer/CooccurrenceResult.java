package org.reactome.cancer;

import org.reactome.r3.util.MathUtilities;

import java.util.*;

public class CooccurrenceResult {
    private List<Long> targetRxns;
    private List<Set<Long>> cooccurringUpstreamRxns;
    private List<Set<String>> cooccurringUpstreamRxnFIs;
    private List<Integer> numSamplesW0MutatedUpstreamRxns;
    private List<Set<String>> samplesW1MutatedUpstreamRxn;
    private List<Set<String>> samplesW2MutatedUpstreamRxns;
    private List<Set<String>> samplesW3plusMutatedUpstreamRxns;
    private List<Set<List<String>>> superIndirectMutations;
    private List<Set<List<String>>> indirectMutations;
    private List<Set<List<String>>> superDirectMutations;
    private List<Set<List<String>>> directMutations;
    private List<Double> pValues;
    private List<Set<String>> superIndirectMutatedGenes;
    private List<Set<String>> indirectMutatedGenes;
    private List<Set<String>> superDirectMutatedGenes;
    private List<Set<String>> directMutatedGenes;
    private Map<Double, Double> pValue2BHAdjustedPValueMap;
    private Map<Double, Double> pValue2EmpiricalPValueMap;

    public CooccurrenceResult(
            List<Long> targetRxns,
            List<Set<Long>> cooccurringUpstreamRxns,
            List<Set<String>> cooccurringUpstreamRxnFIs,
            List<Integer> numSamplesW0MutatedUpstreamRxns,
            List<Set<String>> samplesW1MutatedUpstreamRxn,
            List<Set<String>> samplesW2MutatedUpstreamRxns,
            List<Set<String>> samplesW3plusMutatedUpstreamRxns,
            List<Set<List<String>>> superIndirectMutations,
            List<Set<List<String>>> indirectMutations,
            List<Set<List<String>>> superDirectMutations,
            List<Set<List<String>>> directMutations,
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
            Set<String> superIndirectTargetMutatatedGenes = new HashSet<>();
            Set<String> indirectTargetMutatedGenes = new HashSet<>();
            Set<String> superDirectTargetMutatedGenes = new HashSet<>();
            Set<String> directTargetMutatedGenes = new HashSet<>();
            for (List<String> mutation : this.superIndirectMutations.get(i)) {
                superIndirectTargetMutatatedGenes.add(mutation.get(0));
            }
            for (List<String> mutation : this.indirectMutations.get(i)) {
                indirectTargetMutatedGenes.add(mutation.get(0));
            }
            for (List<String> mutation : this.superDirectMutations.get(i)) {
                superDirectTargetMutatedGenes.add(mutation.get(0));
            }
            for (List<String> mutation : this.directMutations.get(i)) {
                directTargetMutatedGenes.add(mutation.get(0));
            }
            this.superIndirectMutatedGenes.add(superIndirectTargetMutatatedGenes);
            this.indirectMutatedGenes.add(indirectTargetMutatedGenes);
            this.superDirectMutatedGenes.add(superDirectTargetMutatedGenes);
            this.directMutatedGenes.add(directTargetMutatedGenes);
        }
    }

    public List<Set<String>> getSuperIndirectMutatedGenes() {
        if (this.superIndirectMutatedGenes == null) {
            FillMutatedGeneSets();
        }
        return this.superIndirectMutatedGenes;
    }

    public List<Set<String>> getIndirectMutatedGenes() {
        if (this.indirectMutatedGenes == null) {
            FillMutatedGeneSets();
        }
        return this.indirectMutatedGenes;
    }

    public List<Set<String>> getSuperDirectMutatedGenes() {
        if (this.superDirectMutatedGenes == null) {
            FillMutatedGeneSets();
        }
        return this.superDirectMutatedGenes;
    }

    public List<Set<String>> getDirectMutatedGenes() {
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

    public List<Double> getpValues() {
        return pValues;
    }

    public List<Set<String>> getCooccurringUpstreamRxnFIs() {
        return cooccurringUpstreamRxnFIs;
    }

    public List<Long> getTargetRxns() {
        return targetRxns;
    }

    public List<Set<Long>> getCooccurringUpstreamRxns() {
        return cooccurringUpstreamRxns;
    }

    public Map<Double, Double> getpValue2BHAdjustedPValueMap() {
        return pValue2BHAdjustedPValueMap;
    }

    public Map<Double, Double> getpValue2EmpiricalPValueMap() {
        return pValue2EmpiricalPValueMap;
    }

    public List<Integer> getNumSamplesW0MutatedUpstreamRxns() {
        return numSamplesW0MutatedUpstreamRxns;
    }

    public List<Set<String>> getSamplesW1MutatedUpstreamRxn() {
        return samplesW1MutatedUpstreamRxn;
    }

    public List<Set<String>> getSamplesW3plusMutatedUpstreamRxns() {
        return samplesW3plusMutatedUpstreamRxns;
    }

    public List<Set<List<String>>> getSuperIndirectMutations() {
        return superIndirectMutations;
    }

    public List<Set<List<String>>> getIndirectMutations() {
        return indirectMutations;
    }

    public List<Set<List<String>>> getSuperDirectMutations() {
        return superDirectMutations;
    }

    public List<Set<List<String>>> getDirectMutations() {
        return directMutations;
    }

    public List<Set<String>> getSamplesW2MutatedUpstreamRxns() {
        return samplesW2MutatedUpstreamRxns;
    }
}
