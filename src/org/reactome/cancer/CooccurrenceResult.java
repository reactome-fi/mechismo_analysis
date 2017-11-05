package org.reactome.cancer;

import org.reactome.cancer.driver.MechismoAnalyzer;
import org.reactome.r3.util.MathUtilities;

import java.util.*;

public class CooccurrenceResult {
    private List<Double> pValues;
    private List<Set<String>> fIs;
    private List<Long> targetRxns;
    private List<Set<Long>> upstreamRxns;
    private List<Set<String>> upstreamPairSamples;
    private List<Set<String>> excludedUpstreamPairSamples;
    private List<Set<List<String>>> upstreamRxnMutationsIncluded;
    private List<Set<List<String>>> upstreamRxnMutationsExlcuded;
    private Map<Double, Double> pValue2BHAdjustedPValueMap;
    private Map<Double, Double> pValue2EmpiricalPValueMap;

    public CooccurrenceResult(
            List<Double> allPValues,
            List<Set<String>> allFIs,
            List<Long> allTargetRxns,
            List<Set<Long>> allUpstreamRxns,
            List<Set<String>> allUpstreamPairSamples,
            List<Set<String>> allExcludedUpstreamPairSamples,
            List<Set<List<String>>> allUpstreamRxnMutationsIncluded,
            List<Set<List<String>>> allUpstreamRxnMutationsExlcuded) {
        this.pValues = allPValues;
        this.fIs = allFIs;
        this.targetRxns = allTargetRxns;
        this.upstreamRxns = allUpstreamRxns;
        this.upstreamPairSamples = allUpstreamPairSamples;
        this.excludedUpstreamPairSamples = allExcludedUpstreamPairSamples;
        this.upstreamRxnMutationsIncluded = allUpstreamRxnMutationsIncluded;
        this.upstreamRxnMutationsExlcuded = allUpstreamRxnMutationsExlcuded;
    }

    public int getTargetRxnCount() {
        return this.targetRxns.size();
    }

    public void CalculateBHAdjustedPValues() {
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


        List<Double> rewiredNetworkPvalues = new ArrayList<>();

        double p5 = this.pValues.size() * 1.5;
        double m5 = this.pValues.size() * 0.5;
        //loop over rewirings
        for (CooccurrenceResult rewiredNetworkResult : rewiredNetworkResults) {

            if (rewiredNetworkResult.pValues.size() > p5 ||
                    rewiredNetworkResult.pValues.size() < m5) {
                System.out.println(String.format(
                        "The number of p values from this permutation was %d but should be closer to %d",
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
        if (fdr > 1.0)
            fdr = 1.0;
        return fdr;
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
        //leave only pValues
        this.fIs.clear();
        this.targetRxns.clear();
        this.upstreamRxns.clear();
        this.upstreamPairSamples.clear();
        this.excludedUpstreamPairSamples.clear();
        this.upstreamRxnMutationsIncluded.clear();
        this.upstreamRxnMutationsExlcuded.clear();
        System.gc();
    }

    public List<Double> getpValues() {
        return pValues;
    }

    public List<Set<String>> getfIs() {
        return fIs;
    }

    public List<Long> getTargetRxns() {
        return targetRxns;
    }

    public List<Set<Long>> getUpstreamRxns() {
        return upstreamRxns;
    }

    public List<Set<String>> getUpstreamPairSamples() {
        return upstreamPairSamples;
    }

    public List<Set<String>> getExcludedUpstreamPairSamples() {
        return excludedUpstreamPairSamples;
    }

    public List<Set<List<String>>> getUpstreamRxnMutationsIncluded() {
        return upstreamRxnMutationsIncluded;
    }

    public List<Set<List<String>>> getUpstreamRxnMutationsExlcuded() {
        return upstreamRxnMutationsExlcuded;
    }

    public Map<Double, Double> getpValue2BHAdjustedPValueMap() {
        return pValue2BHAdjustedPValueMap;
    }

    public Map<Double, Double> getpValue2EmpiricalPValueMap() {
        return pValue2EmpiricalPValueMap;
    }
}
