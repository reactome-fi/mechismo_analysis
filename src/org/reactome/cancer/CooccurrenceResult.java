package org.reactome.cancer;

import org.reactome.r3.util.MathUtilities;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;

public class CooccurrenceResult {
    private List<BigDecimal> pValues;
    private List<Set<String>> fIs;
    private List<Long> targetRxns;
    private List<Set<Long>> upstreamRxns;
    private List<Set<String>> upstreamPairSamples;
    private List<Set<String>> excludedUpstreamPairSamples;
    private List<Set<List<String>>> upstreamRxnMutationsIncluded;
    private List<Set<List<String>>> upstreamRxnMutationsExlcuded;
    private Map<BigDecimal, BigDecimal> pValue2BHAdjustedPValueMap;
    private Map<BigDecimal, BigDecimal> pValue2EmpiricalPValueMap;

    public CooccurrenceResult(
            List<BigDecimal> allPValues,
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
        List<BigDecimal> pValuesSorted = new ArrayList<>(this.pValues);

        List<BigDecimal> BHAdjustedPValues =
                MathUtilities.calculateFDRWithBenjaminiHochbergBD(pValuesSorted);

        this.pValue2BHAdjustedPValueMap = new HashMap<>();
        for (int i = 0; i < pValues.size(); i++) {
            this.pValue2BHAdjustedPValueMap.put(pValuesSorted.get(i), BHAdjustedPValues.get(i));
        }
    }

    //GSEA FDR Calculation
    public void CalculateEmpiricalPValues(List<CooccurrenceResult> rewiredNetworkResults) {


        List<BigDecimal> rewiredNetworkPvalues = new ArrayList<>();

        double p5 = this.pValues.size() * 1.5;
        double m5 = this.pValues.size() * 0.5;
        //loop over rewirings
        for (CooccurrenceResult rewiredNetworkResult : rewiredNetworkResults) {

            if (rewiredNetworkResult.pValues.size() > p5 ||
                    rewiredNetworkResult.pValues.size() < m5) {
                System.out.println(String.format(
                        "Permutation had %d results but real had %d",
                        rewiredNetworkResult.pValues.size(),
                        this.pValues.size()));
            }

            rewiredNetworkPvalues.addAll(rewiredNetworkResult.pValues);
        }

        this.pValue2EmpiricalPValueMap = new HashMap<>();

        //loop over pValues
        for (BigDecimal pValue : this.pValues) {
            this.pValue2EmpiricalPValueMap.put(pValue,
                    CalculateEmpiricalPValue(
                            CountValuesLTEThresh(rewiredNetworkPvalues, pValue),
                            CountValuesLTEThresh(this.pValues, pValue),
                            rewiredNetworkPvalues.size(),
                            this.pValues.size()));
        }
    }

    private BigDecimal CalculateEmpiricalPValue(int randomLessThan,
                                                int realLessThan,
                                                int randomTotal,
                                                int realTotal) {

        BigDecimal randomLessThanBD = new BigDecimal(randomLessThan);
        BigDecimal realLessThanBD = new BigDecimal(realLessThan);
        BigDecimal randomTotalBD = new BigDecimal(randomTotal);
        BigDecimal realTotalBD = new BigDecimal(realTotal);
        BigDecimal fdr = null;
        try {
            BigDecimal num = randomLessThanBD.divide(realLessThanBD, RoundingMode.HALF_UP);
            BigDecimal denom = randomTotalBD.divide(realTotalBD, RoundingMode.HALF_UP);
            fdr = num.divide(denom, RoundingMode.HALF_UP);
        }catch(ArithmeticException ae){
            int debug = 1;
        }
        if (fdr.compareTo(BigDecimal.ONE) > 0)
            fdr = BigDecimal.ONE;
        return fdr;
    }

    private int CountValuesLTEThresh(List<BigDecimal> vals, BigDecimal threshold) {
        int count = 0;
        for (BigDecimal val : vals) {
            count = val.compareTo(threshold) <= 0 ?
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

    public List<BigDecimal> getpValues() {
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

    public Map<BigDecimal, BigDecimal> getpValue2BHAdjustedPValueMap() {
        return pValue2BHAdjustedPValueMap;
    }

    public Map<BigDecimal, BigDecimal> getpValue2EmpiricalPValueMap() {
        return pValue2EmpiricalPValueMap;
    }
}
