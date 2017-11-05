package org.reactome.cancer;

import org.reactome.cancer.driver.MechismoAnalyzer;

import java.util.Map;
import java.util.Objects;
import java.util.Set;

public class TargetReactionSummary {
    private static final String headerLine =
            "RxnId," +
                    "Num Sup Dn/Up Rxns," +
                    "Num Sup Dn/Up Combos," +
                    "Num All Dn/Up Rxns," +
                    "Support Ratio," +
                    "Sup Dn/Up Rxns," +
                    "Sup Dn/Up Combos," +
                    "All Dn/Up Rxns," +
                    "Supporting FIs";
    private Long rxnId;
    private Integer numSupportedUpstreamRxns;
    private Integer numUpstreamRxns;
    private Double supportedUpstreamRxnRatio;
    private Set<Long> supportedUpstreamRxns;
    private Set<Long> upstreamRxns;
    private Set<String> supportingFIs;
    private Map<Long, String> longRxnDbIdToName;

    public TargetReactionSummary(Long rxnId,
                                 Integer numSupportedUpstreamRxns,
                                 Integer numUpstreamRxns,
                                 Double supportedUpstreamRxnRatio,
                                 Set<Long> supportedUpstreamRxns,
                                 Set<Long> upstreamRxns,
                                 Set<String> supportingFIs,
                                 Map<Long, String> longRxnDbIdToName) {
        this.rxnId = rxnId;
        this.numSupportedUpstreamRxns = numSupportedUpstreamRxns;
        this.numUpstreamRxns = numUpstreamRxns;
        this.supportedUpstreamRxnRatio = supportedUpstreamRxnRatio;
        this.supportedUpstreamRxns = supportedUpstreamRxns;
        this.upstreamRxns = upstreamRxns;
        this.supportingFIs = supportingFIs;
        this.longRxnDbIdToName = longRxnDbIdToName;
    }

    public static String getHeaderLine() {
        return headerLine;
    }

    @Override
    public boolean equals(Object o) {
        if (o == this) return true;
        if (!(o instanceof TargetReactionSummary)) {
            return false;
        }
        TargetReactionSummary targetReactionSummary = (TargetReactionSummary) o;
        return Objects.equals(this.rxnId, targetReactionSummary.rxnId) &&
                Objects.equals(this.numSupportedUpstreamRxns, targetReactionSummary.numSupportedUpstreamRxns) &&
                Objects.equals(this.numUpstreamRxns, targetReactionSummary.numUpstreamRxns) &&
                Objects.equals(this.supportedUpstreamRxnRatio, targetReactionSummary.supportedUpstreamRxnRatio) &&
                Objects.equals(this.supportedUpstreamRxns, targetReactionSummary.supportedUpstreamRxns) &&
                Objects.equals(this.upstreamRxns, targetReactionSummary.upstreamRxns) &&
                Objects.equals(this.supportingFIs, targetReactionSummary.supportingFIs);
    }

    @Override
    public int hashCode() {
        return Objects.hash(this.rxnId,
                this.numSupportedUpstreamRxns,
                this.numUpstreamRxns,
                this.supportedUpstreamRxnRatio,
                this.supportedUpstreamRxns,
                this.supportingFIs);
    }

    public Long getRxnId() {
        return rxnId;
    }

    public Integer getNumSupportedUpstreamRxns() {
        return numSupportedUpstreamRxns;
    }

    public Integer getNumUpstreamRxns() {
        return numUpstreamRxns;
    }

    public Double getSupportedUpstreamRxnRatio() {
        return supportedUpstreamRxnRatio;
    }

    public Set<Long> getSupportedUpstreamRxns() {
        return supportedUpstreamRxns;
    }

    public Set<Long> getUpstreamRxns() {
        return upstreamRxns;
    }

    public Set<String> getSupportingFIs() {
        return supportingFIs;
    }

    public Map<Long, String> getLongRxnDbIdToName() {
        return longRxnDbIdToName;
    }
}
