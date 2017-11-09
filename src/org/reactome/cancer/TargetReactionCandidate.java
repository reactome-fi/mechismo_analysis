package org.reactome.cancer;

import java.util.Objects;
import java.util.Set;

public class TargetReactionCandidate {
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
    private Reaction targetReaction;
    private Integer numUpstreamRxns;
    private Set<Reaction> supportedUpstreamRxns;
    private Set<Reaction> upstreamRxns;
    private Set<FI> supportingFIs;

    public TargetReactionCandidate(Reaction targetReaction,
                                   Set<Reaction> supportedUpstreamRxns,
                                   Set<Reaction> upstreamRxns,
                                   Set<FI> supportingFIs) {
        this.targetReaction = targetReaction;
        this.supportedUpstreamRxns = supportedUpstreamRxns;
        this.upstreamRxns = upstreamRxns;
        this.supportingFIs = supportingFIs;
    }

    public static String getHeaderLine() {
        return headerLine;
    }

    @Override
    public boolean equals(Object o) {
        if (o == this) return true;
        if (!(o instanceof TargetReactionCandidate)) {
            return false;
        }
        TargetReactionCandidate targetReactionCandidate = (TargetReactionCandidate) o;
        return Objects.equals(this.targetReaction, targetReactionCandidate.targetReaction) &&
                Objects.equals(this.numUpstreamRxns, targetReactionCandidate.numUpstreamRxns) &&
                Objects.equals(this.supportedUpstreamRxns, targetReactionCandidate.supportedUpstreamRxns) &&
                Objects.equals(this.upstreamRxns, targetReactionCandidate.upstreamRxns) &&
                Objects.equals(this.supportingFIs, targetReactionCandidate.supportingFIs);
    }

    @Override
    public int hashCode() {
        return Objects.hash(this.targetReaction,
                this.numUpstreamRxns,
                this.supportedUpstreamRxns,
                this.supportingFIs);
    }

    public Reaction getTargetReaction() {
        return targetReaction;
    }

    public Double getSupportedUpstreamRxnRatio() {
        return (double) this.supportedUpstreamRxns.size() / (double) this.upstreamRxns.size();
    }

    public Set<Reaction> getSupportedUpstreamRxns() {
        return supportedUpstreamRxns;
    }

    public Set<Reaction> getUpstreamRxns() {
        return upstreamRxns;
    }

    public Set<FI> getSupportingFIs() {
        return supportingFIs;
    }
}
