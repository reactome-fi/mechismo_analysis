package org.reactome.cancer;

import org.apache.commons.math.MathException;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.reactome.r3.util.FisherExact;
import org.reactome.r3.util.MathUtilities;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class ReactionGraphAnalyzer {
    private Integer depth;
    private FisherExact fisherExact;
    private ReactomeMechismoDataMap reactomeMechismoDataMap;
    private boolean ignoreDependentUpstreamReactions;
    private boolean ignoreIndependentUpstreamReactions;
    private boolean excludeMultipleImmediateUpstreamReactions;

    public ReactionGraphAnalyzer(Integer depth,
                                 ReactomeMechismoDataMap reactomeMechismoDataMap,
                                 boolean ignoreDependentUpstreamReactions,
                                 boolean ignoreIndependentUpstreamReactions,
                                 boolean excludeMultipleImmediateUpstreamReactions) {
        this.depth = depth;
        this.fisherExact = new FisherExact(reactomeMechismoDataMap.getNumPatients());
        this.reactomeMechismoDataMap = reactomeMechismoDataMap;
        this.ignoreDependentUpstreamReactions = ignoreDependentUpstreamReactions;
        this.ignoreIndependentUpstreamReactions = ignoreIndependentUpstreamReactions;
        this.excludeMultipleImmediateUpstreamReactions = excludeMultipleImmediateUpstreamReactions;
    }

    public CooccurrenceResult SearchRxnNetworkAndCalculateCooccurrencePValues(
            DirectedGraph<Long, DefaultEdge> reactionGraph,
            Boolean useRxnDist,
            Integer minNumTargetRxnPatients) throws MathException {

        Set<TargetReactionCandidate> targetReactionCandidates = SearchRxnNetworkForTargetReactionCandidates(
                reactionGraph);

        System.out.println(String.format("Found %d supported Dn/Up reactions",
                targetReactionCandidates.size()));

        return CalculateCooccurrencePValues(
                targetReactionCandidates,
                reactionGraph,
                useRxnDist,
                minNumTargetRxnPatients);
    }

    private Set<TargetReactionCandidate> SearchRxnNetworkForTargetReactionCandidates(
            DirectedGraph<Long, DefaultEdge> reactionGraph) {
        Set<DefaultEdge> outgoingEdges = new HashSet<>();
        Set<Reaction> downstreamReactions = new HashSet<>();
        Set<TargetReactionCandidate> targetReactionCandidates = new HashSet<>();

        for (Reaction reaction : reactomeMechismoDataMap.getSupportedReactions()) {
            //http://jgrapht.org/javadoc/org/jgrapht/graph/DefaultDirectedGraph.html
            //go one reaction downstream
            outgoingEdges.addAll(reactionGraph.outgoingEdgesOf(reaction.getReactionID()));
        }
        for (DefaultEdge outgoingEdge : outgoingEdges) {
            downstreamReactions.add(
                    reactomeMechismoDataMap.getReaction(
                            reactionGraph.getEdgeTarget(outgoingEdge)));
        }
        for (Reaction downstreamReaction : downstreamReactions) {
            targetReactionCandidates.add(CreateTargetReactionCandidate(
                    reactionGraph,
                    downstreamReaction));
        }
        return targetReactionCandidates;
    }

    private TargetReactionCandidate CreateTargetReactionCandidate(
            DirectedGraph<Long, DefaultEdge> reactionGraph,
            Reaction downstreamReaction) {
        Set<DefaultEdge> incomingEdgesCpy = new HashSet<>(reactionGraph.incomingEdgesOf(downstreamReaction.getReactionID()));
        Set<Reaction> upstreamReactions = new HashSet<>();
        Set<Reaction> allUpstreamReactions = new HashSet<>();
        List<Set<Reaction>> upstreamReactionsByDepth = new ArrayList<>();
        for (int i = 1; i <= depth; i++) {
            for (DefaultEdge incomingEdge : incomingEdgesCpy) {
                upstreamReactions.add(reactomeMechismoDataMap.getReaction(reactionGraph.getEdgeSource(incomingEdge)));
            }
            upstreamReactions.removeAll(allUpstreamReactions);
            allUpstreamReactions.addAll(upstreamReactions);
            if (i < depth) {
                upstreamReactionsByDepth.add(new HashSet<>(upstreamReactions));
                upstreamReactions.clear();
                incomingEdgesCpy.clear();
                for (Reaction upstreamReaction : upstreamReactionsByDepth.get(i - 1)) {
                    incomingEdgesCpy.addAll(reactionGraph.incomingEdgesOf(upstreamReaction.getReactionID()));
                }
            }
        }
        incomingEdgesCpy.clear();
        upstreamReactions.clear();
        upstreamReactionsByDepth.clear();

        allUpstreamReactions.remove(downstreamReaction);

        //find intersection of upstream reactions and supported reactions
        Set<Reaction> supportedUpstreamReactions = new HashSet<>(allUpstreamReactions);
        supportedUpstreamReactions.retainAll(reactomeMechismoDataMap.getSupportedReactions());

        //remove upstream reactions directly activated by target
        supportedUpstreamReactions.removeIf(
                upstreamReaction -> !reactionGraph.getAllEdges(
                        downstreamReaction.getReactionID(),
                        upstreamReaction.getReactionID()).isEmpty());

        Set<FI> supportedFIs = new HashSet<>();
        for (Reaction supportedUpstreamReaction : supportedUpstreamReactions) {
            supportedFIs.addAll(reactomeMechismoDataMap.getFIs(supportedUpstreamReaction));
        }

        return new TargetReactionCandidate(
                downstreamReaction,
                supportedUpstreamReactions,
                supportedFIs);
    }

    private CooccurrenceResult CalculateCooccurrencePValues(
            Set<TargetReactionCandidate> targetReactionCandidates,
            DirectedGraph<Long, DefaultEdge> reactionGraph,
            Boolean useRxnDist,
            Integer minNUmTargetRxnPatients) throws MathException {

        List<Reaction> allTargetRxns = new ArrayList<>();
        List<Set<Reaction>> allCooccurringUpstreamRxns = new ArrayList<>();
        List<Set<FI>> allCooccurringUpstreamReactionFIs = new ArrayList<>();
        List<Set<Patient>> allSamplesWTargetRxnMutations = new ArrayList<>();
        List<Integer> allSamplesW0MutatedUpstreamRxns = new ArrayList<>();
        List<Set<Patient>> allSamplesW1MutatedUpstreamRxn = new ArrayList<>();
        List<Set<Patient>> allSamplesW2MutatedUpstreamRxns = new ArrayList<>();
        List<Set<Patient>> allSamplesW3plusMutatedUpstreamRxns = new ArrayList<>();
        List<Set<Mutation>> allSuperIndirectMutatations = new ArrayList<>();
        List<Set<Mutation>> allIndirectMutations = new ArrayList<>();
        List<Set<Mutation>> allSuperDirectMutations = new ArrayList<>();
        List<Set<Mutation>> allDirectMutations = new ArrayList<>();
        List<Double> allPValues = new ArrayList<>();

        int iterationCounter = 0;
        for (TargetReactionCandidate targetReactionCandidate : targetReactionCandidates) {
            Reaction targetRxn = targetReactionCandidate.getTargetReaction();
            List<Reaction[]> targetUpstreamReactionPairs =
                    new ArrayList<>(GenerateReactionSetPairs(targetReactionCandidate.getSupportedUpstreamRxns()));
            Set<Reaction> targetCooccurringUpstreamRxns = new HashSet<>();
            Set<FI> targetCooccurringUpstreamReactionFIs = new HashSet<>();
            Set<Patient> targetSamplesWTargetRxnMutations = new HashSet<>();
            Set<Patient> targetSamplesW0MutatedUpstreamRxns = new HashSet<>();
            Set<Patient> targetSamplesW1MutatedUpstreamRxn = new HashSet<>();
            Set<Patient> targetSamplesW2MutatedUpstreamRxns = new HashSet<>();
            Set<Patient> targetSamplesW3plusMutatedUpstreamRxns = new HashSet<>();
            Set<Mutation> targetSuperIndirectMutations = new HashSet<>();
            Set<Mutation> targetIndirectMutations = new HashSet<>();
            Set<Mutation> targetSuperDirectMutations = new HashSet<>();
            Set<Mutation> targetDirectMutations = new HashSet<>();
            List<Double> targetUpstreamReactionPValues = new ArrayList<>();
            List<List<Double>> targetUpstreamReactionABCDs = new ArrayList<>();

            if (!targetUpstreamReactionPairs.isEmpty()) {
                boolean targetReactionDetected = false;
                for (Reaction[] upstreamReactionPair : targetUpstreamReactionPairs) {
                    Reaction upstreamReaction1 = upstreamReactionPair[0];
                    Reaction upstreamReaction2 = upstreamReactionPair[1];

                    if (ignoreDependentUpstreamReactions &&
                            (!reactionGraph.getAllEdges(upstreamReaction1.getReactionID(),
                                    upstreamReaction2.getReactionID()).isEmpty() ||
                                    !reactionGraph.getAllEdges(upstreamReaction2.getReactionID(),
                                            upstreamReaction1.getReactionID()).isEmpty())) {
                        continue;
                    } else if (ignoreIndependentUpstreamReactions &&
                            (reactionGraph.getAllEdges(upstreamReaction1.getReactionID(),
                                    upstreamReaction2.getReactionID()).isEmpty() &&
                                    reactionGraph.getAllEdges(upstreamReaction2.getReactionID(),
                                            upstreamReaction1.getReactionID()).isEmpty())) {
                        continue;
                    } else if (excludeMultipleImmediateUpstreamReactions &&
                            (!reactionGraph.getAllEdges(upstreamReaction1.getReactionID(),
                                    targetRxn.getReactionID()).isEmpty() &&
                                    !reactionGraph.getAllEdges(upstreamReaction2.getReactionID(),
                                            targetRxn.getReactionID()).isEmpty())) {
                        continue;
                    } else {
                        targetReactionDetected = true;

                        Set<Patient> A = new HashSet<>(reactomeMechismoDataMap.getPatients());
                        Set<Patient> B = new HashSet<>(reactomeMechismoDataMap.getPatients());
                        Set<Patient> C = new HashSet<>(reactomeMechismoDataMap.getPatients());
                        Set<Patient> D = new HashSet<>(reactomeMechismoDataMap.getPatients());
                        Set<Patient> upstreamRxn1Samples = new HashSet<>(reactomeMechismoDataMap.getPatients(upstreamReaction1));
                        Set<Patient> upstreamRxn2Samples = new HashSet<>(reactomeMechismoDataMap.getPatients(upstreamReaction2));

                        D.retainAll(upstreamRxn1Samples);
                        D.retainAll(upstreamRxn2Samples);

                        C.removeAll(D);
                        C.retainAll(upstreamRxn2Samples);

                        B.removeAll(D);
                        B.retainAll(upstreamRxn1Samples);

                        A.removeAll(B);
                        A.removeAll(C);
                        A.removeAll(D);

                        if ((A.size() + B.size() + C.size() + D.size()) != reactomeMechismoDataMap.getPatients().size()) {
                            throw new IllegalStateException(String.format(
                                    "A(%d) + B(%d) + C(%d) + D(%d) sum to %d but should sum to %d",
                                    A.size(),
                                    B.size(),
                                    C.size(),
                                    D.size(),
                                    (A.size() + B.size() + C.size() + D.size()),
                                    reactomeMechismoDataMap.getPatients().size()));
                        }

                        Set<Mutation> targetRxnmuts = new HashSet<>();
                        Set<Mutation> upstreamRxn1Mutations = new HashSet<>();
                        Set<Mutation> upstreamRxn2Mutations = new HashSet<>();

                        //Search D samples for indirect, super-indirect, and super-direct mutations
                        if (D.size() > 0) {
                            targetCooccurringUpstreamRxns.add(upstreamReaction1);
                            targetCooccurringUpstreamRxns.add(upstreamReaction2);

                            targetCooccurringUpstreamReactionFIs.addAll(
                                    reactomeMechismoDataMap.getReactionPatientsFIs(upstreamReaction1, D));
                            targetCooccurringUpstreamReactionFIs.addAll(
                                    reactomeMechismoDataMap.getReactionPatientsFIs(upstreamReaction2, D));

                            upstreamRxn1Mutations.addAll(
                                    reactomeMechismoDataMap.getReactionPatientsMutations(
                                            upstreamReaction1,
                                            D));
                            upstreamRxn2Mutations.addAll(
                                    reactomeMechismoDataMap.getReactionPatientsMutations(
                                            upstreamReaction2,
                                            D));
                            targetRxnmuts.addAll(
                                    reactomeMechismoDataMap.getReactionPatientsMutations(
                                            targetRxn,
                                            D));
                        }

                        targetSamplesWTargetRxnMutations.addAll(
                                reactomeMechismoDataMap.getPatients(targetRxn));

                        //upstream union
                        Set<Mutation> patientRxnMutUnion = new HashSet<>(upstreamRxn1Mutations);
                        patientRxnMutUnion.addAll(upstreamRxn2Mutations);

                        //upstream intersection
                        Set<Mutation> sampleRxnMutIntersection = new HashSet<>(upstreamRxn1Mutations);
                        sampleRxnMutIntersection.retainAll(upstreamRxn2Mutations);

                        //upstream union & target intersection
                        Set<Mutation> sampleRxnMutUnionCpy = new HashSet<>(patientRxnMutUnion);
                        sampleRxnMutUnionCpy.retainAll(targetRxnmuts);
                        targetSuperDirectMutations.addAll(sampleRxnMutUnionCpy);

                        //target - upstream union
                        Set<Mutation> targetRxnmutsCpy = new HashSet<>(targetRxnmuts);
                        targetRxnmutsCpy.removeAll(patientRxnMutUnion);
                        targetDirectMutations.addAll(targetRxnmutsCpy);

                        //upstream intersection - target
                        Set<Mutation> sampleRxnMutIntersectionCpy = new HashSet<>(sampleRxnMutIntersection);
                        sampleRxnMutIntersectionCpy.removeAll(targetRxnmuts);
                        targetSuperIndirectMutations.addAll(sampleRxnMutIntersectionCpy);

                        //upstream union - upstream intersection - target
                        sampleRxnMutUnionCpy = new HashSet<>(patientRxnMutUnion);
                        sampleRxnMutUnionCpy.removeAll(sampleRxnMutIntersection);
                        sampleRxnMutUnionCpy.removeAll(targetRxnmuts);
                        targetIndirectMutations.addAll(sampleRxnMutUnionCpy);

                        //samplesW2 D intersection
                        Set<Patient> samplesW2Dintersection = new HashSet<>(D);
                        samplesW2Dintersection.retainAll(targetSamplesW2MutatedUpstreamRxns);
                        targetSamplesW3plusMutatedUpstreamRxns.addAll(samplesW2Dintersection);

                        targetSamplesW0MutatedUpstreamRxns.addAll(A);
                        targetSamplesW1MutatedUpstreamRxn.addAll(B);
                        targetSamplesW1MutatedUpstreamRxn.addAll(C);
                        targetSamplesW2MutatedUpstreamRxns.addAll(D);//remove samplesW3plus later

                        //FisherExact uses numerical approximation and is sometimes > 1.0000000000
                        targetUpstreamReactionPValues.add(
                                MathUtilities.boundDouble01(
                                        fisherExact.getRightTailedP(
                                                A.size(),
                                                B.size(),
                                                C.size(),
                                                D.size())));
                        List<Double> abcds = new ArrayList<>();
                        abcds.add((double) A.size());
                        abcds.add((double) B.size());
                        abcds.add((double) C.size());
                        abcds.add((double) D.size());
                        targetUpstreamReactionABCDs.add(abcds);
                    }
                }

                if (targetReactionDetected) {

                    if (allTargetRxns.size() != allPValues.size() ||
                            allTargetRxns.size() != allCooccurringUpstreamRxns.size() ||
                            allTargetRxns.size() != allCooccurringUpstreamReactionFIs.size()) {
                        throw new IllegalStateException(String.format("These should all have size == %d",
                                targetUpstreamReactionPValues.size()));
                    }

                    //we lose track of sample support in upstream reaction pair context
                    targetSamplesW0MutatedUpstreamRxns.removeAll(targetSamplesW1MutatedUpstreamRxn);
                    targetSamplesW0MutatedUpstreamRxns.removeAll(targetSamplesW2MutatedUpstreamRxns);
                    targetSamplesW0MutatedUpstreamRxns.removeAll(targetSamplesW3plusMutatedUpstreamRxns);

                    targetSamplesW1MutatedUpstreamRxn.removeAll(targetSamplesW2MutatedUpstreamRxns);
                    targetSamplesW1MutatedUpstreamRxn.removeAll(targetSamplesW3plusMutatedUpstreamRxns);

                    targetSamplesW2MutatedUpstreamRxns.removeAll(targetSamplesW3plusMutatedUpstreamRxns);

                    //we lose track of super relationships in upstream reaction pair context
                    targetIndirectMutations.removeAll(targetSuperIndirectMutations);
                    targetIndirectMutations.removeAll(targetSuperDirectMutations);

                    targetDirectMutations.removeAll(targetSuperDirectMutations);

                    targetSuperIndirectMutations.removeAll(targetSuperDirectMutations);

                    Double combinedPValue;
                    if (targetUpstreamReactionABCDs.size() > 1) {
                        combinedPValue = MathUtilities.calculatePValuesWithEmpiricalBrownsMethod(
                                targetUpstreamReactionABCDs,
                                targetUpstreamReactionPValues);
                    } else {
                        //exactly one upstream reaction pair results in exactly one p-value
                        combinedPValue = targetUpstreamReactionPValues.get(0);
                    }

                    allTargetRxns.add(targetRxn);
                    allCooccurringUpstreamRxns.add(targetCooccurringUpstreamRxns);
                    allCooccurringUpstreamReactionFIs.add(targetCooccurringUpstreamReactionFIs);
                    allSamplesWTargetRxnMutations.add(targetSamplesWTargetRxnMutations);
                    allSamplesW0MutatedUpstreamRxns.add(targetSamplesW0MutatedUpstreamRxns.size());
                    allSamplesW1MutatedUpstreamRxn.add(targetSamplesW1MutatedUpstreamRxn);
                    allSamplesW2MutatedUpstreamRxns.add(targetSamplesW2MutatedUpstreamRxns);
                    allSamplesW3plusMutatedUpstreamRxns.add(targetSamplesW3plusMutatedUpstreamRxns);
                    allSuperIndirectMutatations.add(targetSuperIndirectMutations);
                    allIndirectMutations.add(targetIndirectMutations);
                    allSuperDirectMutations.add(targetSuperDirectMutations);
                    allDirectMutations.add(targetDirectMutations);
                    allPValues.add(MathUtilities.boundDouble01(combinedPValue));
                }
            }

            iterationCounter++;
            if (iterationCounter % 500 == 0) {
                System.out.println(String.format("Processed %d of %d Target Reactions...",
                        iterationCounter, targetReactionCandidates.size()));
                System.out.flush();
            }
        }
        return new CooccurrenceResult(
                allTargetRxns,
                allCooccurringUpstreamRxns,
                allCooccurringUpstreamReactionFIs,
                allSamplesWTargetRxnMutations,
                allSamplesW0MutatedUpstreamRxns,
                allSamplesW1MutatedUpstreamRxn,
                allSamplesW2MutatedUpstreamRxns,
                allSamplesW3plusMutatedUpstreamRxns,
                allSuperIndirectMutatations,
                allIndirectMutations,
                allSuperDirectMutations,
                allDirectMutations,
                allPValues,
                useRxnDist,
                minNUmTargetRxnPatients);
    }

    private Set<Reaction[]> GenerateReactionSetPairs(Set<Reaction> reactionSet) {
        Set<Reaction[]> reactionSetPairs = new HashSet<>();
        List<Reaction> reactionList = new ArrayList<>(reactionSet);
        for (int i = 0; i < reactionList.size() - 1; i++) { //first to second-last
            for (int j = i + 1; j < reactionList.size(); j++) { //second to last
                if (j == i) {
                    continue;
                }
                Reaction[] pair = {reactionList.get(i), reactionList.get(j)};
                reactionSetPairs.add(pair);
            }
        }
        return reactionSetPairs;
    }
}
