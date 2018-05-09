package org.reactome.cancer.driver;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.reactome.fi.util.InteractionUtilities;
import org.reactome.r3.graph.BreadthFirstSearch;

/**
 * This helper class is used to perform BFS related stuff for mechismo data.
 * @author wug
 *
 */
public class MechismoBFSHelper {
    
    public MechismoBFSHelper() {
    }
    
    /**
     * Average values between two sets of reactions. Two steps are used: from set1
     * to set2 and then from set2 to set1 for case like this: one reaction in set1
     * but multiple reactions in set2.
     * @param set1
     * @param set2
     * @param pairToDist
     * @return
     */
    public double calculateMinShortestPath(Set<String> set1,
                                           Set<String> set2,
                                           Map<String, Integer> pairToDist) {
        int total = 0;
        int count = 0;
        // From set1 to set2
        for (String reaction1 : set1) {
            // Find the shortest path
            int shortest = Integer.MAX_VALUE;
            for (String reaction2 : set2) {
                String pair = InteractionUtilities.generateFIFromGene(reaction1,
                        reaction2);
                Integer dist = pairToDist.get(pair);
                if (dist == null) {
                    throw new IllegalStateException(pair + " doesn't have a distance!");
                }
                if (dist < shortest)
                    shortest = dist;
            }
            total += shortest;
            count ++;
        }
        // From set2 to set1
        for (String reaction2 : set2) {
            // Find the shortest path
            int shortest = Integer.MAX_VALUE;
            for (String reaction1 : set1) {
                String pair = InteractionUtilities.generateFIFromGene(reaction1,
                                                                      reaction2);
                Integer dist = pairToDist.get(pair);
                if (dist == null) {
                    throw new IllegalStateException(pair + " doesn't have a distance!");
                }
                if (dist < shortest)
                    shortest = dist;
            }
            total += shortest;
            count ++;
        }
        return (double) total / count;
    }

    /**
     * Calculate all possible pairs of distances for quick performance.
     * @param sampleToRections
     * @param reactionToPartners
     * @param bfs
     * @return
     */
    public Map<String, Integer> calculateShortestPath(Map<String, Set<String>> sampleToRections,
                                                      Map<String, Set<String>> reactionToPartners,
                                                      BreadthFirstSearch bfs) {
        // Collection all reactions 
        Set<String> allReactions = sampleToRections.values()
                                                   .stream()
                                                   .flatMap(reactions -> reactions.stream())
                                                   .collect(Collectors.toSet());
        // Get a sorted list of reactions so that we can get distances among all reactions
        List<String> reactionList = new ArrayList<String>(allReactions);
        reactionList.sort(Comparator.naturalOrder());
        
        Map<String, Integer> pairToDist = new HashMap<>();
        for (int i = 0; i < reactionList.size(); i++) {
            String rxt1 = reactionList.get(i);
            Map<String, Integer> rxtToDist = bfs.getDistances(rxt1, 
                                                              reactionList.subList(i, reactionList.size()), // Don't exclude itself since we need this value too
                                                              reactionToPartners);
            List<String> list = new ArrayList<>(rxtToDist.keySet());
            list.sort(Comparator.naturalOrder());
            for (String rxt2 : list)
                pairToDist.put(rxt1 + "\t" + rxt2, rxtToDist.get(rxt2));
        }
        return pairToDist;
    }

    public double calculateAvgShortestFIPath(Set<String> patient1FISet,
                                             Set<String> patient2FISet,
                                             Map<String, Integer> fiPairToDist) {
        int sumOfShortestPathDistances = 0;
        int numOfShortestPathDistances = 0;
        // From patient1 to patient2
        for (String patient1FI : patient1FISet) {
            // Find the shortest path
            int shortestPathDist = Integer.MAX_VALUE;
            for (String patient2FI : patient2FISet) {
                String fiPair = trimCapsAlphaSort("\t",patient1FI,patient2FI);
                Integer fiDist = fiPairToDist.get(fiPair);
                if (fiDist == null) {
                    throw new IllegalStateException(fiPair + " doesn't have a distance!");
                }
                if (fiDist < shortestPathDist)
                    shortestPathDist = fiDist;
            }
            sumOfShortestPathDistances += shortestPathDist;
            numOfShortestPathDistances ++;
        }
        // From patient2 to patient1
        for (String patient2FI : patient2FISet) {
            // Find the shortest path
            int shortestPathDist = Integer.MAX_VALUE;
            for (String patient1FI : patient1FISet) {
                String fiPair = trimCapsAlphaSort("\t",patient1FI,patient2FI);
                Integer fiDist = fiPairToDist.get(fiPair);
                if (fiDist == null) {
                    throw new IllegalStateException(fiPair + " doesn't have a distance!");
                }
                if (fiDist < shortestPathDist)
                    shortestPathDist = fiDist;
            }
            sumOfShortestPathDistances += shortestPathDist;
            numOfShortestPathDistances ++;
        }
        return (double) sumOfShortestPathDistances / numOfShortestPathDistances;
    }

    public Map<String, Integer> calculateShortestFIPath(Map<String, Set<String>> patientToSigFIsForCancerType,
                                                        Map<String, Set<String>> fiToFIsSharingGene,
                                                        BreadthFirstSearch bfs) {
        Set<String> allSigFIsForCancerType = patientToSigFIsForCancerType.values()
                .stream()
                .flatMap(fis -> fis.stream())
                .collect(Collectors.toSet());
        List<String> allSigFIsForCancerTypeList = new ArrayList<>(allSigFIsForCancerType);
        allSigFIsForCancerTypeList.sort(Comparator.naturalOrder());

        Map<String, Integer> fiPairToDist = new HashMap<>();
        for (int i = 0; i < allSigFIsForCancerTypeList.size(); i++) {
            String fi1 = allSigFIsForCancerTypeList.get(i);
            Map<String, Integer> fiToDistanceFromFI1 = bfs.getDistances(fi1,
                    allSigFIsForCancerTypeList.subList(i,
                            allSigFIsForCancerTypeList.size()), // Don't exclude itself since we need this value too
                    fiToFIsSharingGene);
            List<String> fisWithDistancesToFI1 = new ArrayList<>(fiToDistanceFromFI1.keySet());
            fisWithDistancesToFI1.sort(Comparator.naturalOrder());
            for (String fi2 : fisWithDistancesToFI1)
                fiPairToDist.put(trimCapsAlphaSort("\t",fi1,fi2), fiToDistanceFromFI1.get(fi2));
        }
        return fiPairToDist;
    }

    private String trimCapsAlphaSort(String delim, String s1, String s2){
        List<String> stringList = new ArrayList<>();
        stringList.add(s1.trim().toUpperCase());
        stringList.add(s2.trim().toUpperCase());
        Collections.sort(stringList);
        return String.join(delim,stringList);
    }

}
