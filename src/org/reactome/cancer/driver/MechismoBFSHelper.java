package org.reactome.cancer.driver;

import org.reactome.fi.util.InteractionUtilities;
import org.reactome.r3.graph.BreadthFirstSearch;

import java.util.*;
import java.util.stream.Collectors;

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

    public double calculateMinShortestFIPath(Set<String> patient1FISet,
                                           Set<String> patient2FISet,
                                           Map<String, Integer> pairToDist) {
        int total = 0;
        int count = 0;
        // From set1 to set2
        for (String fi1 : patient1FISet) {
            // Find the shortest path
            int shortest = Integer.MAX_VALUE;
            for (String fi2 : patient2FISet) {
                String pair = trimCapsAlphaSort("\t",fi1,fi2);
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
        for (String fi2 : patient2FISet) {
            // Find the shortest path
            int shortest = Integer.MAX_VALUE;
            for (String fi1 : patient1FISet) {
                String pair = trimCapsAlphaSort("\t",fi1,fi2);
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

    public Map<String, Integer> calculateShortestFIPath(Map<String, Set<String>> patientToFIs,
                                                      Map<String, Set<String>> fiToPartners,
                                                      BreadthFirstSearch bfs) {
        Set<String> allFIs = patientToFIs.values()
                .stream()
                .flatMap(fis -> fis.stream())
                .collect(Collectors.toSet());
        List<String> allFIList = new ArrayList<>(allFIs);
        allFIList.sort(Comparator.naturalOrder());

        Map<String, Integer> pairToDist = new HashMap<>();
        for (int i = 0; i < allFIList.size(); i++) {
            String fi1 = allFIList.get(i);
            Map<String, Integer> fiToDist = bfs.getDistances(fi1,
                    allFIList.subList(i, allFIList.size()), // Don't exclude itself since we need this value too
                    fiToPartners);
            List<String> fiToDistKeysetList = new ArrayList<>(fiToDist.keySet());
            fiToDistKeysetList.sort(Comparator.naturalOrder());
            for (String fi2 : fiToDistKeysetList)
                pairToDist.put(trimCapsAlphaSort("\t",fi1,fi2), fiToDist.get(fi2));
        }
        return pairToDist;
    }

    private String trimCapsAlphaSort(String delim, String s1, String s2){
        List<String> stringList = new ArrayList<>();
        stringList.add(s1.trim().toUpperCase());
        stringList.add(s2.trim().toUpperCase());
        Collections.sort(stringList);
        return String.join(delim,stringList);
    }

}
