package org.reactome.cancer;

import org.junit.Assert;
import org.junit.Test;

import java.util.*;
import java.util.stream.Collectors;

public class OptimizedBFS {
    private final double DELTA = 0.00001d;
    private Set<String> testSet1;
    private Set<String> testSet2;
    private Map<String, Set<String>> testIdToPartners;

    public double calculateMinShortestPath(Set<String> set1,
                                           Set<String> set2,
                                           Map<String, Set<String>> idToPartners) {
        int total = 0;
        int count = 0;
        // From list1 to list2
        for (String gene1 : set1) {
            Map<String, Integer> geneToDistance = getDistances(gene1,
                    set2,
                    idToPartners);
            // Find the shortest path
            total += geneToDistance.values().stream().min(Integer::compare).orElse(Integer.MAX_VALUE);
            count++;
        }
        // From list2 to list1
        for (String gene2 : set2) {
            Map<String, Integer> geneToDistance = getDistances(gene2,
                    set1,
                    idToPartners);
            // Find the shortest path
            total += geneToDistance.values().stream().min(Integer::compare).orElse(Integer.MAX_VALUE);
            count++;
        }
        return count > 0
                ? (double) total / (double) count
                : -1.0d;
    }

    public double calculateAverageDistance(Set<String> set1,
                                           Set<String> set2,
                                           Map<String, Set<String>> idToPartners) {
        int total = 0;
        int count = 0;
        for (String gene1 : set1) {
            Map<String, Integer> geneToDistances = getDistances(gene1,
                    set2,
                    idToPartners);
            count += geneToDistances.size();
            total += geneToDistances.values().stream().mapToInt(Integer::intValue).sum();
        }
        return count > 0
                ? (double) total / (double) count
                : -1.0d;
    }

    private Map<String, Integer> getDistances(String anchor,
                                             Set<String> targets,
                                             Map<String, Set<String>> idToPartners) {
        Map<String, Integer> targetToDist = new HashMap<>();
        // Do a Breadth-First-Search
        // Consider anchor is contained by targets
        int dist = -1;
        Set<String> searched = new HashSet<>();
        Set<String> current = new HashSet<>();
        current.add(anchor);
        Set<String> next = new HashSet<>();
        while (!current.isEmpty()) {
            dist++;
            for (String target : targets.stream().filter(current::contains).collect(Collectors.toList())) {
                targetToDist.put(target, dist);
            }
            if (targetToDist.size() == targets.size())
                break; // All have been found
            for (String currentID : current) {
                next.addAll(idToPartners.get(currentID));
            }
            searched.addAll(current);
            current.clear();
            next.removeAll(searched);
            current.addAll(next);
            next.clear();
        }
        // Do a sanity check
        if (targetToDist.size() < targets.size()) {
            for (String target : targets.stream().filter(target -> targetToDist.get(target) != null).collect(Collectors.toList())) {
                targetToDist.put(target, Integer.MAX_VALUE);
            }
        }
        return targetToDist;
    }

    @Test
    public void TestcalculateMinShortestPath() {
        setUp();
        Assert.assertEquals("minimum distance calculation incorrect",
                2.5d,
                calculateMinShortestPath(this.testSet1, this.testSet2, this.testIdToPartners),
                DELTA);
    }

    @Test
    public void TestcalculateAverageDistance() {
        setUp();
        Assert.assertEquals("average distance calculation incorrect",
                3.0d,
                calculateAverageDistance(this.testSet1, this.testSet2, this.testIdToPartners),
                DELTA);
    }

    @Test
    public void TestgetDistances() {
        setUp();
        Map<String, Integer> expected1a = new HashMap<>();
        expected1a.put("b", 1);
        expected1a.put("i", 2);
        Map<String, Integer> result1a = getDistances("a", testSet1, testIdToPartners);
        Assert.assertEquals(expected1a, result1a);

        Map<String, Integer> expected1b = new HashMap<>();
        expected1b.put("b", 0);
        expected1b.put("i", 3);
        Map<String, Integer> result1b = getDistances("b", testSet1, testIdToPartners);
        Assert.assertEquals(expected1b, result1b);

        Map<String, Integer> expected1c = new HashMap<>();
        expected1c.put("b", 2);
        expected1c.put("i", 1);
        Map<String, Integer> result1c = getDistances("c", testSet1, testIdToPartners);
        Assert.assertEquals(expected1c, result1c);

        Map<String, Integer> expected1d = new HashMap<>();
        expected1d.put("b", 2);
        expected1d.put("i", 3);
        Map<String, Integer> result1d = getDistances("d", testSet1, testIdToPartners);
        Assert.assertEquals(expected1d, result1d);

        Map<String, Integer> expected1e = new HashMap<>();
        expected1e.put("b", 2);
        expected1e.put("i", 3);
        Map<String, Integer> result1e = getDistances("e", testSet1, testIdToPartners);
        Assert.assertEquals(expected1e, result1e);

        Map<String, Integer> expected1f = new HashMap<>();
        expected1f.put("b", 3);
        expected1f.put("i", 4);
        Map<String, Integer> result1f = getDistances("f", testSet1, testIdToPartners);
        Assert.assertEquals(expected1f, result1f);

        Map<String, Integer> expected1g = new HashMap<>();
        expected1g.put("b", 3);
        expected1g.put("i", 4);
        Map<String, Integer> result1g = getDistances("g", testSet1, testIdToPartners);
        Assert.assertEquals(expected1g, result1g);

        Map<String, Integer> expected1h = new HashMap<>();
        expected1h.put("b", 1);
        expected1h.put("i", 4);
        Map<String, Integer> result1h = getDistances("h", testSet1, testIdToPartners);
        Assert.assertEquals(expected1h, result1h);

        Map<String, Integer> expected1i = new HashMap<>();
        expected1i.put("b", 3);
        expected1i.put("i", 0);
        Map<String, Integer> result1i = getDistances("i", testSet1, testIdToPartners);
        Assert.assertEquals(expected1i, result1i);

        Map<String, Integer> expected2a = new HashMap<>();
        expected2a.put("e", 1);
        expected2a.put("f", 2);
        Map<String, Integer> result2a = getDistances("a", testSet2, testIdToPartners);
        Assert.assertEquals(expected2a, result2a);

        Map<String, Integer> expected2b = new HashMap<>();
        expected2b.put("e", 2);
        expected2b.put("f", 3);
        Map<String, Integer> result2b = getDistances("b", testSet2, testIdToPartners);
        Assert.assertEquals(expected2b, result2b);

        Map<String, Integer> expected2c = new HashMap<>();
        expected2c.put("e", 2);
        expected2c.put("f", 3);
        Map<String, Integer> result2c = getDistances("c", testSet2, testIdToPartners);
        Assert.assertEquals(expected2c, result2c);

        Map<String, Integer> expected2d = new HashMap<>();
        expected2d.put("e", 2);
        expected2d.put("f", 1);
        Map<String, Integer> result2d = getDistances("d", testSet2, testIdToPartners);
        Assert.assertEquals(expected2d, result2d);

        Map<String, Integer> expected2e = new HashMap<>();
        expected2e.put("e", 0);
        expected2e.put("f", 3);
        Map<String, Integer> result2e = getDistances("e", testSet2, testIdToPartners);
        Assert.assertEquals(expected2e, result2e);

        Map<String, Integer> expected2f = new HashMap<>();
        expected2f.put("e", 3);
        expected2f.put("f", 0);
        Map<String, Integer> result2f = getDistances("f", testSet2, testIdToPartners);
        Assert.assertEquals(expected2f, result2f);

        Map<String, Integer> expected2g = new HashMap<>();
        expected2g.put("e", 1);
        expected2g.put("f", 4);
        Map<String, Integer> result2g = getDistances("g", testSet2, testIdToPartners);
        Assert.assertEquals(expected2g, result2g);

        Map<String, Integer> expected2h = new HashMap<>();
        expected2h.put("e", 3);
        expected2h.put("f", 4);
        Map<String, Integer> result2h = getDistances("h", testSet2, testIdToPartners);
        Assert.assertEquals(expected2h, result2h);

        Map<String, Integer> expected2i = new HashMap<>();
        expected2i.put("e", 3);
        expected2i.put("f", 4);
        Map<String, Integer> result2i = getDistances("i", testSet2, testIdToPartners);
        Assert.assertEquals(expected2i, result2i);
    }

    public void setUp() {

        /*

        f-d   b-h
           \ /
            a
           / \
        g-e   c-i

        */

        Set<String> tmp;
        this.testSet1 = new HashSet<>();
        this.testSet2 = new HashSet<>();
        this.testIdToPartners = new HashMap<>();

        testSet1.add("b");
        testSet1.add("i");

        testSet2.add("e");
        testSet2.add("f");

        tmp = new HashSet<>();
        tmp.add("b");
        tmp.add("c");
        tmp.add("d");
        tmp.add("e");
        testIdToPartners.put("a", tmp);

        tmp = new HashSet<>();
        tmp.add("a");
        tmp.add("h");
        testIdToPartners.put("b", tmp);

        tmp = new HashSet<>();
        tmp.add("a");
        tmp.add("i");
        testIdToPartners.put("c", tmp);

        tmp = new HashSet<>();
        tmp.add("a");
        tmp.add("f");
        testIdToPartners.put("d", tmp);

        tmp = new HashSet<>();
        tmp.add("a");
        tmp.add("g");
        testIdToPartners.put("e", tmp);

        tmp = new HashSet<>();
        tmp.add("d");
        testIdToPartners.put("f", tmp);

        tmp = new HashSet<>();
        tmp.add("e");
        testIdToPartners.put("g", tmp);

        tmp = new HashSet<>();
        tmp.add("b");
        testIdToPartners.put("h", tmp);

        tmp = new HashSet<>();
        tmp.add("c");
        testIdToPartners.put("i", tmp);

    }
}
