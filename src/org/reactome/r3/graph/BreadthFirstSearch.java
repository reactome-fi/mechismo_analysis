/*
 * Created on Jan 10, 2008
 *
 */
package org.reactome.r3.graph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.reactome.r3.util.InteractionUtilities;

/**
 * Breadth first search is used to find shortest path among two proteins in the functional
 * interaction network.
 * @author guanming
 *
 */
public class BreadthFirstSearch {
    private Comparator<Edge> edgeSorter;
    
    public BreadthFirstSearch() {
    }
    
    public void setEdgeSorter(Comparator<Edge> edgeSorter) {
        this.edgeSorter = edgeSorter;
    }
    
    public Double calculateCliqueness(String id,
                                      Map<String, Set<String>> idToPartners) {
        Set<String> partners = idToPartners.get(id);
        List<String> partnerList = new ArrayList<String>(partners);
        int size = partners.size();
        if (size == 1)
            return 1.0;
        int total = size * (size - 1) / 2;
        int c = 0;
        for (int i = 0; i < size - 1; i++) {
            String id1 = partnerList.get(i);
            Set<String> partners1 = idToPartners.get(id1);
            for (int j = i + 1; j < size; j++) {
                String id2 = partnerList.get(j);
                if (partners1.contains(id2))
                    c ++;
            }
        }
        return (double) c / total;
    }
    
    public Map<String, Set<String>> generateIdToPartnersMap(Set<String> interactions) {
        return InteractionUtilities.generateProteinToPartners(interactions);
    }
    
    public Map<String, Integer> getDistances(String anchor,
                                             Collection<String> targets,
                                             Map<String, Set<String>> idToPartners) {
        Map<String, Integer> targetToDist = new HashMap<String, Integer>();
        // Do a Breadth-First-Search
        // Consider anchor is contained by targets
        int dist = -1;
        Set<String> searched = new HashSet<String>();
        Set<String> current = new HashSet<String>();
        //Set<String> partners1 = idToPartners.get(anchor);
        //current.addAll(partners1);
        current.add(anchor);
        Set<String> next = new HashSet<String>();
        while (current.size() > 0) {
            dist ++;
            for (String target : targets) {
                if (current.contains(target)) {
                    targetToDist.put(target, dist);
                }
            }
            if (targetToDist.size() == targets.size())
                break; // All have been found
            for (String currentID : current) {
                Set<String> currentPartners = idToPartners.get(currentID);
                next.addAll(currentPartners);
            }
            searched.addAll(current);
            current.clear();
            next.removeAll(searched);
            current.addAll(next);
            next.clear();
        }
        // Do a sanity check
        if (targetToDist.size() < targets.size()) {
            for (String target : targets) {
                Integer dist1 = targetToDist.get(target);
                if (dist1 == null)
                    targetToDist.put(target, Integer.MAX_VALUE);
            }
        }
        return targetToDist;
    }
    
    public int getDistance(String id1,
                           String id2,
                           Map<String, Set<String>> idToPartners) {
        // Do a Breadth-First-Search
        int dist = 0;
        Set<String> searched = new HashSet<String>();
        Set<String> current = new HashSet<String>();
        Set<String> partners1 = idToPartners.get(id1);
        current.addAll(partners1);
        Set<String> next = new HashSet<String>();
        boolean isFound = false;
        while (current.size() > 0) {
            dist ++;
            if (current.contains(id2)) {
                isFound = true;
                break;
            }
            for (String currentID : current) {
                Set<String> currentPartners = idToPartners.get(currentID);
                next.addAll(currentPartners);
            }
            searched.addAll(current);
            current.clear();
            next.removeAll(searched);
            current.addAll(next);
            next.clear();
        }
        if (!isFound)
            return Integer.MAX_VALUE;
        return dist;
    }
    
    /**
     * Generate a shortest path between two provided ids.
     * @param id1
     * @param id2
     * @param idToPartners
     * @return
     */
    public List<String> getShortestPath(String id1,
                                        String id2,
                                        Map<String, Set<String>> idToPartners) {
        Map<TreeNode, List<Edge>> nodeToEdges = initGraph(idToPartners);
        List<String> path = generateShortestPath(id1, id2, nodeToEdges);
        return path;
    }

    /**
     * Generate a shortest path between two provided ids.
     * @param id1
     * @param id2
     * @param nodeToEdges
     * @return
     */
    public List<String> generateShortestPath(String id1,
                                             String id2,
                                             Map<TreeNode, List<Edge>> nodeToEdges) {
        // Get the tree node
        TreeNode start = null;
        TreeNode end = null;
        for (TreeNode node : nodeToEdges.keySet()) {
            if (node.id.equals(id1))
                start = node;
            else if (node.id.equals(id2))
                end = node;
        }
        buildBFSTree(start, nodeToEdges, null);
        return readShortestPath(end);
    }

    /**
     * A helper method to read a shortest path from a spanning tree.
     * @param end
     * @return
     */
    private List<String> readShortestPath(TreeNode end) {
        // Get the path
        List<String> path = new ArrayList<String>();
        path.add(end.id);
        TreeNode parent = end.parent;
        while (parent != null) {
            path.add(parent.id);
            parent = parent.parent;
        }
        return path;
    }
    
    
    /**
     * A faster method to generate a pair-wise shortest path among a list of ids. This method
     * used one spanning tree for one starting id, so it should be faster by avoiding creating
     * multiple spanning tree for the same starting id.
     * @param ids
     * @param nodeToEdges
     * @return
     */
    public Map<String, List<String>> generateShortestPath(List<String> ids,
                                                          Map<TreeNode, List<Edge>> nodeToEdges) {
        // For fast search
        Map<String, TreeNode> idToNode = new HashMap<String, TreeNode>();
        for (String id : ids) {
            for (TreeNode node : nodeToEdges.keySet()) {
                if (node.id.equals(id)) {
                    idToNode.put(id, node);
                    break;
                }
            }
        }
        // Do a batch search so that only one tree is needed to be constructed for one start node
        Map<String, List<String>> pairToPath = new HashMap<String, List<String>>();
        for (int i = 0; i < ids.size() - 1; i++) {
            String startId = ids.get(i);
            TreeNode startNode = idToNode.get(startId);
            buildBFSTree(startNode, 
                         nodeToEdges, 
                         edgeSorter);
            for (int j = i + 1; j < ids.size(); j++) {
                String endId = ids.get(j);
                TreeNode endNode = idToNode.get(endId);
                List<String> path = readShortestPath(endNode);
                int compare = startId.compareTo(endId);
                if (compare < 0)
                    pairToPath.put(startId + "\t" + endId,
                                   path);
                else
                    pairToPath.put(endId + "\t" + startId,
                                   path);
            }
        }
        return pairToPath;
    }
    
    /**
     * This method should not used any more since it cannot find a minimum span of a set of
     * genes. Use methods in CancerResequenceDataSetAnalyzer.
     * @param ids
     * @param nodeToEdges
     * @return
     */
    @Deprecated
    public Set<String> getMiniSpanFIsFromGraph(Collection<String> ids,
                                               Map<TreeNode, List<Edge>> nodeToEdges) {
        // Find a Start node to build the spanning tree
        final List<TreeNode> priorityNodes = new ArrayList<TreeNode>();
        for (TreeNode node : nodeToEdges.keySet()) {
            if (ids.contains(node.id)) {
                priorityNodes.add(node);
            }
        }
        TreeNode start = priorityNodes.get(0);
        //priorityNodes.remove(start);
        Comparator<Edge> edgeSorter = new Comparator<Edge>() {
            public int compare(Edge edge1, Edge edge2) {
                TreeNode node1 = edge1.node1;
                if (node1.label != null)
                    node1 = edge1.node2;
                TreeNode node2 = edge2.node1;
                if (node2.label != null)
                    node2 = edge2.node2;
                if (priorityNodes.contains(node1) &&
                    !priorityNodes.contains(node2))
                    return -1;
                if (!priorityNodes.contains(node1) &&
                    priorityNodes.contains(node2))
                    return 1;
                return 0;
            }
        };
        buildBFSTree(start,
                     nodeToEdges, 
                     edgeSorter);
        // Build the returned FIs
        Set<String> fis = new HashSet<String>();
        for (TreeNode node : priorityNodes) {
            List<String> path = buildPath(node);
            fis.addAll(path);
        }
        return fis;
    }
    
    private List<String> buildPath(TreeNode node) {
        List<String> ids = new ArrayList<String>();
        ids.add(node.id);
        TreeNode parent = node.parent;
        while (parent != null) {
            ids.add(parent.id);
            parent = parent.parent;
        }
        // Convert to FIs
        List<String> fis = new ArrayList<String>();
        if (ids.size() == 1)
            return fis;
        for (int i = 0; i < ids.size() - 1; i++) {
            String id1 = ids.get(i);
            String id2 = ids.get(i + 1);
            fis.add(id1 + "\t" + id2);
        }
        return fis;
    }
    
//    @Test
//    public void testShortestPath() throws IOException {
//        String fileName = R3Constants.RESULT_DIR + "FI73InGene_061008_BigComp.txt";
//        FileUtility fu = new FileUtility();
//        Set<String> fis = fu.loadInteractions(fileName);
//        Map<String, Set<String>> idToPartners = generateIdToPartnersMap(fis);
//        //String id1 = "GNL3L";
//        //String id2 = "PMCH";
//        String id1 = "PIK3R1";
//        String id2 = "CDKN2A";
//        List<String> path = getShortestPath(id1, 
//                                            id2, 
//                                            idToPartners);
//        System.out.println(path + ": " + (path.size() - 1));
//        int length = getDistance(id1, id2, idToPartners);
//        System.out.println("From distance: " + length);
//    }
    
    private Map<TreeNode, List<Edge>> initGraph(Map<String, Set<String>> idToPartners) {
        Map<String, TreeNode> idToNode = new HashMap<String, TreeNode>();
        for (String id : idToPartners.keySet()) {
            TreeNode node = new TreeNode();
            node.id = id;
            idToNode.put(id, node);
        }
        List<Edge> edgeList = new ArrayList<Edge>();
        for (String id : idToPartners.keySet()) {
            Set<String> partners = idToPartners.get(id);
            TreeNode node1 = idToNode.get(id);
            for (String partner : partners) {
                TreeNode node2 = idToNode.get(partner);
                Edge edge = new Edge();
                edge.node1 = node1;
                edge.node2 = node2;
                edgeList.add(edge);
            }
        }
        Map<TreeNode, List<Edge>> nodeToEdges = new HashMap<TreeNode, List<Edge>>();
        for (Edge edge : edgeList) {
            List<Edge> list1 = nodeToEdges.get(edge.node1);
            if (list1 == null) {
                list1 = new ArrayList<Edge>();
                nodeToEdges.put(edge.node1, list1);
            }
            list1.add(edge);
            List<Edge> list2 = nodeToEdges.get(edge.node2);
            if (list2 == null) {
                list2 = new ArrayList<Edge>();
                nodeToEdges.put(edge.node2, list2);
            }
            list2.add(edge);
        }
        return nodeToEdges;
    }
    
    public Map<TreeNode, List<Edge>> initGraph(Set<String> fis) {
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(fis);
        Map<String, TreeNode> idToNode = new HashMap<String, TreeNode>();
        for (String id : ids) {
            TreeNode node = new TreeNode();
            node.id = id;
            idToNode.put(id, node);
        }
        List<Edge> edgeList = new ArrayList<Edge>();
        int index = 0;
        for (String fi : fis) {
            index = fi.indexOf("\t");
            String id1 = fi.substring(0, index);
            String id2 = fi.substring(index + 1);
            TreeNode node1 = idToNode.get(id1);
            TreeNode node2 = idToNode.get(id2);
            Edge edge = new Edge();
            edge.node1 = node1;
            edge.node2 = node2;
            edgeList.add(edge);
        }
        Map<TreeNode, List<Edge>> nodeToEdges = new HashMap<TreeNode, List<Edge>>();
        for (Edge edge : edgeList) {
            List<Edge> list1 = nodeToEdges.get(edge.node1);
            if (list1 == null) {
                list1 = new ArrayList<Edge>();
                nodeToEdges.put(edge.node1, list1);
            }
            list1.add(edge);
            List<Edge> list2 = nodeToEdges.get(edge.node2);
            if (list2 == null) {
                list2 = new ArrayList<Edge>();
                nodeToEdges.put(edge.node2, list2);
            }
            list2.add(edge);
        }
        return nodeToEdges;
    }
    
    private void buildBFSTree(TreeNode startNode,
                              Map<TreeNode, List<Edge>> nodeToEdges,
                              Comparator<Edge> edgeSorter) { // To sort edges for next
        // Need to reset the cached tree
        for (TreeNode node : nodeToEdges.keySet())
            node.reset();
        int label = 0;
        startNode.label = label;
        List<Edge> current = new ArrayList<Edge>();
        List<Edge> edgeList = nodeToEdges.get(startNode);
        current.addAll(edgeList);
        List<Edge> next = new ArrayList<Edge>();
        //int cycle = 0;
        List<TreeNode> checkingNodes = new ArrayList<TreeNode>();
        // Need to sort edge list too
        if (edgeSorter != null)
            Collections.sort(current, edgeSorter);
        while (current.size() > 0) {
            // Split the current into two parts: those in the priority nodes and those not.
            // TreeNodes in the priority nodes should be searched first 
            checkingNodes.clear();
            // Split the checking into two steps to avoid double counting of edges
            for (Edge edge : current) {
                // Do label first to avoid double counting
                if (edge.node1.label == null) {
                    edge.node1.label = label ++;
                    checkingNodes.add(edge.node1);
                    edge.node2.addChildNode(edge.node1);
                    edge.node1.parent = edge.node2;
                }
                else if (edge.node2.label == null) {
                    edge.node2.label = label ++;
                    checkingNodes.add(edge.node2);
                    edge.node1.addChildNode(edge.node2);
                    edge.node2.parent = edge.node1;
                }
            }
            // Want to get all the frontier edges
            for (TreeNode checkingNode : checkingNodes) {
                List<Edge> edges = nodeToEdges.get(checkingNode);
                for (Edge tmp : edges) {
                    if (tmp.node1.label != null && tmp.node2.label != null)
                        continue;
                    next.add(tmp);
                }
            }
            // Sorting edge to get the priority list.
            if (edgeSorter != null)
                Collections.sort(next, edgeSorter);
            current.clear();
            current.addAll(next);
            next.clear();
            //System.out.println(cycle +  " Current size: " + current.size());
            //cycle ++;
        }
    }
    
    /**
     * Used as spanning tree node.
     * @author wgm
     *
     */
    public class TreeNode {
        private Integer label;
        private String id;
        private TreeNode parent;
        private List<TreeNode> children;
        
        public TreeNode() {
            
        }
        
        public void addChildNode(TreeNode node) {
            if (children == null)
                children = new ArrayList<TreeNode>();
            children.add(node);
        }
        
        public void reset() {
            label = null;
            parent = null;
            if (children != null)
                children.clear();
        }
        
        public String getId() {
            return this.id;
        }
        
        public Integer getLabel() {
            return this.label;
        }
    }
    
    public class Edge {
        private TreeNode node1;
        private TreeNode node2;
        
        public Edge() {
        }
        
        public TreeNode getNode1() {
            return this.node1;
        }
        
        public TreeNode getNode2() {
            return this.node2;
        }
    }

}
