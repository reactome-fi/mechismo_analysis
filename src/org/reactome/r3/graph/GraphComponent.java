/*
 * Created on Jun 9, 2008
 *
 */
package org.reactome.r3.graph;

import java.util.HashSet;
import java.util.Set;

import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleDirectedGraph;

/**
 * This class is used to describe a graph component.
 * @author wgm
 *
 */
public class GraphComponent {
    // score of this GraphComponent
    private Double score;
    // The underneath DiGraph to keep the data structure. The reason
    // why a digrah is used is that we can keep track of adding sequence:
    // which one is added after which
    private DirectedGraph<String, DefaultEdge> component;
    // Used to calculate score
    private ScoreCalculator scoreCalculator;
    // Id of this component
    private int id;
    
    public GraphComponent() {
        component = new SimpleDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
    }
    
    public void setId(int id) {
        this.id = id;
    }
    
    public int getId() {
        return this.id;
    }
    
    public void setScoreCalculator(ScoreCalculator calculator) {
        this.scoreCalculator = calculator;
    }
    
    /**
     * Add a GraphNode to this GraphComponent.
     * @param node
     */
    public void addNode(String node) {
        component.addVertex(node);
        // Score is invalidate
        resetScore();
    }
    
    public void addEdge(String source, String target) {
        component.addEdge(source, target);
        resetScore();
    }
    
    /**
     * Remove a node from this GraphComponent.
     * @param node
     */
    public void removeNode(String node) {
        component.removeVertex(node);
        resetScore();
    }
    
    /**
     * Reset the score of this GraphComponent so that a new score can be calculated.
     */
    public void resetScore() {
        score = null;
    }
    
    /**
     * Check if this GraphComponent contains the passed GraphNode.
     * @param node
     * @return
     */
    public boolean containsNode(String node) {
        return component.containsVertex(node);
    }
    
    /**
     * Get the score of this GraphComponent.
     * @return
     */
    public double getScore() {
        if (score != null)
            return score;
        if (scoreCalculator == null)
            throw new IllegalStateException("GraphComponent.getScore(): no ScoreCalculator is specified.");
        if (score == null)
            calculateScore();
        return score;
    }
    
    public void setScore(Double score) {
        this.score = score;
    }
    
    /**
     * This method is used to calculate score of this GraphComponent.
     */
    private void calculateScore() {
        score = scoreCalculator.calculateScore(this);
    }
    
    /**
     * So called removable nodes are leaf nodes that can be removed without
     * effect other nodes in the component.
     * @return
     */
    public Set<String> getRemovableNodes() {
        Set<String> removable = new HashSet<String>();
        Set<String> nodes = component.vertexSet();
        for (String node : nodes) {
            if (component.outDegreeOf(node) == 0)
                removable.add(node);
        }
        return removable;
    }
    
    /**
     * This method returns the source node attached to a removable node.
     * @param node
     * @return
     */
    public String getSourceForRemovableNode(String node) {
        Set<DefaultEdge> inEdges = component.incomingEdgesOf(node);
        // Pick the first edge only
        if (inEdges == null || inEdges.size() == 0)
            return null;
        DefaultEdge edge = inEdges.iterator().next();
        return component.getEdgeSource(edge);
    }
    
    /**
     * Return all contained nodes by this GraphComponent.
     * @return
     */
    public Set<String> getAllNodes() {
        return component.vertexSet();
    }
    
    /**
     * This interface defines a method to calculate the score for a GraphComponent.
     * @author wgm
     *
     */
    public static interface ScoreCalculator {
        
        public double calculateScore(GraphComponent comp);
    }
    
}
