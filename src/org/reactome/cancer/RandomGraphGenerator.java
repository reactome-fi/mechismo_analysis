package org.reactome.cancer;

import org.jgrapht.Graph;
import org.jgrapht.Graphs;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;

import java.util.*;

public class RandomGraphGenerator {
    private final Long randomSeed = 88L;
    private DefaultDirectedGraph<Long, DefaultEdge> reactionGraph;
    private Random prng;

    public RandomGraphGenerator(DefaultDirectedGraph<Long, DefaultEdge> reactionGraph) {
        this.reactionGraph = reactionGraph;
        this.prng = new Random(this.randomSeed);
    }

    public Graph<Long, DefaultEdge> GenerateRandomGraph(boolean rewireLargestComponentOnly) {
        DefaultDirectedGraph<Long, DefaultEdge> rewiredReactionGraph =
                new DefaultDirectedGraph<>(DefaultEdge.class);
        Graphs.addGraph(rewiredReactionGraph, this.reactionGraph);

        RewireComponents(rewiredReactionGraph,rewireLargestComponentOnly);

        return rewiredReactionGraph;
    }

    private void RewireComponents(DefaultDirectedGraph<Long, DefaultEdge> rewiredReactionGraph,
                                  boolean rewireLargestComponentOnly) {
        ConnectivityInspector connectivityInspector = new ConnectivityInspector(rewiredReactionGraph);
        List<Set<Long>> connectedSets = connectivityInspector.connectedSets();

        int initialComponentCount = connectedSets.size();
        int initialEdgeCount = rewiredReactionGraph.edgeSet().size();
        int initialVtxCount = rewiredReactionGraph.vertexSet().size();

        Set<Long> largestComponent = new HashSet<>();
        if(rewireLargestComponentOnly){
           for(Set<Long> connectedSet : connectedSets){
               if(connectedSet.size() > largestComponent.size()){
                   largestComponent = connectedSet;
               }
           }
        }

        for (Set<Long> connectedSet : connectedSets) {
            if(rewireLargestComponentOnly &&
                    !connectedSet.equals(largestComponent)){
                continue;
            }
            if (connectedSet.size() >= 3) {

                if (connectedSet.size() >= initialVtxCount) {
                    throw new IllegalStateException(String.format(
                            "component vtx count (%d) should be less than initial (%d)",
                            connectedSet.size(),
                            initialVtxCount));
                }

                Set<DefaultEdge> componentEdges = new HashSet<>();
                for (Long vertex : connectedSet) {
                    componentEdges.addAll(rewiredReactionGraph.incomingEdgesOf(vertex));
                    componentEdges.addAll(rewiredReactionGraph.outgoingEdgesOf(vertex));
                }

                if (componentEdges.size() >= initialEdgeCount) {
                    throw new IllegalStateException(String.format(
                            "component edge count (%d) should be less than initial (%d)",
                            componentEdges.size(),
                            initialEdgeCount));
                }


                int E = componentEdges.size();
                int V = connectedSet.size();
                DefaultEdge randEdge1, randEdge2;
                ArrayList<DefaultEdge> componentEdgesAry = new ArrayList<>(componentEdges);
                ArrayList<Long> componentVtxAry = new ArrayList<>(connectedSet);
                Long sourceVtx1, sourceVtx2, targetVtx1, targetVtx2;

                //from https://en.wikipedia.org/wiki/Degree-preserving_randomization
                int rewires = (int) ((E / 2.0) * Math.log(1 / Math.pow(10, -7)));
                for (int i = 0; i < rewires; i++) {

                    randEdge1 = componentEdgesAry.get(this.prng.nextInt(E));
                    randEdge2 = componentEdgesAry.get(this.prng.nextInt(E));

                    while (Objects.equals(randEdge1, randEdge2)) {
                        randEdge1 = componentEdgesAry.get(this.prng.nextInt(E));
                        randEdge2 = componentEdgesAry.get(this.prng.nextInt(E));
                    }

                    rewiredReactionGraph.removeEdge(randEdge1);
                    componentEdgesAry.remove(randEdge1);

                    int curEdgeCount1 = rewiredReactionGraph.edgeSet().size();
                    if (curEdgeCount1 != (initialEdgeCount - 1)) {
                        throw new IllegalStateException(String.format(
                                "current edge count (%d) should not differ from initial - 1 (%d)",
                                curEdgeCount1,
                                initialEdgeCount - 1));
                    }

                    rewiredReactionGraph.removeEdge(randEdge2);
                    componentEdgesAry.remove(randEdge2);

                    int curEdgeCount2 = rewiredReactionGraph.edgeSet().size();
                    if (curEdgeCount2 != (initialEdgeCount - 2)) {
                        throw new IllegalStateException(String.format(
                                "current edge count (%d) should not differ from initial - 2 (%d)",
                                curEdgeCount2,
                                initialEdgeCount - 2));
                    }

                    sourceVtx1 = componentVtxAry.get(this.prng.nextInt(V));
                    targetVtx1 = componentVtxAry.get(this.prng.nextInt(V));

                    sourceVtx2 = componentVtxAry.get(this.prng.nextInt(V));
                    targetVtx2 = componentVtxAry.get(this.prng.nextInt(V));

                    while (Objects.equals(sourceVtx1, sourceVtx2) ||
                            Objects.equals(targetVtx1, targetVtx2) ||
                            Objects.equals(sourceVtx1, targetVtx1) ||
                            Objects.equals(sourceVtx2, targetVtx2)) {
                        sourceVtx1 = componentVtxAry.get(this.prng.nextInt(V));
                        targetVtx1 = componentVtxAry.get(this.prng.nextInt(V));

                        sourceVtx2 = componentVtxAry.get(this.prng.nextInt(V));
                        targetVtx2 = componentVtxAry.get(this.prng.nextInt(V));
                    }

                    //randomly flip edge direction
                    //rewire edge 1
                    WireEdge(V,
                            sourceVtx1,
                            targetVtx2,
                            rewiredReactionGraph,
                            componentEdgesAry,
                            componentVtxAry);

                    int curEdgeCount3 = rewiredReactionGraph.edgeSet().size();
                    if (curEdgeCount3 != (initialEdgeCount - 1)) {
                        throw new IllegalStateException(String.format(
                                "current edge count (%d) should not differ from initial - 1 (%d)",
                                curEdgeCount3,
                                initialEdgeCount - 1));
                    }

                    //rewire edge 2
                    WireEdge(V,
                            sourceVtx2,
                            targetVtx1,
                            rewiredReactionGraph,
                            componentEdgesAry,
                            componentVtxAry);

                    int curEdgeCount4 = rewiredReactionGraph.edgeSet().size();
                    if (curEdgeCount4 != initialEdgeCount) {
                        throw new IllegalStateException(String.format(
                                "current edge count (%d) should not differ from initial (%d)",
                                curEdgeCount4,
                                initialEdgeCount));
                    }
                }
            }
        }
        int finalComponentCount = connectivityInspector.connectedSets().size();
        int finalEdgeCount = rewiredReactionGraph.edgeSet().size();
        int finalVtxCount = rewiredReactionGraph.vertexSet().size();

        if (initialComponentCount > finalComponentCount) {
            throw new IllegalStateException(String.format(
                    "Final component count (%d) should be no greater than initial (%d)",
                    finalComponentCount,
                    initialComponentCount));
        }

        if (initialEdgeCount != finalEdgeCount) {
            throw new IllegalStateException(String.format(
                    "Final edge count (%d) should not differ from initial (%d)",
                    finalEdgeCount,
                    initialEdgeCount));
        }

        if (initialVtxCount != finalVtxCount) {
            throw new IllegalStateException(String.format(
                    "Final vertex count (%d) should not differ from initial (%d)",
                    finalVtxCount,
                    initialVtxCount));
        }
    }

    private void WireEdge(int V,
                          Long sourceVtx,
                          Long targetVtx,
                          DefaultDirectedGraph<Long, DefaultEdge> graph,
                          ArrayList<DefaultEdge> edgeAry,
                          ArrayList<Long> vtxAry) {
        int rand;
        while (graph.getEdge(sourceVtx, targetVtx) != null) {
            rand = this.prng.nextInt(V);
            sourceVtx = vtxAry.get(rand);
            rand = this.prng.nextInt(V);
            targetVtx = vtxAry.get(rand);
        }
        graph.addEdge(sourceVtx, targetVtx);
        edgeAry.add(
                graph.getEdge(sourceVtx, targetVtx));
    }
}
