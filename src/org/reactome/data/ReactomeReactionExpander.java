/*
 * Created on Jun 8, 2017
 *
 */
package org.reactome.data;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.jgraph.graph.DefaultEdge;
import org.jgrapht.DirectedGraph;
import org.jgrapht.GraphPath;
import org.jgrapht.alg.AllDirectedPaths;
import org.jgrapht.alg.CycleDetector;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.junit.Test;
import org.reactome.fi.util.ReactomeUtilities;

/**
 * If a reaction uses one or more EntitySet instances as its participants, reactions can be expanded
 * into multiple reactions by permuting members in EntitySet. In this implementation, a reaction instance
 * is converted into a DirectedGraph object using jGraphtT library, and then existing algorithm in the library
 * is used to find paths from inputs to outputs. To avoid cycles, inputs are removed in outputs if they are
 * used. 
 * @author gwu
 *
 */
@SuppressWarnings("unchecked")
public class ReactomeReactionExpander {
    private static final Logger logger = Logger.getLogger(ReactomeReactionExpander.class);
    
    /**
     * Default constructor.
     */
    public ReactomeReactionExpander() {
    }
    
    @Test
    public void testExpandReaction() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_59_plus_i",
                                            "root",
                                            "macmysql01");
        Long dbId = 199443L;
        dbId = 2316434L;
        GKInstance reaction = dba.fetchInstance(dbId);
        
        // Expand the passed reaction into a list of lists of instances.
        List<List<GKInstance>> listOfInstances = expandReaction(reaction);
        for (int i = 0; i < listOfInstances.size(); i++) {
            List<GKInstance> instances = listOfInstances.get(i);
            System.out.println("List " + i + ": " + instances.size());
            instances.stream().forEach(System.out::println);
        }
        
        // Expand the test reaction into a list of sets of genes.
        Set<Set<String>> listOfGenes = extractGenesFromReaction(reaction);
        System.out.println("\nExtracted as genes:");
        int c = 0;
        for (Set<String> genes : listOfGenes) {
            System.out.println("Genes " + c++ + ": " + genes.size());
            System.out.println(String.join(", ", genes));
        }
    }
    
    /**
     * This method expands the passed reaction into multiple implementation if EntitySets are annotated
     * in any layer of reaction participants. Each set of genes in the returned List object is extracted
     * from such an implementation of the passed reaction. Two different lists of EWASes may generate the
     * same set of genes (e.g. different phosphorylation sites). To avoid case like this, Set is used.
     * @param reaction
     * @return
     * @throws Exception
     */
    public Set<Set<String>> extractGenesFromReaction(GKInstance reaction) throws Exception {
        List<List<GKInstance>> listOfInstances = expandReaction(reaction);
        Set<Set<String>> setOfGenes = listOfInstances.stream().map(instances -> {
            Set<String> genes = new HashSet<>();
            instances.stream()
                     .filter(inst -> inst.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity))
                     .forEach(inst -> {
                         try {
                             ReactomeUtilities.grepGenesFromEntity(inst, genes);
                         }
                         catch(Exception e) {
                             logger.error(e.getMessage(), e);
                         }
                     });
            return genes;
        }).collect(Collectors.toSet());
        return setOfGenes;
    }
    
    /**
     * Expand a reaction into a list of sets of GKInstances, each of which should be a set of
     * instances that can carry out the passed reaction.
     * @param reaction
     * @throws Exception
     */
    private List<List<GKInstance>> expandReaction(GKInstance reaction) throws Exception {
        DirectedGraph<InstanceVertex, DefaultEdge> graph = convertToGraph(reaction);
        Set<InstanceVertex> sources = graph
                .vertexSet()
                .stream()
                .filter(inst -> graph.incomingEdgesOf(inst).size() == 0).collect(Collectors.toSet());
        Set<InstanceVertex> sinks = graph
                .vertexSet()
                .stream()
                .filter(inst -> graph.outgoingEdgesOf(inst).size() == 0).collect(Collectors.toSet());
        // Make sure there is only one targets
        if (sinks.size() > 1)
            throw new IllegalStateException("Converted graph from " + reaction + " has more than one sink vertex!");
        
        AllDirectedPaths<InstanceVertex, DefaultEdge> helper = new AllDirectedPaths<>(graph);
        List<GraphPath<InstanceVertex, DefaultEdge>> paths = helper.getAllPaths(sources, 
                                                                                sinks, 
                                                                                true, 
                                                                                null);
        List<List<GKInstance>> rtn = paths.stream().map(path -> 
                    path.getVertexList().stream().map(v -> v.instance).collect(Collectors.toList())
                ).collect(Collectors.toList());
        return rtn;
    }
    
    /**
     * Convert a reaction into a DirectedGraph.
     * @param reaction
     * @return
     * @throws Exception
     */
    private DirectedGraph<InstanceVertex, DefaultEdge> convertToGraph(GKInstance reaction) throws Exception {
        DirectedGraph<InstanceVertex, DefaultEdge> graph = new DefaultDirectedGraph<>(DefaultEdge.class);
        InstanceVertex rxtVertex = new InstanceVertex(reaction);
        graph.addVertex(rxtVertex);
        // To pick up genes involved in a reaction, we don't need to study outputs.
        Set<GKInstance> leftSideEntities = getLeftSideEntities(reaction);
        List<GKInstance> list = new ArrayList<>(leftSideEntities);
        List<InstanceVertex> vertexList = list.stream().map(inst -> new InstanceVertex(inst)).collect(Collectors.toList());
        vertexList.stream().forEach(graph::addVertex);
        
        for (int i = 0; i < list.size() - 1; i++) {
            InstanceVertex v1 = vertexList.get(i);
            InstanceVertex v2 = vertexList.get(i + 1);
            graph.addEdge(v1, 
                          v2);
        }
        graph.addEdge(vertexList.get(list.size() - 1), rxtVertex);
        // Create a linear graph so that all instances should be included in one travesal
        expandEntity(graph);
        // A sanity check to make sure the converted graph should not contain any cycle
        CycleDetector<InstanceVertex, DefaultEdge> cycleDetector = new CycleDetector<>(graph);
        if (cycleDetector.detectCycles())
            throw new IllegalStateException("Converted graph from " + reaction + " has cycles!");
        return graph;
    }
    
    /**
     * Expand the vertices representing EntitySets or Complexes until no more these types of
     * instance existing in the graph.
     * @param graph
     * @throws Exception
     */
    private void expandEntity(DirectedGraph<InstanceVertex, DefaultEdge> graph) throws Exception {
        boolean isUpdated = true; // First time
        while (isUpdated) { // Recursively calling until nothing to be changed 
            isUpdated = false;
            // Need to copy this set so that we can manipulate the graph.
            Set<InstanceVertex> copy = new HashSet<>(graph.vertexSet());
            for (InstanceVertex v : copy) {
                GKInstance inst = v.instance;
                // Escape the last reaction
                if (!inst.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity))
                    continue;
                if (inst.getSchemClass().isa(ReactomeJavaConstants.Complex)) {
                    isUpdated = expandComplex(graph, v);
                }
                else if (inst.getSchemClass().isa(ReactomeJavaConstants.EntitySet)) {
                    isUpdated = expandEntitySet(graph, v);
                }
            }
        }
    }

    private boolean expandEntitySet(DirectedGraph<InstanceVertex, DefaultEdge> graph, 
                                    InstanceVertex v) throws Exception {
        Set<GKInstance> members = new HashSet<>();
        members.addAll(v.instance.getAttributeValuesList(ReactomeJavaConstants.hasMember));
        if (v.instance.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasCandidate))
            members.addAll(v.instance.getAttributeValuesList(ReactomeJavaConstants.hasCandidate));
        if (members.size() == 0)
            return false;
        Set<InstanceVertex> vMembers = members.stream().map(m -> new InstanceVertex(m)).collect(Collectors.toSet());
        vMembers.stream().forEach(graph::addVertex);
        vMembers.stream().forEach(member -> {
            graph.incomingEdgesOf(v).stream().forEach(edge -> 
            graph.addEdge(graph.getEdgeSource(edge), member));
            graph.outgoingEdgesOf(v).stream().forEach(edge -> 
            graph.addEdge(member, graph.getEdgeTarget(edge)));
        });
        graph.removeVertex(v);
        return true;
    }

    private boolean expandComplex(DirectedGraph<InstanceVertex, DefaultEdge> graph, 
                                  InstanceVertex v) throws Exception {
        List<GKInstance> components = v.instance.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
        if (components == null || components.size() == 0)
            return false;
        List<InstanceVertex> vComps = components.stream().map(comp -> new InstanceVertex(comp)).collect(Collectors.toList());
        vComps.stream().forEach(graph::addVertex);
        for (int i = 0; i < components.size() - 1; i++) {
            InstanceVertex comp1 = vComps.get(i);
            InstanceVertex comp2 = vComps.get(i + 1);
            graph.addEdge(comp1, comp2);
        }
        graph.incomingEdgesOf(v).stream().forEach(edge -> 
        graph.addEdge(graph.getEdgeSource(edge), vComps.get(0)));
        graph.outgoingEdgesOf(v).stream().forEach(edge -> 
        graph.addEdge(vComps.get(vComps.size() - 1), graph.getEdgeTarget(edge)));
        // Now we can remove this vertex 
        graph.removeVertex(v);
        return true;
    }
    
    private Set<GKInstance> getLeftSideEntities(GKInstance reaction) throws Exception {
        Set<GKInstance> set = new HashSet<GKInstance>();
        List<GKInstance> inputs = reaction.getAttributeValuesList(ReactomeJavaConstants.input);
        if (inputs != null)
            set.addAll(inputs);
        List<GKInstance> cas = reaction.getAttributeValuesList(ReactomeJavaConstants.catalystActivity);
        if (cas != null && cas.size() > 0) {
            for (GKInstance ca : cas) {
                GKInstance catalyst = (GKInstance) ca.getAttributeValue(ReactomeJavaConstants.physicalEntity);
                if (catalyst != null)
                    set.add(catalyst);
            }
        }
        Collection<GKInstance> regulations = reaction.getReferers(ReactomeJavaConstants.regulatedEntity);
        if (regulations != null && regulations.size() > 0) {
            for (GKInstance regulation : regulations) {
                GKInstance regulator = (GKInstance) regulation.getAttributeValue(ReactomeJavaConstants.regulator);
                if (regulator == null)
                    continue;
                // Only take physical entity
                if (regulator.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity))
                    set.add(regulator);
            }
        }
        return set;
    }
      
    /**
     * Use a wrap around GKInstance to avoid circles or complicated graph structure if one or more instances
     * are reused in different places along the generated graph for the expanded reaction. By using this class,
     * multiple copies of InstanceVertex wrapping the same GKInstance may be generated, which guarantee a simple
     * graph (DAG) having a single output.
     * @author gwu
     *
     */
    private class InstanceVertex {
        private GKInstance instance;
        
        public InstanceVertex(GKInstance inst) {
            this.instance = inst;
        }
        
        public String toString() {
            return instance.toString();
        }
        
    }
    
}
