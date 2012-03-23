/*
 * Created on Mar 21, 2012
 *
 */
package org.reactome.fi;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.jgrapht.Graph;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;

/**
 * This class is used to do some graph related file generation, e.g. generating
 * a largest graph component.
 * @author gwu
 *
 */
public class FIGraphAnalyzer {
    
    public FIGraphAnalyzer() {
    }
    
    /**
     * Analyze graph components and output the biggest graph component.
     * @throws Exception
     */
    @Test
    public void analyzeComponents() throws Exception {
        FileUtility fu = new FileUtility();
        Set<String> interactions = fu.loadInteractions(FIConfiguration.getConfiguration().get("GENE_FI_FILE_NAME"));
        System.out.println("Total interactions: " + interactions.size());
        List<Set<String>> componentList = calculateGraphComponents(interactions);
        System.out.println("Total components: " + componentList.size());
        int index = 0;
        for (Set<String> comp : componentList) {
            System.out.println(index + ": " + comp.size());
            index ++;
        }
        Set<String> biggestComp = componentList.get(0);
        // Want to print out the biggest components
        Set<String> interactionsInBiggest = new HashSet<String>();
        index = 0;
        for (String in : interactions) {
            index = in.indexOf("\t");
            if (index < 0)
                index = in.indexOf(" ");
            String id1 = in.substring(0, index);
            String id2 = in.substring(index + 1);
            if (biggestComp.contains(id1) &&
                biggestComp.contains(id2))
                interactionsInBiggest.add(in);
        }
        fu.saveInteractions(interactionsInBiggest, 
                            FIConfiguration.getConfiguration().get("GENE_FI_BIG_COMP_FILE_NAME"));
    }
    
    /**
     * Calculate linked graph components from a set of FIs.
     * @param interactions
     * @return
     */
    public List<Set<String>> calculateGraphComponents(Set<String> interactions) {
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
        SimpleGraph<String, DefaultEdge> graph = (SimpleGraph<String, DefaultEdge>) createGraph(ids, interactions);
        ConnectivityInspector<String, DefaultEdge> inspector = new ConnectivityInspector<String, DefaultEdge>(graph);
        List<Set<String>> components = inspector.connectedSets();
        List<Set<String>> componentList = new ArrayList<Set<String>>(components);
        Collections.sort(componentList, new Comparator<Set<String>>() {
            public int compare(Set<String> set1, Set<String> set2) {
                return set2.size() - set1.size();
            }
        });
        return componentList;
    }
    
    public Graph<String, DefaultEdge> createGraph(Set<String> reactomeIds, 
                                                  Set<String> interactions) {
        int index = 0;
        // Just want to use String as edge class
        Graph<String, DefaultEdge> graph = new SimpleGraph<String, DefaultEdge>(DefaultEdge.class);
        for (String id : reactomeIds)
            graph.addVertex(id);
        // Need to find the delimit index
        String tmp = interactions.iterator().next();
        String delimit = " ";
        index = tmp.indexOf("\t");
        if (index > 0)
            delimit = "\t"; 
        for (String pair : interactions) {
            index = pair.indexOf(delimit);
            //System.out.println(pair);
            graph.addEdge(pair.substring(0, index), pair.substring(index + 1));
        }
//        System.out.printf("Graph: vertics %d edges %d%n", 
//                          graph.vertexSet().size(), 
//                          graph.edgeSet().size());        
        return graph;
    }
}
