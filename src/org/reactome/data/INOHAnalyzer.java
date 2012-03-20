/*
 * Created on Aug 8, 2006
 *
 */
package org.reactome.data;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.PersistenceAdaptor;
import org.gk.model.ReactomeJavaConstants;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;
import org.reactome.fi.util.FileUtility;


public class INOHAnalyzer extends PantherAnalyzer {
    private final String DIR_NAME = "/Users/wgm/Documents/caBIG_R3/datasets/INOH/";
    
    public INOHAnalyzer() {
        dataSourceId = 239582L; // v2
    }
    
    public void extractLocationOntology() throws IOException {
        String fileName = DIR_NAME + "LocationOntology_100/LocationOntology_100.obo";
        Map<String, String> inoh2go = new HashMap<String, String>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        String id = null;
        String goId = null;
        while ((line = fu.readLine()) != null) {
            if (line.equals("[Term]")) {
                // Start a new term
                id = goId = null;
            }
            else if (line.startsWith("id: ")) {
                // ID Line
                id = line.substring(4);
            }
            else if (line.startsWith("xref_analog: GO:")) { // GO dbxref
                goId = line.substring(13);
                inoh2go.put(id, goId);
            }
        }
        fu.close();
        fu.exportMap(inoh2go, DIR_NAME + "INOHLocation2GO.txt");
    }
    
    protected Set<GKInstance> grepRefPepSeqs(GKInstance interactor) throws Exception {
        Set<GKInstance> refPepSeqs = super.grepRefPepSeqs(interactor);
        // Have to exclude ReferencePeptideSequences from other species. In INOH, all
        // RefPepSeqs are grouped together from different species
        GKInstance human = dba.fetchInstance(ReactomeJavaConstants.Species, 48887L);
        // Take only human species
        for (Iterator<GKInstance> it = refPepSeqs.iterator(); it.hasNext();) {
            GKInstance refPepSeq = it.next();
            GKInstance species = (GKInstance) refPepSeq.getAttributeValue(ReactomeJavaConstants.species);
            if (species != human)
                it.remove();
        }
        return refPepSeqs;
    }
    
    protected List<GKInstance> getTopics() throws Exception {
        PersistenceAdaptor dba = getMySQLAdaptor();
        GKInstance dataSource = getDataSource();
        Collection collection = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Pathway, 
                                                             ReactomeJavaConstants.dataSource,
                                                             "=",
                                                             dataSource);
        List<GKInstance> topics = new ArrayList<GKInstance>();
        for (Iterator it = collection.iterator(); it.hasNext();) {
            GKInstance gkInstance = (GKInstance) it.next();
            Collection containers = gkInstance.getReferers(ReactomeJavaConstants.hasComponent);
            if (containers == null || containers.size() == 0)
                topics.add(gkInstance);
        }
        return topics;
    }
    
    public void extractUniProtFromMolecularRoleOntology() throws IOException {
        String fileName = DIR_NAME + "MoleculeRoleOntology_211/MoleculeRoleOntology_211.obo";
        // File for reading
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        // File for writing
        FileUtility outputFu = new FileUtility();
        String outFileName = DIR_NAME + "INOHRole2UniProt.txt";
        outputFu.setOutput(outFileName);
        String line = null;
        String uniId = null;
        String parentId = null;
        String name = null;
        while ((line = fu.readLine()) != null) {
            if (line.equals("[Term]")) {
                // Want to get human uniprot only
                if (name != null &&
                    name.endsWith("_HUMAN") &&
                    uniId != null && 
                    parentId != null)
                    outputFu.printLine(parentId + " " + uniId);
                name = uniId = parentId = null;
            }
            else if (line.startsWith("name: ")) {
                name = line.substring(6);
            }
            else if (line.startsWith("xref_analog: UniProt:")) {
                uniId = line.substring(21);
            }
            else if (line.startsWith("relationship: sequence_of ")) {
                parentId = line.substring(26, 37);
            }
        }
        fu.close();
        outputFu.close();
    }
    
    public void testINOHFile() throws Exception {
        String fileName = DIR_NAME + "AllDiagram_060217/EGF.inoh";
        SAXBuilder builder = new SAXBuilder();
        Document document = builder.build(new File(fileName));
        Set<String> nodeTypes = new HashSet<String>();
        Element root = document.getRootElement();
        List nodes = root.getChildren("node");
        for (Iterator it = nodes.iterator(); it.hasNext();) {
            Element node = (Element) it.next();
            nodeTypes.add(node.getAttributeValue("Type"));
        }
        System.out.println("Node Types: " + nodeTypes.size() + "\n" + nodeTypes);
        Set<String> edgeTypes = new HashSet<String>();
        List edges = root.getChildren("edge");
        for (Iterator it = edges.iterator(); it.hasNext();) {
            Element edge = (Element) it.next();
            edgeTypes.add(edge.getAttributeValue("Type"));
        }
        System.out.println("Edge types: " + edgeTypes.size() + "\n" + edgeTypes);
    }
    
}
