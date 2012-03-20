/*
 * Created on Dec 1, 2006
 *
 */
package org.reactome.panther;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.reactome.convert.common.PostProcessHelper;
import org.reactome.fi.util.FileUtility;

/**
 * This helper class is used to handle mapping from small molecules in Panther
 * to ChEBI.
 * @author guanming
 *
 */
public class SmallMoleculeHandler {
    private String mappingFile = "resources/SpeciesToChEBIId.txt";
    private Map<String, List<String>> nameToIds;
    // cache the found instance
    private Map<String, GKInstance> cached = new HashMap<String, GKInstance>();
    
    public SmallMoleculeHandler() throws IOException {
        loadMapping();
    }
    
    private void loadMapping() throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(mappingFile);
        nameToIds = new HashMap<String, List<String>>();
        // 10-formyl-THF    15637   11304   19108   698
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            List<String> list = new ArrayList<String>();
            nameToIds.put(tokens[0], list);
            for (int i = 1; i < tokens.length; i++) {
                list.add(tokens[i]);
            }
        }
        fu.close();
    }
    
    @SuppressWarnings("unchecked")
    public void processSimpleEntities(MySQLAdaptor dbAdaptor,
                                      XMLFileAdaptor fileAdaptor) throws Exception {
        Collection seList = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.SimpleEntity);
        if (seList == null || seList.size() == 0)
            return;
        GKInstance se = null;
        String name = null;
        List<String> chebiIds;
        GKInstance refMol = null;
        for (Iterator it = seList.iterator(); it.hasNext();) {
            se = (GKInstance) it.next();
            // display name might contain compartment information. use
            // name instead
            name = (String) se.getAttributeValue(ReactomeJavaConstants.name);
            refMol = getReferenceMolecule(name, dbAdaptor);
            if (refMol != null) {
                refMol = PostProcessHelper.downloadDBInstance(refMol, fileAdaptor);
                se.addAttributeValue(ReactomeJavaConstants.referenceEntity, refMol);
            }
        }
    }
    
    private GKInstance getReferenceMolecule(String name,
                                            MySQLAdaptor dbAdaptor) throws Exception {
        if (cached.containsKey(name)) // Have to check this first. Since the value might be null
            return cached.get(name);
        List<String> chebiIds = nameToIds.get(name);
        if (chebiIds == null || chebiIds.size() == 0) {
            cached.put(name, null);
            return null;
        }
        GKInstance refMol = fetchReferenceMolecule(dbAdaptor, chebiIds);
        cached.put(name, refMol);
        return refMol;
    }
    
    private GKInstance fetchReferenceMolecule(MySQLAdaptor dbAdaptor,
                                              List<String> chebiIds) throws Exception {
        GKInstance rtn = null;
        for (String id : chebiIds) {
            Collection c = dbAdaptor.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceMolecule,
                                                              ReactomeJavaConstants.identifier, 
                                                              "=", 
                                                              id);
            if (c != null && c.size() > 0) {
                // Need one only
                rtn = (GKInstance) c.iterator().next();
                return rtn;
            }
        }
        return rtn;
    }
}
