package org.reactome.panther;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.gk.model.GKInstance;

/**
 * An internal Protein class to keep the information in celldesigner:protein elements.
 */
class Protein {
    String type;
    String proteinId;
    String name;
    private Map<String, String> positionMap;
    // One protein might point to multiple PhysicalEntities.
    List<GKInstance> entities;
    
    public Protein(String proteinId, String name, String type) {
        this.proteinId = proteinId;
        this.name = name;
        this.type = type;
    }
    
    public String getProteinId() {
        return this.proteinId;
    }
    
    public String getName() {
        return name;
    }
    
    public String getType() {
        return type;
    }
    
    public void addPositions(String id, String position) {
        if (positionMap == null)
            positionMap = new HashMap<String, String>();
        positionMap.put(id, position);
    }
    
    public String getPosition(String id) {
        if (positionMap == null)
            return null;
        return positionMap.get(id);
    }
    
    public void addEntity(GKInstance entity) {
        if (entities == null)
            entities = new ArrayList<GKInstance>();
        entities.add(entity);
    }
    
    public List<GKInstance> getEntities() {
        return entities;
    }
    
    public String toString() {
    	return this.getName();
    }
}