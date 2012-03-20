/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

import java.util.ArrayList;
import java.util.List;

public class Interactor {
    private Long dbId;
    private int id;
    private Names names;
    private Xref xref;
    private OpenCV interactorType;
    private BioSource organism;
    private String sequence;
    // No Sequence is needed
    private List<Attribute> attributeList;
    
    public Interactor() {
        
    }
    
    public void setSequence(String sequence) {
        this.sequence = sequence;
    }
    
    public String getSequence() {
        return this.sequence;
    }
    
    public void setId(int id) {
        this.id = id;
    }
    
    public int getId() {
        return this.id;
    }
    
    public void addAttribute(Attribute att) {
        if (attributeList == null)
            attributeList = new ArrayList<Attribute>();
        attributeList.add(att);
    }
    
    public void setDbId(Long id) {
        this.dbId = id;
    }
    
    public Long getDbId() {
        return this.dbId;
    }

    public List<Attribute> getAttributeList() {
        return attributeList;
    }

    public void setAttributeList(List<Attribute> attributeList) {
        this.attributeList = attributeList;
    }

    public OpenCV getInteractorType() {
        return interactorType;
    }

    public void setInteractorType(OpenCV interactorType) {
        this.interactorType = interactorType;
    }

    public Names getNames() {
        return names;
    }

    public void setNames(Names names) {
        this.names = names;
    }

    public BioSource getOrganism() {
        return organism;
    }

    public void setOrganism(BioSource organism) {
        this.organism = organism;
    }

    public Xref getXref() {
        return xref;
    }

    public void setXref(Xref xref) {
        this.xref = xref;
    }

}
