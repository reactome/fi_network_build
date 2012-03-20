/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

import java.util.ArrayList;
import java.util.List;

public class Bibref {
    private Long dbId;
    private Xref xref;
    private List<Attribute> attributeList;
    
    public Bibref() {
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
    
    public void addAttribute(Attribute attribute) {
        if (attributeList == null)
            attributeList = new ArrayList<Attribute>();
        attributeList.add(attribute);
    }

    public void setAttributeList(List<Attribute> attribute) {
        this.attributeList = attribute;
    }

    public Xref getXref() {
        return xref;
    }

    public void setXref(Xref xref) {
        this.xref = xref;
    }
}
