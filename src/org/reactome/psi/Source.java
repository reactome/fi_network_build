/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

import java.sql.Date;
import java.util.ArrayList;
import java.util.List;

public class Source {
    // This is for persistence purpose, and is not required by the PSI XML schema
    private Long dbId;
    private String release;
    private Date releaseDate;
    private Names names;
    private Bibref bibref;
    private Xref xref;
    private List<Attribute> attributeList;
    
    public Source() {        
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
    
    public void setRelease(String release) {
        this.release = release;
    }
    
    public String getRelease() {
        return release;
    }
    
    public void setReleaseDate(Date date) {
        this.releaseDate = date;
    }
    
    public Date getReleaseDate() {
        return releaseDate;
    }

    public Bibref getBibref() {
        return bibref;
    }

    public void setBibref(Bibref bibref) {
        this.bibref = bibref;
    }

    public Names getNames() {
        return names;
    }

    public void setNames(Names names) {
        this.names = names;
    }

    public Xref getXref() {
        return xref;
    }

    public void setXref(Xref xref) {
        this.xref = xref;
    }

    public List<Attribute> getAttributeList() {
        return attributeList;
    }

    public void setAttributeList(List<Attribute> attributeList) {
        this.attributeList = attributeList;
    }
}
