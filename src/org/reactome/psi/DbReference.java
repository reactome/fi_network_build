/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

import java.util.ArrayList;
import java.util.List;

public class DbReference {
    private Long dbId;
    private String db;
    private String dbAc;
    private String id;
    private String secondary;
    private String version;
    private String refType;
    private String refTypeAc;
    private List<Attribute> attributeList;
    
    public DbReference() {
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
    
    public void addAttribute(Attribute att) {
        if (attributeList == null)
            attributeList = new ArrayList<Attribute>();
        attributeList.add(att);
    }

    public void setAttributeList(List<Attribute> attributeList) {
        this.attributeList = attributeList;
    }

    public String getDb() {
        return db;
    }

    public void setDb(String db) {
        this.db = db;
    }

    public String getDbAc() {
        return dbAc;
    }

    public void setDbAc(String dbAc) {
        this.dbAc = dbAc;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getRefType() {
        return refType;
    }

    public void setRefType(String refType) {
        this.refType = refType;
    }

    public String getRefTypeAc() {
        return refTypeAc;
    }

    public void setRefTypeAc(String refTypeAc) {
        this.refTypeAc = refTypeAc;
    }

    public String getSecondary() {
        return secondary;
    }

    public void setSecondary(String secondary) {
        this.secondary = secondary;
    }

    public String getVersion() {
        return version;
    }

    public void setVersion(String version) {
        this.version = version;
    }
    
    public boolean equals(Object obj) {
        if (!(obj instanceof DbReference))
            return false;
        DbReference other = (DbReference) obj;
        String key = db + "|DbReference|" + id + "|DbReference|" + secondary + 
        "|DbReference|" + version + "|DbReference|" + refType;
        String otherKey = other.db + "|DbReference|" + other.id + "|DbReference|" + other.secondary +
        "|DbReference|" + other.version + "|DbReference|" + other.refType;
        return key.equals(otherKey);
    }
    
    public int hashCode() {
        String key = db + "|DbReference|" + id + "|DbReference|" + secondary + 
        "|DbReference|" + version + "|DbReference|" + refType;
        return key.hashCode();
    }

}
