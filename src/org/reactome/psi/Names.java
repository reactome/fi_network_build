/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

import java.util.ArrayList;
import java.util.List;

public class Names {
    private Long dbId;
    private String shortLabel;
    private String fullName;
    private List<String> alias;
    
    public Names() {
    }
    
    public void setDbId(Long id) {
        this.dbId = id;
    }
    
    public Long getDbId() {
        return this.dbId;
    }
    
    public void addAlias(String name) {
        if (alias == null) {
            alias = new ArrayList<String>();
        }
        alias.add(name);
    }

    public List<String> getAlias() {
        return alias;
    }

    public void setAlias(List<String> alias) {
        this.alias = alias;
    }

    public String getFullName() {
        return fullName;
    }

    public void setFullName(String fullName) {
        this.fullName = fullName;
    }

    public String getShortLabel() {
        return shortLabel;
    }

    public void setShortLabel(String shortLabel) {
        this.shortLabel = shortLabel;
    }
    
    public boolean equals(Object obj) {
        if (obj == this)
            return true;
        if (!(obj instanceof Names))
            return false;
        Names other = (Names) obj;
        if (!(other.getShortLabel().equals(shortLabel)))
            return false;
        if (!(other.getFullName().equals(fullName)))
            return false;
        return true;
    }
    
    public int hashCode() {
        if (fullName != null && shortLabel != null)
            return fullName.hashCode() * 29 + shortLabel.hashCode();
        if (fullName != null)
            return fullName.hashCode() * 29;
        if (shortLabel != null)
            return shortLabel.hashCode();
        return super.hashCode();
    }
    
}
