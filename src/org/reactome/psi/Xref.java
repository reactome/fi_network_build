/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

import java.util.ArrayList;
import java.util.List;

public class Xref {
    private Long dbId;
    private DbReference primaryRef;
    // secondaryRef should be a list
    private List<DbReference> secondaryRefList;
    
    public Xref() {
    }
    
    public void setDbId(Long id) {
        this.dbId = id;
    }
    
    public Long getDbId() {
        return this.dbId;
    }

    public DbReference getPrimaryRef() {
        return primaryRef;
    }

    public void setPrimaryRef(DbReference primaryRef) {
        this.primaryRef = primaryRef;
    }
    
    public void addSecondaryRef(DbReference secondaryRef) {
        if (secondaryRefList == null)
            secondaryRefList = new ArrayList<DbReference>();
        secondaryRefList.add(secondaryRef);
    }

    public List<DbReference> getSecondaryRefList() {
        return secondaryRefList;
    }

    public void setSecondaryRefList(List<DbReference> secondaryRefList) {
        this.secondaryRefList = secondaryRefList;
    }

}
