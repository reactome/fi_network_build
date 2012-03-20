/*
 * Created on Mar 27, 2006
 *
 */
package org.reactome.psi;

public class Availability {
    private Long dbId;
    private int id;
    private String text;
    
    public Availability() {
        
    }
    
    public void setText(String text) {
        this.text = text;
    }
    
    public String getText() {
        return this.text;
    }

    public Long getDbId() {
        return dbId;
    }

    public void setDbId(Long dbId) {
        this.dbId = dbId;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

}
