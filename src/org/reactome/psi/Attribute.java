/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

public class Attribute {
    private Long dbId;
    private String name;
    private String nameAc;
    private String text;
    
    public Attribute() {
    }
    
    public void setText(String text) {
        this.text = text;
    }
    
    public String getText() {
        return this.text;
    }
    
    public void setDbId(Long id) {
        this.dbId = id;
    }
    
    public Long getDbId() {
        return this.dbId;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getNameAc() {
        return nameAc;
    }

    public void setNameAc(String nameAc) {
        this.nameAc = nameAc;
    }
    
}
