/*
 * Created on Apr 27, 2006
 *
 */
package org.reactome.psi.data;

/**
 * A simplified version of Interactor compared to Interactor in PSI-MI package. Only UniProt
 * ID, name, type are used.
 * @author guanming
 *
 */
public class NTInteractor {
    private Long dbId;
    private String id;
    private String name;
    private NTInteractorType type;
    
    public NTInteractor() {
    }

    public Long getDbId() {
        return dbId;
    }

    public void setDbId(Long dbId) {
        this.dbId = dbId;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public NTInteractorType getType() {
        return type;
    }

    public void setType(NTInteractorType type) {
        this.type = type;
    }
    
    public void exportToXML(StringBuilder xmlString) {
        xmlString.append("<interactor dbId=\"").append(dbId).append("\"");
        xmlString.append(" id=\"").append(id).append("\"");
        xmlString.append(" name=\"").append(name).append("\"");
        xmlString.append(" type=\"").append(type).append("\" />");
    }
}
