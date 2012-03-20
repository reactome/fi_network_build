/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

/**
 * This class is used for both ConfidenceType and ConfidenceListType. ConfidenceListType has a 
 * new element called ExperimentRef, which seems rather strange for use. Just ignore it so that
 * Confidence can be attached to Experiment objects directly.
 * @author guanming
 *
 */
public class Confidence {
    private Long dbId;
    private OpenCV unit;
    private String value;
    
    public Confidence() {
    }
    
    public void setDbId(Long id) {
        this.dbId = id;
    }
    
    public Long getDbId() {
        return this.dbId;
    }

    public OpenCV getUnit() {
        return unit;
    }

    public void setUnit(OpenCV unit) {
        this.unit = unit;
    }

    public String getValue() {
        return value;
    }

    public void setValue(String value) {
        this.value = value;
    }
    
}
