package org.reactome.psi;

import java.util.ArrayList;
import java.util.List;

/**
 * All inner classes in this class are defined as static so that they can be instantiated
 * without the encosing class, Participant.
 * @author guanming
 *
 */
public class ExperimentalConfidence {
    Long dbId;
    Confidence confidence;
    List<Experiment> experiments;
    
    public ExperimentalConfidence() {
        
    }
    
    public void setUnit(OpenCV unit) {
        if (confidence == null)
            confidence = new Confidence();
        confidence.setUnit(unit);
    }
    
    public OpenCV getUnit() {
        return (confidence == null) ? null : confidence.getUnit();
    }
    
    public void setValue(String value) {
        if (confidence == null)
            confidence = new Confidence();
        confidence.setValue(value);
    }
    
    public String getValue() {
        return (confidence == null) ? null : confidence.getValue();
    }
    public void addExperiment(Experiment exp) {
        if (experiments == null)
            experiments = new ArrayList<Experiment>();
        experiments.add(exp);
    }
    
    public void setDbId(Long id) {
        this.dbId = id;
    }
    
    public Long getDbId() {
        return this.dbId;
    }

    public Confidence getConfidence() {
        return confidence;
    }

    public void setConfidence(Confidence confidence) {
        this.confidence = confidence;
    }

    public List<Experiment> getExperiments() {
        return experiments;
    }

    public void setExperiments(List<Experiment> experiments) {
        this.experiments = experiments;
    }
}