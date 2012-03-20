package org.reactome.psi;

import java.util.ArrayList;
import java.util.List;

/**
 * This inner class is used to link an OpenCV object to a list of Experiments.
 * @author guanming
 */
public class OpenCVExperimentalWrapper {
    // It seems that dbId is required by hibernate OR mapping
    Long dbId;
    OpenCV openCV;
    List<Experiment> experiments;
    
    public OpenCVExperimentalWrapper() {
    }
    
    public void setNames(Names names) {
        if (openCV == null)
            openCV = new OpenCV();
        openCV.setNames(names);
    }
    
    public void setXref(Xref xref) {
        if (openCV == null)
            openCV = new OpenCV();
        openCV.setXref(xref);
    }
    
    public void setAttributeList(List<Attribute> attributeList) {
        if (openCV == null)
            openCV = new OpenCV();
        openCV.setAttributeList(attributeList);
    }
    
    public void addAttribute(Attribute att) {
        if (openCV == null)
            openCV = new OpenCV();
        openCV.addAttribute(att);
    }
    
    public void addExperiment(Experiment experiment) {
        if (experiments == null)
            experiments = new ArrayList<Experiment>();
        experiments.add(experiment);
    }
    
    public void setDbId(Long id) {
        this.dbId = id;
    }
    
    public Long getDbId() {
        return this.dbId;
    }
    
    public void setOpenCV(OpenCV cv) {
        this.openCV = cv;
    }
    
    public OpenCV getOpenCV() {
       return this.openCV;
    }
    
    
    
    public void setExperiments(List<Experiment> experiments) {
        this.experiments = experiments;
    }
    
    public List<Experiment> getExperiments() {
        return this.experiments;
    }
}