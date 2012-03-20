/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

import java.util.ArrayList;
import java.util.List;

public class Feature {
    private Long dbId;
    private int id;
    private Names names;
    private Xref xref;
    private OpenCV featureType;
    private OpenCV featureDetectionMethod;
    private List<Experiment> experimentList;
    private List<BaseLocation> featureRangeList;
    private List<Attribute> attributeList;
    
    public Feature() {
    }
    
    public void addExperiment(Experiment experiment) {
        if (experimentList == null)
            experimentList = new ArrayList<Experiment>();
        experimentList.add(experiment);
    }
    
    public void addFeatureRange(BaseLocation location) {
        if (featureRangeList == null)
            featureRangeList = new ArrayList<BaseLocation>();
        featureRangeList.add(location);
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

    public List<Attribute> getAttributeList() {
        return attributeList;
    }

    public void setAttributeList(List<Attribute> attributeList) {
        this.attributeList = attributeList;
    }

    public List<Experiment> getExperimentList() {
        return experimentList;
    }

    public void setExperimentList(List<Experiment> experimentList) {
        this.experimentList = experimentList;
    }

    public OpenCV getFeatureDetectionMethod() {
        return featureDetectionMethod;
    }

    public void setFeatureDetectionMethod(OpenCV featureDetectionMethod) {
        this.featureDetectionMethod = featureDetectionMethod;
    }

    public List<BaseLocation> getFeatureRangeList() {
        return featureRangeList;
    }

    public void setFeatureRangeList(List<BaseLocation> featureRangeList) {
        this.featureRangeList = featureRangeList;
    }

    public OpenCV getFeatureType() {
        return featureType;
    }

    public void setFeatureType(OpenCV featureType) {
        this.featureType = featureType;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
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
    
}
