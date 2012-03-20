/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

import java.util.ArrayList;
import java.util.List;

public class Experiment {
    // For database identifier, internally.
    private Long dbId;
    private String id;
    private Names names;
    private Bibref bibref;
    private Xref xref;
    private List<BioSource> hostOrganismList;
    // controlled by the PSI-MI controlled vocabulary 
    // "interaction detection method", root term id MI:0001.
    private OpenCV interactionDetectionMethod;
    // controlled by the PSI-MI controlled vocabulary 
    // "participant identification method", root term id MI:0002.
    private OpenCV participantIdentificationMethod;
    // controlled by the PSI-MI controlled vocabulary 
    // "feature detection method", root term id MI:0003.
    private OpenCV featureDetectionMethod;
    private List<ExperimentalConfidence> confidenceList;
    private List<Attribute> attributeList;
    
    public Experiment() {
    }
    
    public void addHostOrganism(BioSource biosource) {
        if (hostOrganismList == null)
            hostOrganismList = new ArrayList<BioSource>();
        hostOrganismList.add(biosource);
    }
    
    public void addConfidence(ExperimentalConfidence confidence) {
        if (confidenceList == null)
            confidenceList = new ArrayList<ExperimentalConfidence>();
        confidenceList.add(confidence);
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

    public void setAttributeList(List<Attribute> attribute) {
        this.attributeList = attribute;
    }

    public Bibref getBibref() {
        return bibref;
    }

    public void setBibref(Bibref bibref) {
        this.bibref = bibref;
    }

    public List<ExperimentalConfidence> getConfidenceList() {
        return confidenceList;
    }

    public void setConfidenceList(List<ExperimentalConfidence> confidenceList) {
        this.confidenceList = confidenceList;
    }

    public OpenCV getFeatureDetectionMethod() {
        return featureDetectionMethod;
    }

    public void setFeatureDetectionMethod(OpenCV featureDetectionMethod) {
        this.featureDetectionMethod = featureDetectionMethod;
    }

    public List<BioSource> getHostOrganismList() {
        return hostOrganismList;
    }

    public void setHostOrganismList(List<BioSource> hostOrganismList) {
        this.hostOrganismList = hostOrganismList;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public OpenCV getInteractionDetectionMethod() {
        return interactionDetectionMethod;
    }

    public void setInteractionDetectionMethod(OpenCV interactiondetectionMethod) {
        this.interactionDetectionMethod = interactiondetectionMethod;
    }

    public Names getNames() {
        return names;
    }

    public void setNames(Names names) {
        this.names = names;
    }

    public OpenCV getParticipantIdentificationMethod() {
        return participantIdentificationMethod;
    }

    public void setParticipantIdentificationMethod(
            OpenCV participantIdentificationMethod) {
        this.participantIdentificationMethod = participantIdentificationMethod;
    }

    public Xref getXref() {
        return xref;
    }

    public void setXref(Xref xref) {
        this.xref = xref;
    }

}
