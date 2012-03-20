/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

import java.util.ArrayList;
import java.util.List;

public class Participant {
    private Long dbId;
    private int id;
    private Names names;
    private Xref xref;
    // This is for the choice element: interactorRef, interactor, interactionRef.
    // InteractionRef will NOT be supported here.
    private Interactor interactor;
    // Original name in the XML schema is participantIdentificationMethodList.
    private List<OpenCVExperimentalWrapper> identificationMethodList;
    private OpenCV biologicalRole;
    private List<OpenCVExperimentalWrapper> experimentalRoleList;
    private List<OpenCVExperimentalWrapper> experimentalPreparationList;
    private List<ExperimentalInteractor> experimentalInteractorList;
    private List<Feature> featureList;
    private List<BioSource> hostOrganismList;
    private List<ExperimentalConfidence> confidenceList;
    private List<Attribute> attributeList;
    
    public Participant() {
    }
    
    public void addConfidence(ExperimentalConfidence confidence) {
        if (confidenceList == null)
            confidenceList = new ArrayList<ExperimentalConfidence>();
        confidenceList.add(confidence);
    }
    
    public void addHostOrganism(BioSource source) {
        if (hostOrganismList == null)
            hostOrganismList = new ArrayList<BioSource>();
        hostOrganismList.add(source);
    }
    
    public void addParticipantIdentificationMethod(OpenCVExperimentalWrapper wrapper) {
        if (identificationMethodList == null)
            identificationMethodList = new ArrayList<OpenCVExperimentalWrapper>();
        identificationMethodList.add(wrapper);
    }
    
    public void addExperimentalRole(OpenCVExperimentalWrapper wrapper) {
        if (experimentalRoleList == null)
            experimentalRoleList = new ArrayList<OpenCVExperimentalWrapper>();
        experimentalRoleList.add(wrapper);
    }
    
    public void addExperimentalPreparation(OpenCVExperimentalWrapper wrapper) {
        if (experimentalPreparationList == null)
            experimentalPreparationList = new ArrayList<OpenCVExperimentalWrapper>();
        experimentalPreparationList.add(wrapper);
    }
    
    public void addExperimentalInteractor(ExperimentalInteractor interactor) {
        if (experimentalInteractorList == null)
            experimentalInteractorList = new ArrayList<ExperimentalInteractor>();
        experimentalInteractorList.add(interactor);
    }
    
    public void addFeature(Feature feature) {
        if (featureList == null)
            featureList = new ArrayList<Feature>();
        featureList.add(feature);
    }
    
    public void setDbId(Long id) {
        this.dbId = id;
    }
    
    public Long getDbId() {
        return this.dbId;
    }
   
    public OpenCV getBiologicalRole() {
        return biologicalRole;
    }

    public void setBiologicalRole(OpenCV biologicalRole) {
        this.biologicalRole = biologicalRole;
    }

    public List<ExperimentalConfidence> getConfidenceList() {
        return confidenceList;
    }

    public void setConfidenceList(List<ExperimentalConfidence> confidenceList) {
        this.confidenceList = confidenceList;
    }

    public List<ExperimentalInteractor> getExperimentalInteractorList() {
        return experimentalInteractorList;
    }

    public void setExperimentalInteractorList(
            List<ExperimentalInteractor> experimentalInteractorList) {
        this.experimentalInteractorList = experimentalInteractorList;
    }

    public List<OpenCVExperimentalWrapper> getExperimentalPreparationList() {
        return experimentalPreparationList;
    }

    public void setExperimentalPreparationList(
            List<OpenCVExperimentalWrapper> experimentalPreparationList) {
        this.experimentalPreparationList = experimentalPreparationList;
    }

    public List<OpenCVExperimentalWrapper> getExperimentalRoleList() {
        return experimentalRoleList;
    }

    public void setExperimentalRoleList(List<OpenCVExperimentalWrapper> experimentalRoleList) {
        this.experimentalRoleList = experimentalRoleList;
    }

    public List<Feature> getFeatureList() {
        return featureList;
    }

    public void setFeatureList(List<Feature> featureList) {
        this.featureList = featureList;
    }

    public List<BioSource> getHostOrganismList() {
        return hostOrganismList;
    }

    public void setHostOrganismList(List<BioSource> hostOrganismList) {
        this.hostOrganismList = hostOrganismList;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public Interactor getInteractor() {
        return interactor;
    }

    public void setInteractor(Interactor interactor) {
        this.interactor = interactor;
    }

    public Names getNames() {
        return names;
    }

    public void setNames(Names names) {
        this.names = names;
    }

    public List<OpenCVExperimentalWrapper> getParticipantIdentificationMethodList() {
        return identificationMethodList;
    }

    public void setParticipantIdentificationMethodList(
            List<OpenCVExperimentalWrapper> identificationMethodList) {
        this.identificationMethodList = identificationMethodList;
    }

    public Xref getXref() {
        return xref;
    }

    public void setXref(Xref xref) {
        this.xref = xref;
    }
    
    public List<Attribute> getAttributeList() {
        return attributeList;
    }

    public void setAttributeList(List<Attribute> attributeList) {
        this.attributeList = attributeList;
    }
    
    public void addAttribute(Attribute attribute) {
        if (attributeList == null)
            attributeList = new ArrayList<Attribute>();
        attributeList.add(attribute);
    }
    
}
