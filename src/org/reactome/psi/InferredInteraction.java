/*
 * Created on Mar 27, 2006
 *
 */
package org.reactome.psi;

import java.util.ArrayList;
import java.util.List;

public class InferredInteraction {
    private Long dbId;
    private List<Participant> participants;
    private List<Feature> participantFeatures;
    private List<Experiment> experimentRefList;
    
    public InferredInteraction() {
        
    }
    
    public Long getDbId() {
        return dbId;
    }

    public void setDbId(Long dbId) {
        this.dbId = dbId;
    }

    public List<Experiment> getExperimentRefList() {
        return experimentRefList;
    }

    public void setExperimentRefList(List<Experiment> experimentRefList) {
        this.experimentRefList = experimentRefList;
    }
    
    public void addExperimentRef(Experiment experiment) {
        if (experimentRefList == null)
            experimentRefList = new ArrayList<Experiment>();
        experimentRefList.add(experiment);
    }

    public List<Participant> getParticipants() {
        return participants;
    }

    public void addParticipant(Participant newValue) {
        if (participants == null)
            participants = new ArrayList<Participant>();
        participants.add(newValue);
    }
    
    public void setParticipants(List<Participant> participant) {
        this.participants = participant;
    }
    
    public List<Feature> getParticipantFeatures() {
        return this.participantFeatures;
    }
    
    public void setParticipantFeatures(List<Feature> list) {
        this.participantFeatures = list;
    }
    
    public void addParticipantFeature(Feature newFeature) {
        if (participantFeatures == null)
            participantFeatures = new ArrayList<Feature>();
        participantFeatures.add(newFeature);
    }
    
    
}
