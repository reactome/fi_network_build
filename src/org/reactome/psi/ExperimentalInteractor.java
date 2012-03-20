/*
 * Created on Mar 24, 2006
 *
 */
package org.reactome.psi;

import java.util.ArrayList;
import java.util.List;

public class ExperimentalInteractor {
    Long dbId;
    Interactor interactor;
    List<Experiment> experiments;
    
    public ExperimentalInteractor() {
        
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

    public List<Experiment> getExperiments() {
        return experiments;
    }

    public void setExperiments(List<Experiment> experiments) {
        this.experiments = experiments;
    }

    public Interactor getInteractor() {
        return interactor;
    }

    public void setInteractor(Interactor interactor) {
        this.interactor = interactor;
    }
    
}