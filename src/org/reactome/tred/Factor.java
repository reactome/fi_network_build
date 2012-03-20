/*
 * Created on Apr 9, 2009
 *
 */
package org.reactome.tred;

import java.util.Set;

public class Factor {
    private Integer id;
    private String primaryName;
    private String species;
    private String allNames;
    private String geneId;
    private Set<FactorAnnotation> annotations;
    
    public Integer getId() {
        return id;
    }

    public void setId(Integer id) {
        this.id = id;
    }

    public String getPrimaryName() {
        return primaryName;
    }

    public void setPrimaryName(String primaryName) {
        this.primaryName = primaryName;
    }

    public String getSpecies() {
        return species;
    }

    public void setSpecies(String species) {
        this.species = species;
    }

    public String getAllNames() {
        return allNames;
    }

    public void setAllNames(String allNames) {
        this.allNames = allNames;
    }

    public String getGeneId() {
        return geneId;
    }

    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public Set<FactorAnnotation> getAnnotations() {
        return annotations;
    }

    public void setAnnotations(Set<FactorAnnotation> annotations) {
        this.annotations = annotations;
    }
    
    
    
}
