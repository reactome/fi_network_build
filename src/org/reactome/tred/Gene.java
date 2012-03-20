/*
 * Created on Apr 9, 2009
 *
 */
package org.reactome.tred;

import java.util.Set;

public class Gene {
    private Integer id;
    private String primaryName;
    private String allNames;
    private String species;
    private Integer genomeId;
    private String chromosome;
    private String strand;
    private Integer start;
    private Integer end;
    private String chromLocation;
    private String gidZx;
    private Set<GeneAnnotation> annotations;
    
    public Gene() {
    }

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

    public String getAllNames() {
        return allNames;
    }

    public void setAllNames(String allNames) {
        this.allNames = allNames;
    }

    public String getSpecies() {
        return species;
    }

    public void setSpecies(String species) {
        this.species = species;
    }

    public Integer getGenomeId() {
        return genomeId;
    }

    public void setGenomeId(Integer genomeId) {
        this.genomeId = genomeId;
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public Set<GeneAnnotation> getAnnotations() {
        return annotations;
    }

    public void setAnnotations(Set<GeneAnnotation> annotations) {
        this.annotations = annotations;
    }

    public Integer getStart() {
        return start;
    }

    public void setStart(Integer start) {
        this.start = start;
    }

    public Integer getEnd() {
        return end;
    }

    public void setEnd(Integer end) {
        this.end = end;
    }

    public String getChromLocation() {
        return chromLocation;
    }

    public void setChromLocation(String chromLocation) {
        this.chromLocation = chromLocation;
    }

    public String getGidZx() {
        return gidZx;
    }

    public void setGidZx(String gidZx) {
        this.gidZx = gidZx;
    }
    
    
}
