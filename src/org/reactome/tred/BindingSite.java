/*
 * Created on Apr 9, 2009
 *
 */
package org.reactome.tred;

public class BindingSite {
    private Integer id;
    private FactorPromoter factorPromoter;
    private Integer genomeId;
    private String chromosome;
    private String strand;
    private Integer start;
    private Integer end;
    private String siteSequence;
    
    public BindingSite() {
    }

    public Integer getId() {
        return id;
    }

    public void setId(Integer id) {
        this.id = id;
    }

    public FactorPromoter getFactorPromoter() {
        return factorPromoter;
    }

    public void setFactorPromoter(FactorPromoter factorPromoter) {
        this.factorPromoter = factorPromoter;
    }

    public Integer getGenomeId() {
        return genomeId;
    }

    public void setGenomeId(Integer genomeId) {
        this.genomeId = genomeId;
    }

    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
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

    public String getSiteSequence() {
        return siteSequence;
    }

    public void setSiteSequence(String siteSequence) {
        this.siteSequence = siteSequence;
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }
    
}
