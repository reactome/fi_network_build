/*
 * Created on Apr 9, 2009
 *
 */
package org.reactome.tred;

import java.util.Set;

/**
 * Note: There are two tables regarding Promoter: Promoter and FirstExonEndNullifiedPromoter. The second
 * table is more like the modification of the first table for Promoters whose fistExonEnd are null. In the
 * second table, this information is provided. However, this Java data model will not handle this information.
 * @author wgm
 *
 */
public class Promoter {
    private Integer id;
    private String name;
    private Gene gene;
    private String chromosome;
    private String strand;
    private String tss;
    private String sequence;
    private Integer upstreamOffset;
    private Integer downstreamOffset;
    private Integer firstExonEnd;
    private Integer firstIntronEnd;
    private PromoterQuality quality;
    private String pidZX;
    private Integer curated;
    private String cpg;
    private Set<PromoterGeneEvidence> geneEvidences;
    
    
    public Set<PromoterGeneEvidence> getGeneEvidences() {
        return geneEvidences;
    }


    public void setGeneEvidences(Set<PromoterGeneEvidence> geneEvidences) {
        this.geneEvidences = geneEvidences;
    }


    public Promoter() {
    }


    public Integer getId() {
        return id;
    }


    public void setId(Integer id) {
        this.id = id;
    }


    public String getName() {
        return name;
    }


    public void setName(String name) {
        this.name = name;
    }


    public Gene getGene() {
        return gene;
    }


    public void setGene(Gene gene) {
        this.gene = gene;
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


    public String getTss() {
        return tss;
    }


    public void setTss(String tss) {
        this.tss = tss;
    }


    public String getSequence() {
        return sequence;
    }


    public void setSequence(String sequence) {
        this.sequence = sequence;
    }


    public Integer getUpstreamOffset() {
        return upstreamOffset;
    }


    public void setUpstreamOffset(Integer upstreamOffset) {
        this.upstreamOffset = upstreamOffset;
    }


    public Integer getDownstreamOffset() {
        return downstreamOffset;
    }


    public void setDownstreamOffset(Integer downstreamOffset) {
        this.downstreamOffset = downstreamOffset;
    }


    public Integer getFirstExonEnd() {
        return firstExonEnd;
    }


    public void setFirstExonEnd(Integer firstExonEnd) {
        this.firstExonEnd = firstExonEnd;
    }


    public Integer getFirstIntronEnd() {
        return firstIntronEnd;
    }


    public void setFirstIntronEnd(Integer firstIntronEnd) {
        this.firstIntronEnd = firstIntronEnd;
    }


    public PromoterQuality getQuality() {
        return quality;
    }


    public void setQuality(PromoterQuality quality) {
        this.quality = quality;
    }


    public String getPidZX() {
        return pidZX;
    }


    public void setPidZX(String pidZX) {
        this.pidZX = pidZX;
    }


    public Integer getCurated() {
        return curated;
    }


    public void setCurated(Integer curated) {
        this.curated = curated;
    }


    public String getCpg() {
        return cpg;
    }


    public void setCpg(String cpg) {
        this.cpg = cpg;
    }
    
    
    
}
