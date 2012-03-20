/*
 * Created on Apr 10, 2009
 *
 */
package org.reactome.tred;

public class RefGene {
    private String name;
    private String chrom;
    private String strand;
    private Integer txStart;
    private Integer txEnd;
    private Integer cdsStart;
    private Integer cdsEnd;
    private Integer exonCount;
    private String exonStarts;
    private String exonEnds;
    private String speciesCode;
    
    public RefGene() {
        
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getChrom() {
        return chrom;
    }

    public void setChrom(String chrom) {
        this.chrom = chrom;
    }

    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public Integer getTxStart() {
        return txStart;
    }

    public void setTxStart(Integer txStart) {
        this.txStart = txStart;
    }

    public Integer getTxEnd() {
        return txEnd;
    }

    public void setTxEnd(Integer txEnd) {
        this.txEnd = txEnd;
    }

    public Integer getCdsStart() {
        return cdsStart;
    }

    public void setCdsStart(Integer cdsStart) {
        this.cdsStart = cdsStart;
    }

    public Integer getCdsEnd() {
        return cdsEnd;
    }

    public void setCdsEnd(Integer cdsEnd) {
        this.cdsEnd = cdsEnd;
    }

    public Integer getExonCount() {
        return exonCount;
    }

    public void setExonCount(Integer exonCount) {
        this.exonCount = exonCount;
    }

    public String getExonStarts() {
        return exonStarts;
    }

    public void setExonStarts(String exonStarts) {
        this.exonStarts = exonStarts;
    }

    public String getExonEnds() {
        return exonEnds;
    }

    public void setExonEnds(String exonEnds) {
        this.exonEnds = exonEnds;
    }

    public String getSpeciesCode() {
        return speciesCode;
    }

    public void setSpeciesCode(String speciseCode) {
        this.speciesCode = speciseCode;
    }
    
}
