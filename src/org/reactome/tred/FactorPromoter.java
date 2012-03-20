/*
 * Created on Apr 9, 2009
 *
 */
package org.reactome.tred;

import java.util.Set;

public class FactorPromoter {
    
    private Integer id;
    private Factor factor;
    private Promoter promoter;
    private String siteName;
    private String siteSequence;
    private Integer startOffset;
    private Integer endOffset;
    private FactorPromoterQuality quality;
    private Set<FactorPromoterEvidence> evidences;
    
    public FactorPromoter() {
    }

    public Integer getId() {
        return id;
    }

    public void setId(Integer id) {
        this.id = id;
    }

    public Factor getFactor() {
        return factor;
    }

    public void setFactor(Factor factor) {
        this.factor = factor;
    }

    public Promoter getPromoter() {
        return promoter;
    }

    public void setPromoter(Promoter promoter) {
        this.promoter = promoter;
    }

    public String getSiteName() {
        return siteName;
    }

    public void setSiteName(String siteName) {
        this.siteName = siteName;
    }

    public String getSiteSequence() {
        return siteSequence;
    }

    public void setSiteSequence(String siteSequence) {
        this.siteSequence = siteSequence;
    }

    public Integer getStartOffset() {
        return startOffset;
    }

    public void setStartOffset(Integer startOffset) {
        this.startOffset = startOffset;
    }

    public Integer getEndOffset() {
        return endOffset;
    }

    public void setEndOffset(Integer endOffset) {
        this.endOffset = endOffset;
    }

    public FactorPromoterQuality getQuality() {
        return quality;
    }

    public void setQuality(FactorPromoterQuality quality) {
        this.quality = quality;
    }

    public Set<FactorPromoterEvidence> getEvidences() {
        return evidences;
    }

    public void setEvidences(Set<FactorPromoterEvidence> evidences) {
        this.evidences = evidences;
    }
    
}
