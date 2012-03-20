/*
 * Created on Sep 28, 2006
 *
 */
package org.reactome.funcInt;

import javax.xml.bind.annotation.XmlRootElement;

/**
 * This class is used to describe evidence for predicted functional interactions.
 * @author guanming
 *
 */
@XmlRootElement
public class Evidence {
    
    private long dbId;
    private Boolean humanInteraction;
    @Deprecated
    private Boolean orthoInteraction;
    @Deprecated
    private Boolean yeastInteraction;
    @Deprecated
    private GeneExpressionType geneExp;
    @Deprecated
    private Integer goBPSemanticSimilarity;
    // Following are new attributes as of April 15, 2009
    private Boolean dmePPI;
    private Boolean celPPI;
    private Boolean scePPI;
    private Boolean mousePPI;
    private Boolean genewaysPPI;
    private Boolean pfamDomainInt;
    private Boolean goBPSharing;
    private Boolean pavlidisGeneExp;
    private Boolean carlosGeneExp;
    // Probability for true functional interaction
    private double probability;
    
    public Boolean getOrthoInteraction() {
        return orthoInteraction;
    }

    public void setOrthoInteraction(Boolean orthoInteraction) {
        this.orthoInteraction = orthoInteraction;
    }

    public Boolean getYeastInteraction() {
        return yeastInteraction;
    }

    public void setYeastInteraction(Boolean yeastInteraction) {
        this.yeastInteraction = yeastInteraction;
    }

    public GeneExpressionType getGeneExp() {
        return geneExp;
    }

    public void setGeneExp(GeneExpressionType geneExp) {
        this.geneExp = geneExp;
    }

    public Integer getGoBPSemanticSimilarity() {
        return goBPSemanticSimilarity;
    }

    public void setGoBPSemanticSimilarity(Integer goBPSemanticSimilarity) {
        this.goBPSemanticSimilarity = goBPSemanticSimilarity;
    }

    public Evidence() {
    }

    
    public Boolean getMousePPI() {
        return mousePPI;
    }

    public void setMousePPI(Boolean mousePPI) {
        this.mousePPI = mousePPI;
    }
    
    public long getDbId() {
        return dbId;
    }

    public void setDbId(long dbId) {
        this.dbId = dbId;
    }

    public Boolean getHumanInteraction() {
        return humanInteraction;
    }

    public void setHumanInteraction(Boolean humanInteraction) {
        this.humanInteraction = humanInteraction;
    }

    public Boolean getDmePPI() {
        return dmePPI;
    }

    public void setDmePPI(Boolean dmePPI) {
        this.dmePPI = dmePPI;
    }

    public Boolean getCelPPI() {
        return celPPI;
    }

    public void setCelPPI(Boolean celPPI) {
        this.celPPI = celPPI;
    }

    public Boolean getScePPI() {
        return scePPI;
    }

    public void setScePPI(Boolean scePPI) {
        this.scePPI = scePPI;
    }

    public Boolean getGenewaysPPI() {
        return genewaysPPI;
    }

    public void setGenewaysPPI(Boolean genewaysPPI) {
        this.genewaysPPI = genewaysPPI;
    }

    public Boolean getPfamDomainInt() {
        return pfamDomainInt;
    }

    public void setPfamDomainInt(Boolean pfamDomainInt) {
        this.pfamDomainInt = pfamDomainInt;
    }

    public Boolean getGoBPSharing() {
        return goBPSharing;
    }

    public void setGoBPSharing(Boolean goBPSharing) {
        this.goBPSharing = goBPSharing;
    }

    public Boolean getPavlidisGeneExp() {
        return pavlidisGeneExp;
    }

    public void setPavlidisGeneExp(Boolean pavlidisGeneExp) {
        this.pavlidisGeneExp = pavlidisGeneExp;
    }

    public Boolean getCarlosGeneExp() {
        return carlosGeneExp;
    }

    public void setCarlosGeneExp(Boolean carlosGeneExp) {
        this.carlosGeneExp = carlosGeneExp;
    }

    public double getProbability() {
        return probability;
    }

    public void setProbability(double probability) {
        this.probability = probability;
    }

}
