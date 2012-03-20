/*
 * Created on Apr 9, 2009
 *
 */
package org.reactome.tred;

public class FactorPromoterEvidence {
    private Source source;
    private String sourceAccession;
    private String sourceDescription;
    private String sourceExtract;
    private Integer sourceRank;
    private String alert;
    
    public FactorPromoterEvidence() {
    }

    public Source getSource() {
        return source;
    }

    public void setSource(Source source) {
        this.source = source;
    }

    public String getSourceAccession() {
        return sourceAccession;
    }

    public void setSourceAccession(String sourceAccession) {
        this.sourceAccession = sourceAccession;
    }

    public String getSourceDescription() {
        return sourceDescription;
    }

    public void setSourceDescription(String sourceDescription) {
        this.sourceDescription = sourceDescription;
    }

    public String getSourceExtract() {
        return sourceExtract;
    }

    public void setSourceExtract(String sourceExtract) {
        this.sourceExtract = sourceExtract;
    }

    public Integer getSourceRank() {
        return sourceRank;
    }

    public void setSourceRank(Integer sourceRank) {
        this.sourceRank = sourceRank;
    }

    public String getAlert() {
        return alert;
    }

    public void setAlert(String alert) {
        this.alert = alert;
    }

}
