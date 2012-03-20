/*
 * Created on Jan 20, 2009
 *
 */
package org.reactome.funcInt;

/**
 * This model class is used to model miRNA, which will be used in miRNA and target
 * interaction.
 * @author wgm
 *
 */
public class MicroRNA {
    
    private long dbId;
    private String accession;
    private String sequence;
    private String checksum;
    
    public MicroRNA() {
    }

    public long getDbId() {
        return dbId;
    }

    public void setDbId(long dbId) {
        this.dbId = dbId;
    }

    public String getAccession() {
        return accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public String getChecksum() {
        return checksum;
    }

    public void setChecksum(String checksum) {
        this.checksum = checksum;
    }
    
}
