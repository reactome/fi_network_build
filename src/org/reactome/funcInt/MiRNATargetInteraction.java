/*
 * Created on Jan 20, 2009
 *
 */
package org.reactome.funcInt;

/**
 * This class is used to model miRNA and mRNA interaction. If no score has been attached to an 
 * MiRNATargetInteraction, the interaction should be annotated instead of predicted.
 * 
 * @author wgm
 *
 */
public class MiRNATargetInteraction {
    
    private long dbId;
    private MiRNATargetScore score;
    // Some prediction algorithm using two scores. For example,
    // miRanda.
    private MiRNATargetScore secondaryScore;
    private MicroRNA miRNA;
    private MessengerRNA mRNA;
    
    public MiRNATargetInteraction() {
    }

    public long getDbId() {
        return dbId;
    }

    public void setDbId(long dbId) {
        this.dbId = dbId;
    }

    public MiRNATargetScore getScore() {
        return score;
    }

    public void setScore(MiRNATargetScore score) {
        this.score = score;
    }

    public MiRNATargetScore getSecondaryScore() {
        return secondaryScore;
    }

    public void setSecondaryScore(MiRNATargetScore secondaryScore) {
        this.secondaryScore = secondaryScore;
    }

    public MicroRNA getMiRNA() {
        return miRNA;
    }

    public void setMiRNA(MicroRNA miRNA) {
        this.miRNA = miRNA;
    }

    /**
     * Have to use lower-case m in mRNA to make hibernate work.
     * @return
     */
    public MessengerRNA getmRNA() {
        return mRNA;
    }
    
    public void setmRNA(MessengerRNA mrna) {
        mRNA = mrna;
    }
    
}
