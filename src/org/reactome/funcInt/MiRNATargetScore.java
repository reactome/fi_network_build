/*
 * Created on Jan 20, 2009
 *
 */
package org.reactome.funcInt;

/**
 * This class is used to model score from a miRNA target prediction algorithm.
 * @author wgm
 *
 */
public class MiRNATargetScore {
    
    private long dbId;
    // The actual value generated from the algorithm
    private double score;
    // The name of the algorithm (e.g. miRanda)
    private String algorithmName;
    // The name of the score (e.g. align score)
    private String scoreName;
    
    public MiRNATargetScore() {
    }

    public long getDbId() {
        return dbId;
    }

    public void setDbId(long dbId) {
        this.dbId = dbId;
    }

    public double getScore() {
        return score;
    }

    public void setScore(double score) {
        this.score = score;
    }

    public String getAlgorithmName() {
        return algorithmName;
    }

    public void setAlgorithmName(String algorithmName) {
        this.algorithmName = algorithmName;
    }

    public String getScoreName() {
        return scoreName;
    }

    public void setScoreName(String scoreName) {
        this.scoreName = scoreName;
    }
}
