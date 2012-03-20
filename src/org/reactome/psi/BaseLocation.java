/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

/**
 * Have to pay attention for these int properties. Default 0 might mean something in the 
 * actual cases.
 * @author guanming
 *
 */
public class BaseLocation {
    private Long dbId;
    private OpenCV startStatus;
    // Only one of them can be specified
    private int begin;
    private Interval beginInterval;
    private OpenCV endStatus;
    private int end;
    private Interval endInterval;
    private boolean isLink;
    
    public BaseLocation() {
        
    }
    
    public void setDbId(Long id) {
        this.dbId = id;
    }
    
    public Long getDbId() {
        return this.dbId;
    }

    public int getBegin() {
        return begin;
    }

    public void setBegin(int begin) {
        this.begin = begin;
    }

    public Interval getBeginInterval() {
        return beginInterval;
    }

    public void setBeginInterval(Interval beginInterval) {
        this.beginInterval = beginInterval;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public Interval getEndInterval() {
        return endInterval;
    }

    public void setEndInterval(Interval endInterval) {
        this.endInterval = endInterval;
    }

    public OpenCV getEndStatus() {
        return endStatus;
    }

    public void setEndStatus(OpenCV endStatus) {
        this.endStatus = endStatus;
    }

    public boolean getIsLink() {
        return isLink;
    }

    public void setIsLink(boolean isLink) {
        this.isLink = isLink;
    }

    public OpenCV getStartStatus() {
        return startStatus;
    }

    public void setStartStatus(OpenCV startStatus) {
        this.startStatus = startStatus;
    }
    
    /**
     * Used to model base interval. This inner class is defined as static so that an Interval
     * instance can be created without an initialization of the outer instance.
     * @author guanming
     *
     */
    public static class Interval {
        private int begin;
        private int end;
        
        public Interval() {
        }
        
        public void setBegin(int start) {
            this.begin = start;
        }
        
        public int getBegin() {
            return this.begin;
        }
        
        public void setEnd(int end) {
            this.end = end;
        }
        
        public int getEnd() {
            return this.end;
        }
        
    }
    
}
