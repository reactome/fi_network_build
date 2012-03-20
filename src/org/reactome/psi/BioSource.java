/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

public class BioSource {
    private Long dbId;
    private int ncbiTaxId;
    private Names names;
    private OpenCV cellType;
    private OpenCV compartment;
    private OpenCV tissue;
    
    public BioSource() {
        
    }
    
    public void setDbId(Long id) {
        this.dbId = id;
    }
    
    public Long getDbId() {
        return this.dbId;
    }

    public OpenCV getCellType() {
        return cellType;
    }

    public void setCellType(OpenCV cellType) {
        this.cellType = cellType;
    }

    public OpenCV getCompartment() {
        return compartment;
    }

    public void setCompartment(OpenCV compartment) {
        this.compartment = compartment;
    }

    public Names getNames() {
        return names;
    }

    public void setNames(Names names) {
        this.names = names;
    }

    public int getNcbiTaxId() {
        return ncbiTaxId;
    }

    public void setNcbiTaxId(int ncbiTaxId) {
        this.ncbiTaxId = ncbiTaxId;
    }

    public OpenCV getTissue() {
        return tissue;
    }

    public void setTissue(OpenCV tissue) {
        this.tissue = tissue;
    }
    

}
