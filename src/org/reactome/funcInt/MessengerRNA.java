/*
 * Created on Jan 20, 2009
 *
 */
package org.reactome.funcInt;

/**
 * This class is used to model mRNA, which will be used in miRNA and target 
 * interaction.
 * @author wgm
 *
 */
public class MessengerRNA extends GenomeEncodedEntity {
    
    private Protein protein;
    
    public MessengerRNA() {
    }

    public Protein getProtein() {
        return protein;
    }

    public void setProtein(Protein protein) {
        this.protein = protein;
    }
    
    public boolean equals(Object mRNA) {
        if (!(mRNA instanceof MessengerRNA))
            return false;
        MessengerRNA out = (MessengerRNA) mRNA;
        return out.getCheckSum().equals(getCheckSum());
    }
    
}
