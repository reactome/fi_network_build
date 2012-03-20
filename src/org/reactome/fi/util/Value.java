package org.reactome.fi.util;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * This inner class is used to hold values
 * @author guanming
 */
public class Value {
    public Boolean functionalInteraction;
    public Boolean humanInteraction;
    // Human PPIs from three different human PPI databases
    public Boolean intactHumanPPI;
    public Boolean hprdHumanPPI;
    public Boolean biogridHumanPPI;
    // Mapped PPIs from three species
    public Boolean dmePPI;
    public Boolean celPPI;
    public Boolean scePPI;
    public Boolean mousePPI;
    // Gene Co-expression: two data sets
    public Boolean pavlidisGeneExp;
    public Boolean carlosGeneExp;
    
    public Boolean orthoInteraction;
    public Boolean yeastInteraction;
    public Boolean y2h;
    public Boolean nonY2h;
    public Boolean flyNonY2H;
    public Boolean wormY2H;
    public Boolean flyY2H; 
    public Integer goBPOccurence;
    public Integer goBPDepth;
    public Integer goMFOccurence;
    public Integer goMFDepth;
    //Boolean localization;
    public String geneExp;
    public List<Float> geneExpList;
    // For GO BP term sharing
    public Boolean goBPSharing;
    public Boolean goMFSharing;
    public Boolean goCCSharing;
    public Double goMFSemSimilarity;
    public Integer goBPSemSimilarity;
    // For Pfam domain interactions
    public Boolean pfamDomainInt;
    // For geneways interactions
    public Boolean genewaysPPI;
    public Boolean unambGenewaysPPI;
    public Boolean ambGenewaysPPI;
    public Boolean unambGenewaysFlybaseDmePPI;
    // For tissue protein expression
    public Double tissueExp;
    /*
     * The following features stands for
     * Theileria parva:
	286131
	Caenorhabditis elegans:
	216051
	Dictyostelium discoideum:
	310914
	Gallus gallus:
	342698
	Halobacterium sp. NRC-1:
	192607
	Danio rerio:
	188996
	Toxoplasma gondii:
	305940
	Drosophila melanogaster:
	361401
	Encephalitozoon cuniculi:
	216410
	Schizosaccharomyces pombe:
	363515
	Caenorhabditis briggsae:
	349703
	Tetraodon nigroviridis:
	236620
	Homo sapiens:
	272924
	Arabidopsis thaliana:
	298426
	Plasmodium falciparum 3D7:
	323000
	Entamoeba histolytica:
	312340
	Oryza sativa:
	372618
	Mus musculus:
	364033
	Plasmodium knowlesi:
	283609
	Rattus norvegicus:
	347274
	Fugu rubripes:
	348581
	Ashbya gossypii:
	220206
	Ciona intestinalis:
	349741
	Yeast
	239492*/
   
//    public Boolean f342698;
//    public Boolean f364033;
//    public Boolean f347274;
//    public Boolean f239492;

    
    
    
    public Value() {}
    
    public void addGeneExpValue(Float f) {
        if (geneExpList == null)
            geneExpList = new ArrayList<Float>();
        geneExpList.add(f);
    }
    
    /**
     * This method is used to convert a list of feature in String to Field used
     * in the Value class.
     * @return
     */
    public static Map<String, Field> convertValueFeatureToField(List<String> featureList) {
        Class<Value> valueCls = Value.class;
        Field[] fields = valueCls.getDeclaredFields();
        Map<String, Field> featureToField = new HashMap<String, Field>();
        boolean isFound = false;
        for (String feature : featureList) {
            isFound = false;
            for (Field f : fields) {
                if (f.getName().equals(feature)) {
                    featureToField.put(feature, f);
                    isFound = true;
                    break;
                }
            }
            if (!isFound) {
                throw new IllegalStateException("NaviveBayesClassifier.convertFeatureToField(): " + 
                                                feature + " cannot be mapped to field in Class value.");
            }
        }
        return featureToField;
    }
    
    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append("functionInteraction: ").append(functionalInteraction);
        builder.append(",humanInteraction: ").append(humanInteraction);
        builder.append(",orthoInteraction: ").append(orthoInteraction);
        builder.append(",yeastInterction: ").append(yeastInteraction);
        builder.append(",geneExp: ").append(geneExp);
        builder.append(",goBPSemSim: ").append(goBPSemSimilarity);
        builder.append(",unambGenewaysPPI: ").append(unambGenewaysPPI);
        builder.append(",ambGenewaysPPI: ").append(ambGenewaysPPI);
        builder.append(",tissueExp: ").append(tissueExp);
        return builder.toString();
    }
}