/*
 * Created on Dec 2, 2009
 *
 */
package org.reactome.hibernate;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.gk.util.FileUtilities;
import org.hibernate.Query;
import org.hibernate.Session;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.funcInt.DbReference;
import org.reactome.funcInt.Evidence;
import org.reactome.funcInt.Interaction;
import org.reactome.funcInt.Protein;
import org.reactome.funcInt.ReactomeSource;

/**
 * This class is used to export FIs in the PSI MITAB format.
 * @author wgm
 *
 */
public class FIMITabExporter extends HibernateFIReader {
    private final String NA_LABEL = "-";
    
    public FIMITabExporter() {
        
    }
    
    public void export(String fileName) throws Exception {
        FileUtilities fu = new FileUtilities();
        fu.setOutput(fileName);
        exportHeaders(fu);
        // Start query
        initSession();
        Session session = sessionFactory.openSession();
        // Query predicted FIs
        Query query = session.createQuery("FROM Interaction as i WHERE i.evidence.probability >= ?");
        Double cutoff = new Double(FIConfiguration.getConfiguration().get("CUT_OFF_VALUE"));
        query.setDouble(0, cutoff);
        List list = query.list();
        System.out.println("Total interactions from prediction: " + list.size());
        Interaction interaction = null;
        for (Iterator it = list.iterator(); it.hasNext();) {
            interaction = (Interaction) it.next();
            exportFI(interaction, fu);
        }
        // Query extracted FIs
        query = session.createQuery("FROM Interaction as i WHERE i.evidence is null");
        list = query.list();
        System.out.println("Total interactions from pathways: " + list.size());
        for (Iterator it = list.iterator(); it.hasNext();) {
            interaction = (Interaction) it.next();
            exportFI(interaction, fu);
        }
        fu.close();
    }
    
    private String getDbName(DbReference dbRef) {
        String dbName = dbRef.getDbName();
        if (dbName.equals("UniProt"))
            dbName = "uniprotkb";
        return dbName;
    }
    
    private void exportFI(Interaction interaction, FileUtilities fu) throws Exception {
        StringBuilder builder = new StringBuilder();
        // Output primary db information
        Protein firstProtein = interaction.getFirstProtein();
        DbReference firstDb = firstProtein.getPrimaryDbReference();
        String dbName = getDbName(firstDb);
        builder.append(dbName).append(":").append(firstDb.getAccession());
        Protein secondProtein = interaction.getSecondProtein();
        DbReference secondDb = secondProtein.getPrimaryDbReference();
        dbName = getDbName(secondDb);
        builder.append("\t").append(dbName).append(":").append(secondDb.getAccession());
        // Output gene names
        String name = getProteinName(firstProtein);
        builder.append("\t").append(name);
        name = getProteinName(secondProtein);
        builder.append("\t").append(name);
        // No aliases are handled
        builder.append("\t").append(NA_LABEL).append("\t").append(NA_LABEL);
        // Interaction detection method: a rather generic term is used since we cannot find an appropriate term under this
        // category.
        // Have to wrap the second ":" in a quotation mark. Otherwise, it cannot be recognized
        // by the web service app.
        builder.append("\t").append("psi-mi:\"MI:0046\"(experimental knowledge based)");
        // Use myself as the first author
        builder.append("\t").append("Wu et al.(2010)");
        // Not published yet
        // Published: use PMID 20482850
        builder.append("\tpubmed:20482850");
        // taxon: always homo sapiens for both proteins
        builder.append("\ttaxid:9606\ttaxid:9606");
        // Interaction type: no approprate term can be found
        builder.append("\t").append(NA_LABEL);
        // Source DB: use Reactome for all
        builder.append("\tpsi-mi:\"MI:0467\"(reactome)");
        // interaction in source db
        builder.append("\t").append(NA_LABEL);
        // Confidence score
        if (interaction.getEvidence() == null)
            builder.append("\t").append(NA_LABEL);
        else
            builder.append("\tNBC:").append(interaction.getEvidence().getProbability());
        // Optional: pathway source for extracted FIs
        if (interaction.getEvidence() == null) {
            Set<ReactomeSource> sources = interaction.getReactomeSources();
            builder.append("\t");
            for (Iterator<ReactomeSource> it = sources.iterator(); it.hasNext();) {
                ReactomeSource src = it.next();
                builder.append(src.getDataSource()).append(":");
                builder.append(src.getSourceType()).append(":");
                builder.append(src.getReactomeId());
                if (it.hasNext())
                    builder.append("|");
            }
        }
        else {
            builder.append("\t").append(NA_LABEL);
        }
        // Optional: evidences for predicted
        if (interaction.getEvidence() != null) {
            // Want to list positive evidences only
            List<String> list = new ArrayList<String>();
            Evidence evidence = interaction.getEvidence();
            if (evidence.getHumanInteraction())
                list.add("humanInteraction");
            if (evidence.getDmePPI())
                list.add("dmePPI");
            if (evidence.getCelPPI())
                list.add("celPPI");
            if (evidence.getScePPI())
                list.add("scePPI");
            if (evidence.getGenewaysPPI())
                list.add("genewaysPPI");
            if (evidence.getPfamDomainInt())
                list.add("pfamDomaintInt");
            if (evidence.getGoBPSharing())
                list.add("goBPSharing");
            if (evidence.getPavlidisGeneExp())
                list.add("pavlidisGeneExp");
            if (evidence.getCarlosGeneExp())
                list.add("carlosGeneExp");
            builder.append("\t");
            for (Iterator<String> it = list.iterator(); it.hasNext();) {
                String evi = it.next();
                builder.append(evi);
                if (it.hasNext())
                    builder.append("|");
            }
        }
        else {
            builder.append("\t").append(NA_LABEL);
        }
        fu.printLine(builder.toString());
    }
    
    private String getProteinName(Protein protein) {
        String name = protein.getShortName();
        if (name == null)
            name = NA_LABEL;
        else {
            String dbName = getDbName(protein.getPrimaryDbReference());
            name = dbName + ":" + name;
        }
        return name;
    }
    
    /**
     * This helper method is used to generate headers for the MITAB format.
     * @param fu
     * @throws IOException
     */
    private void exportHeaders(FileUtilities fu) throws IOException {
        String header = "unique id A\t" +
        		        "unique id B\t" + 
                        "alternative id A\t" +
                        "alternative id B\t" + 
                        "aliases for A\t" +
                        "aliases for B\t" + 
                        "interaction detection method\t" + 
                        "first author\t" +
                        "Id of publication\t" + 
                        "NCBI taxon A\t" +
                        "NCBI taxon B\t" + 
                        "interaction type\t" +
                        "source db\t" + 
                        "interaction in source db\t" + 
                        "confidence score\t" + // Score from NBC will be used.
                        "pathway source\t" + // Optional: the source for extracted FIs
                        "evidences"; // Optional: evidences for predicted FIs.
        fu.printLine(header);
    }
    
    @Test
    public void generateMITABFormat() throws Exception {
        String fileName = "tmp/FIsInMITTab.txt";
        export(fileName);
    }
    
}
