/*
 * Created on Jul 1, 2010
 *
 */
package org.reactome.fi;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.reactome.fi.util.FileUtility;

/**
 * This class is used to load GO terms.
 * @author wgm
 *
 */
public class GOTermLoader {
    private String goaFileName;
    private String goTermIdFileName;
    private String proteinIdToNameFileName;
    private FileUtility fu;
    
    public GOTermLoader() {
        fu = new FileUtility();
    }
    
    public String getGoaFileName() {
        return goaFileName;
    }

    public String getGoTermIdFileName() {
        return goTermIdFileName;
    }

    public void setGoTermIdFileName(String goTermIdFileName) {
        this.goTermIdFileName = goTermIdFileName;
    }

    public String getProteinIdToNameFileName() {
        return proteinIdToNameFileName;
    }

    public void setProteinIdToNameFileName(String proteinIdToNameFileName) {
        this.proteinIdToNameFileName = proteinIdToNameFileName;
    }

    public void setGoaFileName(String goaFileName) {
        this.goaFileName = goaFileName;
    }

    public Map<String, Set<String>> loadProteinToGOBPTerms() throws IOException {
        String aspect = "P";
        return loadProteinToGOTerms(aspect);
    }
    
    public Map<String, Set<String>> loadProteinToGOMFTerms() throws IOException {
        return loadProteinToGOTerms("F");
    }
    
    public Map<String, Set<String>> loadProteinToGOCCTerms() throws IOException {
        return loadProteinToGOTerms("C");
    }

    private Map<String, Set<String>> loadProteinToGOTerms(String aspect) throws IOException {
        // The following terms should be escaped
        String[] escapedTerms = new String[] {
                "GO:0008150", // biological_process
                "GO:0003674", // molecular_function
                "GO:0005575" // cellular_component
        };
        // For each check
        List<String> escapeList = Arrays.asList(escapedTerms);
        FileUtility fu = new FileUtility();
        fu.setInput(goaFileName);
        String line = null;
        Map<String, Set<String>> proteinToTerms = new HashMap<String, Set<String>>();
        while ((line = fu.readLine()) != null) {
//            System.out.println(line);
            if (line.startsWith("!") || line.startsWith("\"!"))
                continue; // Comments
            String[] tokens = line.split("\t");
            // Check aspect 
            if (!tokens[8].equals(aspect))
                continue;
            String protein = tokens[1];
            String term = tokens[4];
            if (escapeList.contains(term))
                continue;
            Set<String> termSet = proteinToTerms.get(protein);
            if (termSet == null) {
                termSet = new HashSet<String>();
                proteinToTerms.put(protein, termSet);
            }
            termSet.add(term);
        }
        fu.close();
        return proteinToTerms;
    }
    
    public Map<String, Set<String>> convertProteinIdToNameForGO(Map<String, Set<String>> proteinToGOTerms) throws IOException {
        Map<String, String> idToName = fu.importMap(proteinIdToNameFileName);
        Map<String, Set<String>> geneToGoTerms = new HashMap<String, Set<String>>();
        for (String protein : proteinToGOTerms.keySet()) {
            Set<String> goTerms = proteinToGOTerms.get(protein);
            String gene = idToName.get(protein);
            if (gene == null)
                continue;
            geneToGoTerms.put(gene, goTerms);
        }
        return geneToGoTerms;
    }

    public Map<String, Set<String>> convertGOIdsToTerms(Map<String, Set<String>> geneToIds) throws IOException {
        Map<String, String> idToName = loadGOIdToTermMap();
        Map<String, Set<String>> geneToNames = new HashMap<String, Set<String>>();
        for (String gene : geneToIds.keySet()) {
            Set<String> ids = geneToIds.get(gene);
            Set<String> names = new HashSet<String>(ids.size());
            for (String id : ids) {
                String name = idToName.get(id);
                names.add(name);
            }
            geneToNames.put(gene, names);
        }
        return geneToNames;
    }
    
    public Map<String, String> loadGOIdToTermMap() throws IOException {
        fu.setInput(goTermIdFileName);
        String line = null;
        Map<String, String> id2Term = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("!"))
                continue;
            String[] tokens = line.split("\t");
            id2Term.put(tokens[0], tokens[1]);
        }
        fu.close();
        return id2Term;
    }
}
