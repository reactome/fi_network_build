/*
 * Created on Jul 1, 2010
 *
 */
package org.reactome.data;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
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
    
    /**
     * As of 2016, only go.obo is downloaded and then processed to generate GO.terms_and_ids.txt
     * file to keep to update. The original file was generated in 2012 and not udpated.
     * @throws IOException
     */
    @Test
    public void generateGOTermsToIdsFromObo() throws IOException {
        String goDir = FIConfiguration.getConfiguration().get("GO_DIR");
        String srcFile = goDir + "go.obo";
        String targetFile = goDir + "GO.terms_and_ids.txt";
        fu.setInput(srcFile);
        fu.setOutput(targetFile);
        // Output follows the following format as in the original GO.terms_and_ids.txt
        // GO:0000000 [tab] text string [tab] F|P|C
        // where F = molecular function, P = biological process, C = cellular component
        String line = null;
        boolean isInTerm = false;
        StringBuilder builder = new StringBuilder();
        while ((line = fu.readLine()) != null) {
            if (line.trim().length() == 0) {
                isInTerm = false;
            }
            else if (line.equals("[Term]"))
                isInTerm = true;
            else if (isInTerm) {
                if (line.startsWith("id:")) {
                    builder.append(extractValue(line));
                }
                else if (line.startsWith("name:"))
                    builder.append("\t").append(extractValue(line));
                else if (line.startsWith("namespace:")) {
                    builder.append("\t").append(getNameSapce(line));
                    fu.printLine(builder.toString());
                    builder.setLength(0);
                }
            }
        }
        fu.close();
    }
    
    private String extractValue(String line) {
        int index = line.indexOf(":");
        return line.substring(index + 1).trim();
    }
    
    private String getNameSapce(String line) {
        String value = extractValue(line);
        if (value.equals("molecular_function"))
            return "F";
        if (value.equals("biological_process"))
            return "P";
        if (value.equals("cellular_component"))
            return "C";
        throw new IllegalArgumentException("Don't know the namespace: " + value);
    }
}
