/*
 * Created on Jan 31, 2008
 *
 */
package org.reactome.weka;

import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;

/**
 * This class is used to analyze the training positive data set
 * @author guanming
 *
 */
public class TrainingDataAnalyzer {
    
    @Test
    public void getDistributions() throws IOException {
        String arffName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "FourDBInteractions013108.arff";
        FileUtility fu = new FileUtility();
        fu.setInput(arffName);
        String line = fu.readLine();
        Map<String, Integer> typeToCount = new HashMap<String, Integer>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("true,")) {
                // Count the positive dataset only
                Integer c = typeToCount.get(line);
                if (c == null)
                    typeToCount.put(line, 1);
                else
                    typeToCount.put(line, ++c);
            }
        }
        fu.close();
        // Print out
        int total = 0;
        for (Iterator<String> it = typeToCount.keySet().iterator(); it.hasNext();) {
            String type = it.next();
            Integer c = typeToCount.get(type);
            System.out.println(type + ": " + c);
            total += c;
        }
        System.out.println("Total: " + total);
    }
    
}
