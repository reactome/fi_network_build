/*
 * Created on May 17, 2006
 *
 */
package org.reactome.fi.util;

import java.io.*;
import java.nio.channels.FileChannel;
import java.util.*;

/**
 * A utility class to consolidate read/write activities related to simple text processing.
 * @author guanming
 */
public class FileUtility {
    private FileReader fileReader;
    private BufferedReader bufferedReader;
    private FileWriter fileWriter;
    private PrintWriter printWriter;
    
    public FileUtility() {
    }
    
    public PrintStream createPrintStream(String fileName) throws IOException {
        FileOutputStream fos = new FileOutputStream(fileName);
        PrintStream ps = new PrintStream(fos);
        return ps;
    }
    
    public void setInput(String fileName) throws IOException {
        fileReader = new FileReader(fileName);
        bufferedReader = new BufferedReader(fileReader);
    }
    
    public void setOutput(String fileName) throws IOException {
        fileWriter = new FileWriter(fileName);
        printWriter = new PrintWriter(fileWriter);
    }
    
    public PrintWriter getPrintWriter() {
        return printWriter;
    }
    
    public String readLine() throws IOException {
        return bufferedReader.readLine();
    }
    
    public void printLine(String line) throws IOException {
        printWriter.println(line);
    }
    
    public void close() throws IOException {
        if (bufferedReader != null)
            bufferedReader.close();
        if (fileReader != null)
            fileReader.close();
        if (printWriter != null)
            printWriter.close();
        if (fileWriter != null)
            fileWriter.close();
    }
    
    public void exportMap(Map<String, String> map, String fileName) throws IOException {
        setOutput(fileName);
        for (Iterator<String> it = map.keySet().iterator(); it.hasNext();) {
            String key = it.next();
            String value = map.get(key);
            printLine(key + "\t" + value);
        }
        close();
    }
    
    public Map<String, String> importMap(String fileName) throws IOException {
        Map<String, String> map = new HashMap<String, String>();
        setInput(fileName);
        String line = null;
        int index = 0;
        while ((line = readLine()) != null) {
            index = line.indexOf("\t");
            map.put(line.substring(0, index), line.substring(index + 1));
        }
        close();
        return map;
    }
    
    public Set<String> loadInteractions(String interactionFile) throws IOException {
        Set<String> interactions = new HashSet<String>();
        setInput(interactionFile);
        String line = null;
        while ((line = readLine()) != null)
            interactions.add(line);
        close();
        return interactions;
    }
    
    public void outputSet(Set<String> set, String fileName) throws IOException {
        setOutput(fileName);
        for (String tmp : set)
            printLine(tmp);
        close();
    }
    
    public Set<String> loadSet(String fileName) throws IOException {
        setInput(fileName);
        Set<String> set = new HashSet<String>();
        String line = null;
        while ((line = readLine()) != null)
            set.add(line);
        close();
        return set;
    }
    
    public void saveInteractions(Set<String> interactions, String fileName) throws IOException {
        saveCollection(interactions, 
                             fileName);
    }
    
    /**
     * Save a String typed values in any Collection. The order will be set by the
     * client to this method.
     * @param values
     * @param fileName
     * @throws IOException
     */
    public void saveCollection(Collection<String> values, 
                               String fileName) throws IOException {
        setOutput(fileName);
        for (String value : values)
            printLine(value);
        close();
    }
    
    public Object loadObject(String fileName) throws IOException, ClassNotFoundException {
        FileInputStream fis = new FileInputStream(fileName);
        ObjectInputStream ois = new ObjectInputStream(fis);
        Object rtn = ois.readObject();
        ois.close();
        fis.close();
        return rtn;
    }
    
    public void saveObjet(Object obj, String fileName) throws Exception {
        FileOutputStream fos = new FileOutputStream(fileName);
        ObjectOutputStream oos = new ObjectOutputStream(fos);
        oos.writeObject(obj);
        oos.close();
        fos.close();
    }

    public void saveSetMap(Map<String, Set<String>> setMap, 
                           String outputFileName) throws IOException {
        setOutput(outputFileName);
        for (Iterator<String> it = setMap.keySet().iterator(); it.hasNext();) {
            String key = it.next();
            Set<String> set = setMap.get(key);
            for (String tmp : set) 
                printLine(key + "\t" + tmp);
        }
        close();
    }

    
    public void saveSetMapInSort(Map<String, Set<String>> setMap, 
                                 String outputFileName) throws IOException {
        Set<String> keys = setMap.keySet();
        List<String> keyList = new ArrayList<String>(keys);
        Collections.sort(keyList);
        setOutput(outputFileName);
        for (String key : keyList) {
            Set<String> set = setMap.get(key);
            List<String> list = new ArrayList<String>(set);
            Collections.sort(list);
            for (String tmp : list) 
                printLine(key + "\t" + tmp);
        }
        close();
    }
    
    public Map<String, Set<String>> loadSetMap(String fileName) throws IOException {
        Map<String, Set<String>> setMap = new HashMap<String, Set<String>>();
        setInput(fileName);
        String line = null;
        int index = 0;
        Set<String> set = null;
        String term1, term2;
        while ((line = readLine()) != null) {
            index = line.indexOf("\t");
            term1 = line.substring(0, index);
            term2 = line.substring(index + 1);
            set = setMap.get(term1);
            if (set == null) {
                set = new HashSet<String>();
                setMap.put(term1, set);
            }
            set.add(term2);
        }
        close();
        return setMap;
    }
      public Set<String> loadInteractions(String interactionFile, int numberOfPairs) throws IOException {
        Set<String> interactions = new HashSet<String>();
        setInput(interactionFile);
        String line = null;
        while ((line = readLine()) != null && numberOfPairs-->0)
            interactions.add(line);
        close();
        return interactions;
    }
      
      /*
       * Load a set of interactions from a sif file.
       * 
       */
      public Set<String> loadInteractionsInSifFile(String fileName) throws IOException {
          setInput(fileName);
          String line = null;
          Set<String> interactions = new HashSet<String>();
          while ((line = readLine()) != null) {
              String[] tokens = line.split("\t");
              String gene1 = tokens[0];
              String gene2 = tokens[2];
              int compare = gene1.compareTo(gene2);
              if (compare < 0)
                  interactions.add(gene1 + "\t" + gene2);
              else
                  interactions.add(gene2 + "\t" + gene1);
          }
          close();
          return interactions;
      }
      
      /*
       * Self interactions stand for those interactions that contain two same ids
       */
      public Set<String> removeSelfInteractions(Set<String> interactionSet,String delimiter) throws Exception {
    	  Set<String> processedInteractionsSet = new HashSet<String>();
    	  for(String interaction:interactionSet) {
    		  int index = interaction.indexOf(delimiter);
    		  String id1 = interaction.substring(0,index);
    		  String id2 = interaction.substring(index+1);
    		  if( !(id1.equals(id2)) ){
    			  processedInteractionsSet.add(id1+delimiter+id2);
    		  }
    	  }
    	  return processedInteractionsSet;
      }
      public Map<String,String> loadInteractionPairs(String interactionFile, String delimiter) throws IOException {
    	Map<String,String> pairs = new HashMap<String,String>();
    	setInput(interactionFile);
    	String line = null;
    	int index = 0;
    	String term1, term2;
    	while ((line = readLine()) != null) {
    		index = line.indexOf(delimiter);
    		term1 = line.substring(0, index);
    		term2 = line.substring(index + 1);
    		pairs.put(term1, term2);
    	}
    	return pairs;
    }
    
    public Set<String>  loadInteractionTerms(String interactionFile, String delimiter) throws IOException {
    	Set<String> termSet = new HashSet<String>();
    	setInput(interactionFile);
    	String line = null;
    	int index = 0;
    	String term1, term2;
    	while ((line = readLine()) != null) {
    		index = line.indexOf(delimiter);
    		term1 = line.substring(0, index);
    		term2 = line.substring(index + 1);
    		termSet.add(term1);
    		termSet.add(term2);
    	}
    	return termSet;
    }
    
    public Set<String> loadInteractionTerms(String fileName, String delimiter, boolean useString1) throws Exception {
		FileUtility fu = new FileUtility();
		Set<String> pairSet = fu.loadInteractions(fileName);
		Set<String> itemSet = new HashSet<String>();
		
		if(useString1)
			for(String pair:pairSet) {
				int index = 0;
				index = pair.indexOf(delimiter);
				itemSet.add(pair.substring(0,index));
			}
		else
			for(String pair:pairSet) {
				int index = 0;
				index = pair.indexOf(delimiter);
				itemSet.add(pair.substring(index+1));
			}
		
		return pairSet;
	}
    public Map<String,String> loadInteractionPairs(String interactionFile, String delimiter, boolean useTerm1AsKey) throws IOException {
    	if(useTerm1AsKey) {
    		return loadInteractionPairs(interactionFile, delimiter);
    	}
    	else {
    		Map<String,String> pairs = new HashMap<String,String>();
        	setInput(interactionFile);
        	String line = null;
        	int index = 0;
        	String term1, term2;
        	while ((line = readLine()) != null) {
        		index = line.indexOf(delimiter);
        		term1 = line.substring(0, index);
        		term2 = line.substring(index + 1);
        		pairs.put(term2, term1);
        	}
        	return pairs;
    	}
    }
    
    public Map<String,String> loadInteractionPairs(String interactionFile, 
    													   String delimiter, 
    													   boolean useTerm1AsKey,
    													   boolean removeSelfPair) throws IOException {
    	String term1, term2;
    	Map<String,String> pairs = new HashMap<String,String>();
    	String line = null;
    	int index = 0;
    	if(useTerm1AsKey) {
    		if (!removeSelfPair) {
    			return loadInteractionPairs(interactionFile, delimiter);
    		}   
    		else {
    			setInput(interactionFile);
    			while ((line = readLine()) != null) {
    				index = line.indexOf(delimiter);
    				term1 = line.substring(0, index);
    				term2 = line.substring(index + 1);
    				//If P00533 is mapped as P00533, remove this mapping pair
    				if(term1 == term2)
    					continue;
    				else
    					pairs.put(term1, term2);
    			}
    		}
    	}
    	else {
    		if(removeSelfPair) {
    			setInput(interactionFile);
    			while ((line = readLine()) != null) {
    				index = line.indexOf(delimiter);
    				term1 = line.substring(0, index);
    				term2 = line.substring(index + 1);
    				//If P00533 is mapped as P00533, remove this mapping pair
    				if(term1 == term2)
    					continue;
    				else
    					pairs.put(term2, term1);
    			}
    		}
    		else {
    			setInput(interactionFile);
            	while ((line = readLine()) != null) {
            		index = line.indexOf(delimiter);
            		term1 = line.substring(0, index);
            		term2 = line.substring(index + 1);
            		pairs.put(term2, term1);
            	}
    		}
    	}
    	return pairs;
    	
    }
    	
    public Map<String, Set<String>> loadSetMap(String fileName, String delimiter) throws IOException {
        Map<String, Set<String>> setMap = new HashMap<String, Set<String>>();
        setInput(fileName);
        String line = null;
        int index = 0;
        Set<String> set = null;
        String term1, term2;
        while ((line = readLine()) != null) {
            index = line.indexOf(delimiter);
            term1 = line.substring(0, index);
            term2 = line.substring(index + 1);
            set = setMap.get(term1);
            if (set == null) {
                set = new HashSet<String>();
                setMap.put(term1, set);
            }
            set.add(term2);
        }
        close();
        return setMap;
    }
    
    public Map<String, Set<String>> loadSetMap(String fileName, String delimiter,boolean useTerm1AsKey) throws IOException {
        if (useTerm1AsKey){
            return loadSetMap(fileName,delimiter);
        }
        else{
            Map<String, Set<String>> setMap = new HashMap<String, Set<String>>();
            setInput(fileName);
            String line = null;
            int index = 0;
            Set<String> set = null;
            String term1, term2;
            while ((line = readLine()) != null) {
                index = line.indexOf(delimiter);
                term1 = line.substring(0, index);
                term2 = line.substring(index + 1);
                set = setMap.get(term2);
                if (set == null) {
                    set = new HashSet<String>();
                    setMap.put(term2, set);
                }
                set.add(term1);
            }
            close();
            return setMap;
        }
    }
    
    public ArrayList<HashSet<String>> loadClusters(String fileName, int minSize) throws IOException {
        setInput(fileName);
        ArrayList<HashSet<String>> set = new ArrayList<HashSet<String>>();
        String line = null;
        while ((line = readLine()) != null) {
        	List<String> list = Arrays.asList(line.split("\t"));
        	if (list.size() >= minSize) {
        		HashSet<String>tmp = new HashSet<String>(list);
            	set.add(tmp);
        	}
       }
        close();
        return set;
    }
    
    /**
     * Copy a file.
     * @param src
     * @param dest
     * @throws IOException
     */
    public void copy(File src, File dest) throws IOException {
        if (!src.exists())
            throw new IllegalArgumentException(src.getAbsolutePath() + " doesn't exist for copy!");
        if (!dest.exists())
            dest.createNewFile();
        FileChannel sourceChannel = new FileInputStream(src).getChannel();
        FileChannel destChannel = new FileOutputStream(dest).getChannel();
        destChannel.transferFrom(sourceChannel, 0, sourceChannel.size());
        sourceChannel.close();
        destChannel.close();
    }
    
}
