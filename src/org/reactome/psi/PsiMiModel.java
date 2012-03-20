/*
 * Created on Mar 25, 2006
 *
 */
package org.reactome.psi;

import java.beans.BeanInfo;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.hibernate.Session;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;
import org.jdom.xpath.XPath;

public class PsiMiModel {
    
    public static final String alias = "alias";
    public static final String attribute = "attribute";
    public static final String attributeList = "attributeList";
    public static final String availability = "availability";
    public static final String availabilityList = "availabilityList";
    public static final String availabilityRef = "availabilityRef";
    public static final String base = "base";
    public static final String begin = "begin";
    public static final String beginInterval = "beginInterval";
    public static final String bibref = "bibref";
    public static final String biologicalRole = "biologicalRole";
    public static final String cellType = "cellType";
    public static final String compartment = "compartment";
    public static final String confidence = "confidence";
    public static final String confidenceList = "confidenceList";
    public static final String db = "db";
    public static final String dbAc = "dbAc";
    public static final String end = "end";
    public static final String endInterval = "endInterval";
    public static final String endStatus = "endStatus";
    public static final String entry = "entry";
    public static final String entrySet = "entrySet";
    public static final String experimentDescription = "experimentDescription";
    public static final String experimentList = "experimentList";
    public static final String experimentRef = "experimentRef";
    public static final String experimentRefList = "experimentRefList";
    public static final String experimentalInteractor = "experimentalInteractor";
    public static final String experimentalInteractorList = "experimentalInteractorList";
    public static final String experimentalPreparation = "experimentalPreparation";
    public static final String experimentalPreparationList = "experimentalPreparationList";
    public static final String experimentalRole = "experimentalRole";
    public static final String experimentalRoleList = "experimentalRoleList";
    public static final String exponent = "exponent";
    public static final String factor = "factor";
    public static final String feature = "feature";
    public static final String featureDetectionMethod = "featureDetectionMethod";
    public static final String featureList = "featureList";
    public static final String featureRange = "featureRange";
    public static final String featureRangeList = "featureRangeList";
    public static final String featureType = "featureType";
    public static final String fullName = "fullName";
    public static final String hostOrganism = "hostOrganism";
    public static final String hostOrganismList = "hostOrganismList";
    public static final String id = "id";
    public static final String imexId = "imexId";
    public static final String inferredInteraction = "inferredInteraction";
    public static final String inferredInteractionList = "inferredInteractionList";
    public static final String interaction = "interaction";
    public static final String interactionDetectionMethod = "interactionDetectionMethod";
    public static final String interactionList = "interactionList";
    public static final String interactionRef = "interactionRef";
    public static final String interactionType = "interactionType";
    public static final String interactor = "interactor";
    public static final String interactorList = "interactorList";
    public static final String interactorRef = "interactorRef";
    public static final String interactorType = "interactorType";
    public static final String intraMolecular = "intraMolecular";
    public static final String isLink = "isLink";
    public static final String level = "level";
    public static final String minorVersion = "minorVersion";
    public static final String modelled = "modelled";
    public static final String name = "name";
    public static final String nameAc = "nameAc";
    public static final String names = "names";
    public static final String ncbiTaxId = "ncbiTaxId";
    public static final String negative = "negative";
    public static final String organism = "organism";
    public static final String parameter = "parameter";
    public static final String parameterList = "parameterList";
    public static final String participant = "participant";
    public static final String participantFeatureRef = "participantFeatureRef";
    public static final String participantIdentificationMethod = "participantIdentificationMethod";
    public static final String participantIdentificationMethodList = "participantIdentificationMethodList";
    public static final String participantList = "participantList";
    public static final String participantRef = "participantRef";
    public static final String position = "position";
    public static final String primaryRef = "primaryRef";
    public static final String refType = "refType";
    public static final String refTypeAc = "refTypeAc";
    public static final String release = "release";
    public static final String releaseDate = "releaseDate";
    public static final String secondary = "secondary";
    public static final String secondaryRef = "secondaryRef";
    public static final String sequence = "sequence";
    public static final String shortLabel = "shortLabel";
    public static final String source = "source";
    public static final String startStatus = "startStatus";
    public static final String term = "term";
    public static final String termAc = "termAc";
    public static final String tissue = "tissue";
    public static final String type = "type";
    public static final String typeAc = "typeAc";
    public static final String unit = "unit";
    public static final String unitAc = "unitAc";
    public static final String value = "value";
    public static final String version = "version";
    public static final String xref = "xref";
    
    private Map<String, Class> elementTypeMap;
    private Map<Class, Set<String>> beanProperties;
    private Map<String, OpenCV> interactorTypeMap;
    
    public PsiMiModel() {
        init();
    }
    
    public OpenCV getInteractorType(String elementName) {
        if (interactorTypeMap == null) {
            interactorTypeMap = new HashMap<String, OpenCV>();
        }
        OpenCV type = interactorTypeMap.get(elementName);
        if (type == null) {
            type = new OpenCV();
            String typeName = elementName.substring(0, elementName.length() - 10); // 10 for "interactor"
            type = new OpenCV();
            Names names = new Names();
            names.setShortLabel(typeName);
            type.setNames(names);
            interactorTypeMap.put(elementName, type);
        }
        return type;
    }
    
    public void prepareModelForSave(Session session) throws Exception {
        if (interactorTypeMap != null) {
            Collection<OpenCV> types = interactorTypeMap.values();
            for (OpenCV cv : types) {
                session.saveOrUpdate(cv); // re-attach to a session so that
                                          // other objects referring to them
                                          // can be persisted.
            }
        }
    }
    
    private void init() {
        String mapFile = "resources" + File.separator + "ElementMap.txt";
        elementTypeMap = new HashMap<String, Class>();
        try {
            FileReader fileReader = new FileReader(mapFile);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            String line = null;
            int index = 0;
            String elmName = null;
            String typeName = null;
            while ((line = bufferedReader.readLine()) !=  null) {
                index = line.indexOf("->");
                elmName = line.substring(0, index);
                typeName = line.substring(index + 2);
                if (typeName.startsWith("java")) // will be handled later
                    elementTypeMap.put(elmName, null);
                else
                    elementTypeMap.put(elmName, Class.forName("org.reactome.psi." + typeName));
            }
            loadBeanProperties();
        }
        catch(Exception e) {
            System.err.println("PsiMiConstants.init(): " + e);
            e.printStackTrace();
        }
    }
    
    private void loadBeanProperties() {
        beanProperties = new HashMap<Class, Set<String>>();
        Set<String> keys = elementTypeMap.keySet();
        try {
            for (String elmName : keys) {
                Class type = elementTypeMap.get(elmName);
                if (type == null)
                    continue;
                BeanInfo beanInfo = Introspector.getBeanInfo(type);
                PropertyDescriptor[] properties = beanInfo.getPropertyDescriptors();
                Set<String> propNames = new HashSet<String>();
                for (PropertyDescriptor pd : properties) {
                    propNames.add(pd.getName());
                }
                beanProperties.put(type, propNames);
            }
        }
        catch(Exception e) {
            System.err.println("PsiMiModel.loadBeanProperties(): " + e);
            e.printStackTrace();
        }
    }
    
    public boolean containProperty(String typeName, String propName) {
        Set<String> propNames = beanProperties.get(typeName);
        return propName.contains(propName);
    }
    
    /*
     * A partial implementation to map element types to Java types. This mapping is not complete and
     * should NOT be used on the fly in any application.
     */
    public Map<String, String> setPsiMiSchema(String schemaFileName) {
        Map<String, String> elementTypeNameMap = new HashMap<String, String>();
        try {
            SAXBuilder builder = new SAXBuilder();
            Document document = builder.build(schemaFileName);
            Element root = document.getRootElement();
            String xpath = ".//xs:element";
            List elements = XPath.selectNodes(root, xpath);
            Element elm = null;
            for (Iterator it = elements.iterator(); it.hasNext();) {
                elm = (Element) it.next();
                String name = elm.getAttributeValue("name");
                String type = elm.getAttributeValue("type");
                type = cleanupType(name, type);
                elementTypeNameMap.put(name, type);
            }
            // This is for the old schema
            elementTypeNameMap.put("proteinParticipant", "Participant");
            elementTypeNameMap.put("proteinInteractor", "Interactor");
        }
        catch(Exception e) {
            System.err.println("PsiMiContants.setPsiMiSchema(): " + e);
            e.printStackTrace();
        }
        return elementTypeNameMap;
    }
    
    public Class getElementType(String elementName) {
        if (elementTypeMap != null)
            return elementTypeMap.get(elementName);
        return null;
    }
    
    public Map<String, Class> getElementTypeMap() {
        return elementTypeMap;
    }
    
    private String cleanupType(String name, String type) {
        if (type == null)
            type = name;
        if (type.endsWith("Type")) {
            type = type.substring(0, type.length() - 4);
        }
        String rtn = type.substring(0, 1).toUpperCase() + type.substring(1);
        if (rtn.equals("Cv"))
            rtn = "OpenCV";
        else if (rtn.equals("OpenCv"))
            rtn = "OpenCV";
        else if (rtn.equals("Xs:boolean"))
            rtn = "java.lang.Boolean";
        else if (rtn.equals("Xs:string"))
            rtn = "java.lang.String";
        return rtn;
    }
    
    
    public static void main(String[] args) {
        // Used to generate constants from MIF25.xsd
        String schemaFileName = "/Users/wgm/Documents/caBIG_R3/datasets/PSI-MI/MIF25.xsd";
        String regex = "xs:[:element|attribute].*name=\"(\\w*)\"";
        Pattern pattern = Pattern.compile(regex);
        try {
            FileReader fileReader = new FileReader(schemaFileName);
            BufferedReader reader = new BufferedReader(fileReader);
            String line = null;
            Set<String> nameSet = new HashSet<String>();
            while ((line = reader.readLine()) != null) {
                Matcher matcher = pattern.matcher(line);
                if (matcher.find()) {
                    String name = matcher.group(1);
                    if (name.equals("cvType"))
                        System.out.println("Line: " + line);
                    nameSet.add(matcher.group(1));
                }
            }
            System.out.println("Total Names: " + nameSet.size());
            List<String> nameList = new ArrayList<String>(nameSet);
            Collections.sort(nameList);
            for (String name : nameList) {
                System.out.printf("    public static final String %s = \"%s\";%n", name, name);
            }
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }
}
