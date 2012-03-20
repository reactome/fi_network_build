/*
 * Created on Mar 25, 2006
 *
 */
package org.reactome.psi;

import static javax.xml.stream.XMLStreamConstants.CHARACTERS;
import static javax.xml.stream.XMLStreamConstants.END_DOCUMENT;
import static javax.xml.stream.XMLStreamConstants.END_ELEMENT;
import static javax.xml.stream.XMLStreamConstants.START_DOCUMENT;
import static javax.xml.stream.XMLStreamConstants.START_ELEMENT;

import java.io.FileInputStream;
import java.lang.reflect.Method;
import java.sql.Date;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.Characters;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import org.apache.log4j.Logger;

/**
 * This class is used to parse a PSI-MI xml file into the Java objects in this package. The parsing 
 * is based on new StAX XML parser. This is a pull based parser, and should be fast and with very 
 * small memory footprint.
 * @author guanming
 *
 */
public class PsiMiLoader {
    // A stack to keep the current converted object from the current element
    private Stack<Object> objectStack;
    // cache the method list to improve the performance
    private Map<Class, Method[]> clsMethods;
    private PsiMiModel model;
    private PsiMiPostParseProcessor postProcessor;
    // For logging
    static Logger logger = Logger.getLogger(PsiMiLoader.class);
    // For dereference
    private Map<String, Object> idMap;
    private List<Ref> references;
    
    public PsiMiLoader() {
        objectStack = new Stack<Object>();
        clsMethods = new HashMap<Class, Method[]>();
        model = new PsiMiModel();
        idMap = new HashMap<String, Object>();
        references = new ArrayList<Ref>();
    }
        
    public void setPostProcessor(PsiMiPostParseProcessor postProcessor) {
        this.postProcessor = postProcessor;
    }
    
    public void parse(String fileName) throws Exception {
        objectStack.clear(); // In case there is anything left. Should NOT be the case.
        XMLInputFactory inputFactory = XMLInputFactory.newInstance();
        inputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
        inputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
        inputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.TRUE);
        inputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.TRUE);
        FileInputStream fis = new FileInputStream(fileName);
        XMLEventReader eventReader = inputFactory.createXMLEventReader(fis);
        while (eventReader.hasNext()) {
            XMLEvent xmlEvent = eventReader.nextEvent();
            handleXMLEvent(xmlEvent);
        }
    }

    private void handleXMLEvent(XMLEvent event) throws Exception {
        int eventType = event.getEventType();
        switch (eventType) {
            case START_DOCUMENT:
                logger.info("Starting document loading...");
                break;
            case START_ELEMENT:
                StartElement startElement = event.asStartElement();
                handleStartElement(startElement);
                break;               
            case END_ELEMENT :
                handleEndElement(event.asEndElement());
                break;
            case CHARACTERS :
                handleCharacters(event.asCharacters());
                break;
            case END_DOCUMENT :
                logger.info("Finished document parsing!");
                break;
        }
    }
    
    private void handleEndElement(EndElement endElement) throws Exception{
        String elementName = endElement.getName().getLocalPart();
        if (!accept(elementName))
            return;
        Object obj = objectStack.pop(); // Just remove the last handled element (actually converted object).
        if (obj instanceof Entry) {
            postProcessEntry((Entry)obj);
        }
    }
    
    private void postProcessEntry(Entry entry) throws Exception {
        for (Ref ref : references) {
            if (ref.getType() == null)
                System.out.println("Null Type Ref: " + ref);
            Object propValue = idMap.get(ref.getType().getName() + "." + ref.getId());
            String propName = ref.getPropName();
            if (propName.endsWith("Ref")) {
                // remove "Ref"
                propName = propName.substring(0, propName.length() - 3);
            }
            setObjectProperty(ref.getContainer(), propName, propValue);
        }
        references.clear();
        idMap.clear();
        save(entry);
    }
    
    private void save(Object obj) throws Exception {
        if (postProcessor != null)
            postProcessor.postProcess(obj, model);
    }
    
    private void handleCharacters(Characters characters) throws Exception {
        if (characters.isWhiteSpace())
            return; // Nothing to process
        // Trim will remove new line characters too.
        String text = characters.getData().trim();
        // The previous should be an element name
        Object holder = objectStack.peek();
        if (holder.getClass().equals(String.class)) { // Common case
            String elementName = (String) objectStack.pop();
            Object container = objectStack.peek();
            setStringProperty(container, elementName, text);    
           // Add elementName back to the stack so that it can be removed in END_ELEMENT method call
            objectStack.push(elementName);
        }
        else {
            setStringProperty(holder, "text", text);
        }
    }
    
    private void setStringProperty(Object container, String propName, String propValue) throws Exception{
        if (propValue.equals("unassigned5"))
            System.out.println();
        // A special case to avoid handling parameter
        if (propName.equals(PsiMiModel.parameter))
            return;
        // Check if there is addMethod. If true, use add method since the property might be an array
        Method method = getAddMethod(container, propName);
        if (method == null)
            method = getSetMethod(container, propName);
        if (method == null) {
            logger.error("cannot find method for: " + container + " -> " + propName + " (" + propValue + ")");
            return;
        }
        Object paraValue = null;
        Class para = method.getParameterTypes()[0];
        // Only one
        if (para.equals(Integer.class) || para.equals(int.class)) // Cannot use "=" since int.class can be used
            paraValue = new Integer(propValue);
        else if (para.equals(Long.class) || para.equals(long.class))
            paraValue = new Long(propValue);
        else if (para.equals(Date.class)) {
            if (propValue.length() > 10) // Should be in the format: yyyy-MM-dd
                propValue = propValue.substring(0, 10); // Avoid considering time-zone.
//            if (propValue.endsWith("Z")) { // UTC value used by IntAct
//                propValue = propValue.substring(0, propValue.length() - 1);
//            }
            paraValue = Date.valueOf(propValue);
        }
        else if (para.equals(Boolean.class) || para.equals(boolean.class))
            paraValue = Boolean.valueOf(propValue);
        else
            paraValue = propValue;
        logger.info("Set Property: " + propName + " for " + container.getClass() + " with value " + propValue);
        method.invoke(container, new Object[]{paraValue});
        if (propName.equals("id")) {
            idMap.put(container.getClass().getName() + "." + propValue, container);
        }
    }
    
    private Method getAddMethod(Object obj, String propName) {
        if (propName.endsWith("List"))
            propName = propName.substring(0, propName.length() - 4); // Remove "List"
        String addMethodName = "add" + upperFirst(propName);
        Method[] methods = clsMethods.get(obj.getClass());
        if (methods == null) {
            methods = obj.getClass().getMethods();
            clsMethods.put(obj.getClass(), methods);
        }
        for (Method m : methods) {
            if (m.getName().equals(addMethodName))
                return m;
        }
        return null;
    }
    
    private Method getSetMethod(Object obj, String propName) {
        // for old schema
        if (propName.endsWith("Interactor") &&
            !propName.equals(PsiMiModel.experimentalInteractor))
            propName = "interactor";
        String setMethodName = "set" + upperFirst(propName);
        Method[] methods = clsMethods.get(obj.getClass());
        if (methods == null) {
            methods = obj.getClass().getMethods();
            clsMethods.put(obj.getClass(), methods);
        }
        for (Method m : methods) {
            if (m.getName().equals(setMethodName))
                return m;
        }
        return null;
    }
    
    private boolean accept(String elmName) {
        if (elmName.equals(PsiMiModel.entrySet))
            return false;
        // This is a special case
        if (elmName.equals(PsiMiModel.participant)) {
            // Peel the parent
            Object container = objectStack.peek();
            if (container instanceof InferredInteraction) {
                // Escape it
                return false;
            }
        }
        // Another special case: experimentRef in parameter should not be handled since parameter is not handled
        if (elmName.equals(PsiMiModel.experimentRef)) {
            Object container = objectStack.peek();
            if (container instanceof String && container.equals(PsiMiModel.parameter))
                return false;
        }
        return true;
    }
    
    private void handleStartElement(StartElement startElement) throws Exception {
        String elmName = startElement.getName().getLocalPart();
        if (!accept(elmName))
            return;
        Class cls = getClassFromElementName(elmName);
        if (cls == String.class) {
            objectStack.push(elmName);
        }
        else {
            Object obj = cls.newInstance();
            handleAttributes(startElement, obj);
            if (obj instanceof Interactor)
                handleInteractorType((Interactor)obj, elmName);
            if (!objectStack.isEmpty()) {
                Object container = objectStack.peek();
                setObjectProperty(container, elmName, obj);
            }
            objectStack.push(obj);
        }
    }
    
    private void handleInteractorType(Interactor interactor, String elmName) {
        if (elmName.equals(PsiMiModel.interactor))
            return; // Should be handled 
        // For BIND, type is defined in the elementname. Probably it will be updated soon.
        OpenCV interactorType = model.getInteractorType(elmName);
        interactor.setInteractorType(interactorType);
    }
    
    private void handleAttributes(StartElement startElement, Object object) throws Exception {
        for (Iterator it = startElement.getAttributes(); it.hasNext();) {
            Attribute attribute = (Attribute) it.next();
            String attName = attribute.getName().getLocalPart();
            String value = attribute.getValue();
            setStringProperty(object, attName, value);
        }
    }
    
    private String upperFirst(String name) {
        return name.substring(0, 1).toUpperCase() + name.substring(1);
    }
    
    private Class getClassFromElementName(String elementName) {
        // Have to do a little process
        if (elementName.endsWith("Interactor") &&
            !elementName.equals(PsiMiModel.experimentalInteractor))
            elementName = "interactor";
        else if (elementName.endsWith("Participant"))
            elementName = "participant";
        Class cls = model.getElementType(elementName);
        if (cls == null)
            cls = String.class;
        return cls;
    }
    
    private Class getTypeForRef(String elementName) {
        String typeName = elementName.substring(0, elementName.length() - 3); // It ends with ref
        return model.getElementType(typeName);
    }
    
    private void setObjectProperty(Object container, String propName, Object propValue) throws Exception {
        String oldPropName = propName;
        if (container instanceof String) {
            // It should be the one before container
            Object tmp = objectStack.pop();
            container = objectStack.peek();
            objectStack.push(tmp); // push it back
            // Have to change the propName
            propName = tmp.toString();
        }
        if (propValue instanceof Ref) {
            Ref ref = (Ref) propValue;
            ref.setContainer(container);
            ref.setPropName(propName);
            ref.setType(getTypeForRef(oldPropName));
            references.add(ref);
            return;
        }
        Method method = getAddMethod(container, propName);
        if (method == null)
            method = getSetMethod(container, propName);
        if (method == null) {
            logger.error("Cannot find method for: " + container + " <- " + propName);
            return;
        }
        logger.info("invoke " + method + " with value \"" + propValue + "\" on " + container);
        // A special case
        if (container.getClass() == Experiment.class &&
            propValue.getClass() == OpenCVExperimentalWrapper.class)
            propValue = ((OpenCVExperimentalWrapper)propValue).getOpenCV();
        //System.out.println("invoke: " + container + " <- " + propValue + " for " + propName);
        method.invoke(container, new Object[]{propValue});
    }
}
