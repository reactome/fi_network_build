/*
 * Created on Mar 28, 2006
 *
 */
package org.reactome.psi;

public class Ref {
    private Class type;
    private String id;
    private Object container;
    private String propName;
    
    public Ref() {
    }

    public Object getContainer() {
        return container;
    }

    public void setContainer(Object container) {
        this.container = container;
    }

    public String getText() {
        return id;
    }
    
    public String getId() {
        return id;
    }
    
    /**
     * ref in BIND uses attribute "ref". Probably this is based an old schema.
     * @param id
     */
    public void setRef(String id) {
        this.id = id;
    }

    public void setText(String id) {
        this.id = id;
    }

    public String getPropName() {
        return propName;
    }

    public void setPropName(String propName) {
        this.propName = propName;
    }

    public Class getType() {
        return type;
    }

    public void setType(Class type) {
        this.type = type;
    }
    
    public String toString() {
        return container + ", " + id + ", " + propName + ", " + type;
    }
    
}