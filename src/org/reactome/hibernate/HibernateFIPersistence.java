/*
 * Created on Sep 29, 2006
 *
 */
package org.reactome.hibernate;

import java.io.File;

import org.apache.log4j.PropertyConfigurator;
import org.hibernate.SessionFactory;
import org.reactome.fi.util.HibernateUtil;

/**
 * Hibernate based database analysis class for Functional interactions.
 * @author guanming
 */
public class HibernateFIPersistence {
    
    protected SessionFactory sessionFactory;
    
    protected HibernateFIPersistence() {     
        init();
    }
    
    /**
     * For system set up
     */
    private void init() {
        try {
            PropertyConfigurator.configure("resources/log4j.properties");
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }
    
    public SessionFactory initSession() throws Exception {
        String configFileName = "resources/funcIntHibernate.cfg.xml";
        File configFile = new File(configFileName);
        sessionFactory = HibernateUtil.getSessionFactory(configFile);
        return sessionFactory;
    }
    
    public SessionFactory getSessionFactory() {
        return sessionFactory;
    }
    
    public void setSessionFactor(SessionFactory sf) {
        this.sessionFactory = sf;
    }
    
}
