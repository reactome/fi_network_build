/*
 * Created on Mar 21, 2006
 *
 */
package org.reactome.fi.util;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;

public class HibernateUtil {

    private static Map<File, SessionFactory> configToSf = new HashMap<File, SessionFactory>();

    public static SessionFactory getSessionFactory(File configFile) {
        SessionFactory sf = configToSf.get(configFile);
        if (sf != null)
            return sf;
        if (sf == null) {
            try {
                sf = new Configuration().configure(configFile).buildSessionFactory();
            }
            catch(Exception e) {
                System.err.println("HibernateUtil.getSessionFactory(): " + e);
                throw new ExceptionInInitializerError(e);
            }
            configToSf.put(configFile, sf);
        }
        return sf;
    }
    
    public static SessionFactory getFISessionFactory() {
        File file = new File("resources/funcIntHibernate.cfg.xml");
        return getSessionFactory(file);
    }
    
}
