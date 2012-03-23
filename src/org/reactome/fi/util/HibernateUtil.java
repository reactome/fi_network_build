/*
 * Created on Mar 21, 2006
 *
 */
package org.reactome.fi.util;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;

import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;

public class HibernateUtil {

    private static Map<File, SessionFactory> configToSf = new HashMap<File, SessionFactory>();

    /**
     * Get a SessionFactory configured by a configuration file.
     * @param configFile
     * @return
     */
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
    
    /**
     * Get a SessionFactory that can be used to re-generate database schema.
     * @param configFile
     * @return
     */
    public static SessionFactory getSessionFactoryWithCreate(File configFile) {
        try {
            Configuration configuration = new Configuration().configure(configFile);
            configuration.setProperty("hibernate.hbm2ddl.auto", "create");
            SessionFactory sf = configuration.buildSessionFactory();
            return sf;
        }
        catch(Exception e) {
            e.printStackTrace();
        }
        return null;
    }
    
}
