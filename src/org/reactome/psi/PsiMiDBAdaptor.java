/*
 * Created on Mar 28, 2006
 *
 */
package org.reactome.psi;

import org.apache.log4j.Logger;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.Transaction;

/**
 * This class is used to talk to a PSI-MI database based on hibernate.
 * @author guanming
 *
 */
public class PsiMiDBAdaptor implements PsiMiPostParseProcessor {
    private SessionFactory sessionFactory;
    private Logger logger = Logger.getLogger(PsiMiDBAdaptor.class);
    
    public PsiMiDBAdaptor() {
    }
    
    public PsiMiDBAdaptor(SessionFactory sessionFactory) {
        this();
        setSessionFactory(sessionFactory);
    }
    
    public void postProcess(Object obj, PsiMiModel model) throws Exception {
        save(obj, model);
    }
    
    public void setSessionFactory(SessionFactory sessionFactory) {
        this.sessionFactory = sessionFactory;
    }
    
    public SessionFactory getSessionFactory() {
        return this.sessionFactory;
    }
    
    public void saveInSession(Object obj) throws Exception {
        logger.info("Save object: " + obj);
        Session session = sessionFactory.getCurrentSession();
        session.persist(obj);
    }
    
    public Session getCurrentSession() {
        return sessionFactory.getCurrentSession();
    }
    
    public void save(Object obj, PsiMiModel model) throws Exception {
        logger.info("Save object: " + obj);
        Session session = sessionFactory.getCurrentSession();
        Transaction tx = null;
        try {
            tx = session.beginTransaction();
            if (model != null)
                model.prepareModelForSave(session);
            session.persist(obj);
            tx.commit();
            if (session.isOpen())
                session.close();
        }
        catch(Exception e) {
            if (tx != null)
                tx.rollback();
            if (session.isOpen())
                session.close();
            throw e; // rethrow the exception to be caught by the caller.
        }
    }

}
