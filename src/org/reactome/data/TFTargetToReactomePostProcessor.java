package org.reactome.data;

import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.reactome.tred.TREDToReactomePostProcessor;

public class TFTargetToReactomePostProcessor extends TREDToReactomePostProcessor {
    @Override
    protected void attachDataSource(MySQLAdaptor dbAdaptor,
                                    XMLFileAdaptor fileAdaptor) throws Exception {
        super.attachDataSource("ENCODE", "", dbAdaptor, fileAdaptor);
    }
}
