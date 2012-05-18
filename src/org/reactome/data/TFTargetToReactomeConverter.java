package org.reactome.data;

import java.util.Map;
import java.util.HashMap;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.persistence.MySQLAdaptor;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;

/**
 * Converter that converts (TF -> Target files) to Reactome Curator Tool projects
 *
 * @author Adrian Duong
 */
public class TFTargetToReactomeConverter {
	private XMLFileAdaptor fileAdaptor;
	private Map<String, GKInstance> nameToEntityMap;
	// private Map<String, GKInstance> interactionMap;
	private final GKInstance human;

	public TFTargetToReactomeConverter() throws Exception { // should this throw an exception?
		fileAdaptor = new XMLFileAdaptor();
		nameToEntityMap = new HashMap<String, GKInstance>();
		// interactionMap = new HashMap<String, GKInstance>();

		human = fileAdaptor.createNewInstance(ReactomeJavaConstants.Species);
		human.setAttributeValue(ReactomeJavaConstants.name, "Homo sapiens");
	}

    @Test
	public void convert() throws Exception {
		convert(FIConfiguration.getConfiguration().get("ENCODE_TFF_FILE"),
				FIConfiguration.getConfiguration().get("ENCODE_TFF_CONVERTED_FILE"));
	}

	private void convert(String tftargetFileName, String destFileName) throws Exception {
		MySQLAdaptor db = new MySQLAdaptor("localhost",
										   FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"),
										   FIConfiguration.getConfiguration().get("DB_USER"),
										   FIConfiguration.getConfiguration().get("DB_PWD"),
										   3306);

		// Map<String, Set<String>> gnToAccessMap = new UniProtAnalyzer.generateGeneNameToUniAccessMap(true);

		FileUtility fu = new FileUtility();
		Map<String, Set<String>> tfToTargetsMap = fu.loadSetMap(tftargetFileName);

		// make GKInstances of the interactions

		for (String tfName : tfToTargetsMap.keySet()) {
            GKInstance factor = nameToEntityMap.get(tfName);
            if (factor == null) {
                factor = fileAdaptor.createNewInstance(ReactomeJavaConstants.EntityWithAccessionedSequence);
                factor.addAttributeValue(ReactomeJavaConstants.name, tfName);
                factor.setDisplayName(tfName);

                nameToEntityMap.put(tfName, factor);
            }

			Set<String> targetNames = tfToTargetsMap.get(tfName);
			for (String targetName : targetNames) {
				GKInstance target = nameToEntityMap.get(targetName);
                if (target == null) {
                    target = fileAdaptor.createNewInstance(ReactomeJavaConstants.EntityWithAccessionedSequence);
                    target.addAttributeValue(ReactomeJavaConstants.name, targetName);
                    target.setDisplayName(targetName);

                    nameToEntityMap.put(targetName, target);
                }

				String interactionKey = factor.getDBID() + "-" + target.getDBID();
				GKInstance interaction = fileAdaptor.createNewInstance(ReactomeJavaConstants.TargettedInteraction);
				interaction.setAttributeValue(ReactomeJavaConstants.factor, factor);
				interaction.setAttributeValue(ReactomeJavaConstants.target, target);
				interaction.setAttributeValue(ReactomeJavaConstants.species, human);
				interaction.setDisplayName(factor.getDisplayName() + "-" + target.getDisplayName());
			}
		}

        TFTargetToReactomePostProcessor postProcessor = new TFTargetToReactomePostProcessor();
        postProcessor.postProcess(new MySQLAdaptor("localhost",
                                                   FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"),
                                                   FIConfiguration.getConfiguration().get("DB_USER"),
                                                   FIConfiguration.getConfiguration().get("DB_PWD"),
                                                   3306),
                                  fileAdaptor);
        

		fileAdaptor.save(destFileName);
	}
}