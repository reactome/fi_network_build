/*
 * Created on Jun 27, 2006
 *
 */
package org.reactome.weka;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Enumeration;

import javax.swing.BorderFactory;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;

import weka.classifiers.Classifier;
import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;

public class WEKAClassifierGUITest {
    
    /**
     * This method is used to launch an GUI to check predicated values
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {
        FileUtility fu = new FileUtility();
        //Classifier classifier = (Classifier) fu.loadObject(R3Constants.NBC_MODEL_NAME);
        //Classifier classifier = (Classifier) fu.loadObject("results/NoMFDisc20NaiveBayes091506.model");
        //Classifier classifier = (Classifier) fu.loadObject("results/NaiveBayes091106.model");
        //Classifier classifier = (Classifier) fu.loadObject("results/NaiveBayesCombined.model");
        Classifier classifier = (Classifier) fu.loadObject(FIConfiguration.getConfiguration().get("RESULT_DIR") + "NBC_Pos_Neg_02.model");
        JFrame frame = new JFrame("Check Instances");
        WEKAResultAnalyzer resultAnalyzer = new WEKAResultAnalyzer();
        initPanel(frame, resultAnalyzer.createDataSetForV3(), classifier);
        frame.setLocation(465, 80);
        frame.setSize(425, 550);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
    }
    
    private static void initPanel(JFrame frame, 
                                  final Instances dataset, 
                                  final Classifier bayesNet) {
        JPanel predicatorPane = new JPanel();
        predicatorPane.setBorder(BorderFactory.createEtchedBorder());
        predicatorPane.setLayout(new GridBagLayout());
        GridBagConstraints constraints = new GridBagConstraints();
        constraints.anchor = GridBagConstraints.WEST;
        constraints.fill = GridBagConstraints.HORIZONTAL;
        Attribute attribute = null;
        Enumeration attEnum = dataset.enumerateAttributes();
        int y = 0;
        // Values label
        final JLabel firstValueLabel = new JLabel();
        final JLabel secondValueLabel = new JLabel();
        final JComponent[] valueBoxes = new JComponent[dataset.numAttributes()];
        ActionListener l = new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Instance instance = new Instance(dataset.numAttributes());
                instance.setDataset(dataset);
                instance.setClassMissing();
                for (int i = 0; i < valueBoxes.length - 1; i++) {
                    if (valueBoxes[i] instanceof JComboBox) {
                        JComboBox box = (JComboBox) valueBoxes[i];
                        //if (box.getSelectedIndex() == 0)
                        //    continue;
                        instance.setValue(i + 1, box.getSelectedIndex());
                    }
                    else if (valueBoxes[i] instanceof JTextField) {
                        JTextField tf = (JTextField) valueBoxes[i];
                        String text = tf.getText().trim();
                        if (text.length() == 0 || text.equalsIgnoreCase("unknown"))
                            continue;
                        //instance.setValue(i + 1, Double.parseDouble(text));
                        instance.setValue(i + 1, Integer.parseInt(text));
                    }
                }
                try {
                    double[] prob = bayesNet.distributionForInstance(instance);
                    firstValueLabel.setText(prob[0] + "");
                    secondValueLabel.setText(prob[1] + "");
                }
                catch(Exception e1) {
                    e1.printStackTrace();
                }
            }
        };
        int index = 0;
        while (attEnum.hasMoreElements()) {
            attribute = (Attribute) attEnum.nextElement();
            constraints.gridy = y;
            constraints.gridx = 0;
            constraints.insets = new Insets(4, 4, 4, 4);
            predicatorPane.add(new JLabel(attribute.name()), constraints);
            if (attribute.isNominal()) {
                JComboBox valueBox = new JComboBox();
                valueBoxes[index] = valueBox;
                //valueBox.addItem("Unknown");
                for (int i = 0; i < attribute.numValues(); i++) {
                    valueBox.addItem(attribute.value(i));
                }
                valueBox.setSelectedIndex(1); // Select false first
                valueBox.addActionListener(l);
            }
            else if (attribute.isNumeric()) {
                JTextField tf = new JTextField();
                valueBoxes[index] = tf;
                tf.addActionListener(l);
                
            }
            constraints.gridx = 1;
            predicatorPane.add(valueBoxes[index], constraints);
            index ++;
            y ++;
        }
        frame.getContentPane().add(predicatorPane, BorderLayout.CENTER);
        // Class Pane
        JPanel clsPane = new JPanel();
        clsPane.setBorder(BorderFactory.createEtchedBorder());
        clsPane.setLayout(new GridBagLayout());
        constraints.gridx = 0;
        constraints.gridy = 0;
        constraints.gridheight = 2;
        attribute = dataset.classAttribute();
        JLabel clsLabel = new JLabel(attribute.name());
        clsPane.add(clsLabel, constraints);
        constraints.gridx = 1;
        constraints.gridheight = 1;
        for (int i = 0; i < attribute.numValues(); i++) {
            JLabel label = new JLabel();
            label.setText(attribute.value(i));
            constraints.gridy = i;
            clsPane.add(label, constraints);
        }
        constraints.gridx = 2;
        constraints.gridy = 0;
        clsPane.add(firstValueLabel, constraints);
        constraints.gridy = 1;
        clsPane.add(secondValueLabel, constraints);
        frame.getContentPane().add(clsPane, BorderLayout.SOUTH);
    }
    
}
