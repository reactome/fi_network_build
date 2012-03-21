/*
 * Created on Apr 7, 2009
 *
 */
package org.reactome.fi;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import org.reactome.fi.util.Value;
import org.reactome.weka.NaiveBayesClassifier;

/**
 * This GUI test is used to check the results from NaiveBayesClassifier.
 * @author wgm
 *
 */
public class NBCGUITest extends JFrame {
    private NaiveBayesClassifier classifier;
    // A list of JComboBox
    private List<JComboBox> listBoxes;
    private JTextField scoreTF;
    
    public NBCGUITest() {
        init();
    }
    
    private void init() {
        loadNBC();
        initGUI();
    }
    
    private void loadNBC() {
        try {
            NBCAnalyzer analyzer = new NBCAnalyzer();
            classifier = analyzer.loadSavedNBC();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }
    
    private void initGUI() {
        List<String> featureList = classifier.getFeatureList();
        // Create JPanel
        JPanel contentPane = new JPanel();
        contentPane.setLayout(new GridBagLayout());
        GridBagConstraints constraints = new GridBagConstraints();
        constraints.insets = new Insets(2, 4, 2, 4);
        int gridY = 0;
        ActionListener l = createActionListener();
        listBoxes = new ArrayList<JComboBox>();
        for (String feature : featureList) {
            JLabel label = new JLabel(feature);
            Boolean[] values = new Boolean[] {Boolean.FALSE, Boolean.TRUE};
            JComboBox list = new JComboBox(values);
            list.addActionListener(l);
            list.setEditable(false);
            listBoxes.add(list);
            constraints.gridx = 0;
            constraints.gridy = gridY;
            contentPane.add(label, constraints);
            constraints.gridx = 1;
            contentPane.add(list, constraints);
            gridY ++;
        }
        // Add a result
        JLabel scoreLabel = new JLabel("Score: ");
        scoreTF = new JTextField();
        scoreTF.setEditable(false);
        scoreTF.setColumns(15);
        constraints.insets = new Insets(10, 4, 2, 4);
        constraints.gridy = gridY;
        constraints.gridx = 0;
        contentPane.add(scoreLabel, constraints);
        constraints.gridx = 1;
        contentPane.add(scoreTF, constraints);
        getContentPane().add(contentPane, BorderLayout.CENTER);
    }
    
    private ActionListener createActionListener() {
        ActionListener l = new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                updateScore();
            }
        };
        return l;
    }
    
    private void updateScore() {
        Map<String, Field> featureToField = Value.convertValueFeatureToField(classifier.getFeatureList());
        Value value = new Value();
        try {
            for (int i = 0; i < listBoxes.size(); i++) {
                JComboBox box = listBoxes.get(i);
                Boolean featureValue = (Boolean) box.getSelectedItem();
                String feature = classifier.getFeatureList().get(i);
                Field field = featureToField.get(feature);
                field.set(value, featureValue);
            }
        }
        catch(IllegalAccessException e) {
            e.printStackTrace();
        }
        double score = classifier.calculateScore(value);
        scoreTF.setText(score + "");
    }
    
    public static void main(String[] args) {
        NBCGUITest test = new NBCGUITest();
        test.setLocation(465, 80);
        test.setSize(425, 550);
        test.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        test.setVisible(true);
    }
    
}
