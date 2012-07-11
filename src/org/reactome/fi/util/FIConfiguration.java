/*
 * Created on Mar 28, 2007
 *
 */
package org.reactome.fi.util;

import java.io.File;
import java.io.FileInputStream;
import java.sql.Connection;
import java.sql.DriverManager;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**
 * This class holds configuration for building a new version of FI network.
 * @author guanming
 *
 */
public class FIConfiguration {
    private static FIConfiguration configuration;
    private Properties properties;
    
    private FIConfiguration() {
        init();
    }
    
    public static Connection getConnection(String dbName) throws Exception {
        Class.forName("com.mysql.jdbc.Driver");
        String url = "jdbc:mysql://localhost:3306/" + dbName;
        Properties info = new Properties();
        info.setProperty("user", getConfiguration().get("DB_USER"));
        info.setProperty("password", getConfiguration().get("DB_PWD"));
        Connection connection = DriverManager.getConnection(url, info);
        return connection;
    }
    
    public static FIConfiguration getConfiguration() {
        if (configuration == null)
            configuration = new FIConfiguration();
        return configuration;
    }
    
    public String get(String name) {
        return properties.getProperty(name);
    }
    
    public Properties getAll() {
        return this.properties;
    }
    
    private void init() {
        // Load configurations from a configuration file
        properties = new Properties();
        try {
            File file = new File("resources/configuration.prop");
            properties.load(new FileInputStream(file));
            normalizeProperties();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Make sure all replaceable properties have been replaced. For example,
     * ${RESULT_DIR} should be replaced by actual values as following
     * 
     */
    private void normalizeProperties() throws Exception {
        String regExp = "\\$\\{(\\w+)\\}";
        Pattern pattern = Pattern.compile(regExp);
        boolean hasSymbols = true;
        while (hasSymbols) {
            hasSymbols = false;
            for (Object name : properties.keySet()) {
                String value = properties.getProperty(name.toString());
                Matcher matcher = pattern.matcher(value);
//                System.out.println(name + "=" + value);
                if (matcher.find()) {
                    hasSymbols = true;
                    String variable = matcher.group(1);
                    String variableValue = properties.getProperty(variable);
                    matcher = pattern.matcher(variableValue);
                    // Don't replace right now
                    if (matcher.find())
                        continue;
                    value = value.replaceAll(regExp, variableValue);
                    properties.setProperty(name.toString(), value);
                }
            }
        }
    }
    
//    @Test
//    public void testLoading() {
//        FIConfigruation configuration = FIConfigruation.getConfigruation();
//        Properties prop = configuration.getProperties();
//        for (Object key : prop.keySet()) {
//            Object value = prop.get(key);
//            System.out.println(key + "=" + value);
//        }
//    }
}
