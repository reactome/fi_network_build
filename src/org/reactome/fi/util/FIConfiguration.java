/*
 * Created on Mar 28, 2007
 *
 */
package org.reactome.fi.util;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
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
        try {
            init(new FileInputStream("resources/configuration.prop")); // Default
        }
        catch(IOException e) {
            e.printStackTrace();
        }
    }
    
    private FIConfiguration(InputStream is) {
        init(is);
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
    
    /**
     * Provide a configuration based on the passed configuration file. The original configuration
     * object will be overwritten. 
     * @param configFileName
     * @return
     */
    public static FIConfiguration getConfiguration(InputStream is) {
        if (configuration == null)
            configuration = new FIConfiguration(is);
        return configuration;
    }
    
    public String get(String name) {
        return properties.getProperty(name);
    }
    
    public Properties getAll() {
        return this.properties;
    }
    
    private void init(InputStream is) {
        // Load configurations from a configuration file
        properties = new Properties();
        try {
            properties.load(is);
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
                    value = value.replace("${" + variable + "}", variableValue);
                    properties.setProperty(name.toString(), value);
                }
            }
        }
    }
    
//    @Test
//    public void testLoading() {
//        FIConfiguration configuration = new FIConfiguration();
//        for (Object key : properties.keySet()) {
//            Object value = properties.get(key);
//            System.out.println(key + "=" + value);
//        }
//    }
}
