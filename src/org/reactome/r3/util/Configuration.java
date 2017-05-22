/*
 * Created on May 17, 2017
 *
 */
package org.reactome.r3.util;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

import org.apache.log4j.Logger;
import org.gk.persistence.MySQLAdaptor;

/**
 * Some generation configurations externalized into a property file in the resources directory.
 * Contents in that file should NOT be placed into git.
 * @author gwu
 *
 */
public class Configuration {
    private static final Logger logger = Logger.getLogger(Configuration.class);
    private static final String CONFIG_FILE_NAME = "resources/config.prop";
    private Properties properties;
    private static Configuration configuration;
    
    /**
     * Default constructor.
     */
    private Configuration() {
        load();
    }
    
    public static Configuration getConfiguration() {
        if (configuration == null)
            configuration = new Configuration();
        return configuration;
    }
    
    public MySQLAdaptor getReactomeDBA() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor(properties.getProperty("reactome.dbHost"), 
                                            properties.getProperty("reactome.dbName"), 
                                            properties.getProperty("reactome.dbUser"), 
                                            properties.getProperty("reactome.dbPwd"));
        return dba;
    }
    
    private void load() {
        properties = new Properties();
        try {
            FileInputStream fis = new FileInputStream(CONFIG_FILE_NAME);
            properties.load(fis);
        }
        catch(IOException e) {
            logger.error(e.getMessage(), e);
        }
    }
}
