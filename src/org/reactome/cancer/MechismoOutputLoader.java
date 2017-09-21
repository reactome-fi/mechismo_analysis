package org.reactome.cancer;

import org.apache.log4j.Logger;
import org.reactome.cancer.driver.Interactome3dDriverAnalyzer;
import org.reactome.r3.util.FileUtility;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class MechismoOutputLoader {
    private final static Logger logger = Logger.getLogger(MechismoOutputLoader.class);
    private final int primaryIdA1Idx = 1;
    private final int primaryIdB1Idx = 19;
    private final int rxnKlassIdx = 21;
    private final String rxnKlass = "hetero";
    private String mechismoOuputFilePath = null;

    public MechismoOutputLoader(){

    }

    public MechismoOutputLoader(String mechismoOutputFilePath){
        this.mechismoOuputFilePath = mechismoOutputFilePath;
    }

    public Set<String> extractMechismoFIs(){
        return this.extractMechismoFIs(this.mechismoOuputFilePath);
    }

    public Set<String> extractMechismoFIs(String mechismoOuputFilePath)  {
        Set<String> fis = new HashSet<>();
        FileUtility fileUtility = new FileUtility();
        try {
            fileUtility.setInput(mechismoOuputFilePath);
            String line = null;
            String[] tokens = null;
            while ((line = fileUtility.readLine()) != null) {
                tokens = line.split("\t");
                if(tokens.length > rxnKlassIdx &&
                        rxnKlass.equals(tokens[rxnKlassIdx])) {
                    fis.add(String.format("%s\t%s",
                            tokens[primaryIdA1Idx],
                            tokens[primaryIdB1Idx]));
                    fis.add(String.format("%s\t%s",
                            tokens[primaryIdB1Idx],
                            tokens[primaryIdA1Idx]));
                }
            }
            fileUtility.close();
        }catch(IOException ioe){
            logger.error(String.format("Couldn't use %s, %s: %s",
                    mechismoOuputFilePath.toString(),
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
        return fis;
    }
}
