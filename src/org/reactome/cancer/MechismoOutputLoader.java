package org.reactome.cancer;

import org.apache.log4j.Logger;
import org.reactome.r3.util.FileUtility;

import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MechismoOutputLoader {
    private final static Logger logger = Logger.getLogger(MechismoOutputLoader.class);
    private final int primaryIdA1Idx = 1;
    private final int userInputIdx = 6;
    private final int primaryIdB1Idx = 19;
    private final int rxnKlassIdx = 21;
    private final String rxnKlass = "hetero";
    private final String sampleIDPatternString = "(TCGA-[0-9A-Z-]+)";
    private String mechismoOuputFilePath = null;
    private Set<String> fis = null;
    private Map<String,Set<String>> samples2fis = null;
    private Map<String,Set<String>> fis2Samples = null;
    private Pattern sampleIDPattern = null;
    private Matcher sampleIDMatcher = null;

    public MechismoOutputLoader(){
        this(null);
    }

    public MechismoOutputLoader(String mechismoOutputFilePath){
        this.mechismoOuputFilePath = mechismoOutputFilePath;
        this.sampleIDPattern = Pattern.compile(this.sampleIDPatternString);
    }

    private void ParseMechismoOutputFile(){
        this.ParseMechismoOutputFile(this.mechismoOuputFilePath);
    }

    private void ParseMechismoOutputFile(String mechismoOuputFilePath){
        this.fis = new HashSet<>();
        this.samples2fis = new HashMap<>();
        this.fis2Samples = new HashMap<>();
        FileUtility fileUtility = new FileUtility();
        try {
            fileUtility.setInput(mechismoOuputFilePath);
            String line = null;
            String[] tokens = null;
            while ((line = fileUtility.readLine()) != null) {
                tokens = line.split("\t");
                if(tokens.length > rxnKlassIdx &&
                        rxnKlass.equals(tokens[rxnKlassIdx])) {
                    String forwardFI = String.format("%s\t%s",
                            tokens[primaryIdA1Idx],
                            tokens[primaryIdB1Idx]);
                    String backwardFI = String.format("%s\t%s",
                            tokens[primaryIdB1Idx],
                            tokens[primaryIdA1Idx]);
                    String userInput = tokens[userInputIdx];
                    this.sampleIDMatcher = this.sampleIDPattern.matcher(userInput);
                    while(this.sampleIDMatcher.find()){
                        String sampleID = this.sampleIDMatcher.group();

                        //update samples2fis
                        Set<String> sampleFIs;
                        if(samples2fis.containsKey(sampleID)){
                            sampleFIs = samples2fis.get(sampleID);
                        }else {
                            sampleFIs = new HashSet<>();
                        }
                        sampleFIs.add(forwardFI);
                        sampleFIs.add(backwardFI);
                        samples2fis.put(sampleID,sampleFIs);

                        //update fis2samples
                        Set<String> fiSamples;
                        if(fis2Samples.containsKey(forwardFI) ||
                                fis2Samples.containsKey(backwardFI)){
                            if(!fis2Samples.get(forwardFI)
                                    .equals(fis2Samples.get(backwardFI))){
                                throw new IllegalStateException("These should always match");
                            }
                            fiSamples = fis2Samples.get(forwardFI);
                        }else{
                            fiSamples = new HashSet<>();
                        }
                        fiSamples.add(sampleID);
                        fis2Samples.put(forwardFI,fiSamples);
                        fis2Samples.put(backwardFI,fiSamples);
                    }
                    fis.add(forwardFI);
                    fis.add(backwardFI);
                }
            }
            fileUtility.close();
        }catch(IOException ioe){
            logger.error(String.format("Couldn't use %s, %s: %s",
                    mechismoOuputFilePath.toString(),
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    public Map<String,Set<String>> ExtractFIs2Samples(){
        return this.ExtractFIs2Samples(this.mechismoOuputFilePath);
    }

    public Map<String,Set<String>> ExtractFIs2Samples(String mechismoOuputFilePath){
        if(this.fis2Samples == null){
            this.ParseMechismoOutputFile(mechismoOuputFilePath);
        }
        return this.fis2Samples;
    }

    public Map<String,Set<String>> ExtractSamples2FIs(){
        return this.ExtractSamples2FIs(this.mechismoOuputFilePath);
    }

    public Map<String,Set<String>> ExtractSamples2FIs(String mechismoOuputFilePath){
        if(this.samples2fis == null){
            this.ParseMechismoOutputFile(mechismoOuputFilePath);
        }
        return this.samples2fis;
    }

    public Set<String> ExtractMechismoFIs(){
        return this.ExtractMechismoFIs(this.mechismoOuputFilePath);
    }

    public Set<String> ExtractMechismoFIs(String mechismoOuputFilePath)  {
        if(this.fis == null){
            this.ParseMechismoOutputFile(mechismoOuputFilePath);
        }
        return this.fis;
    }
}
