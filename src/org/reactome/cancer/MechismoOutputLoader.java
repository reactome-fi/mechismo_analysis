package org.reactome.cancer;

import org.apache.log4j.Logger;
import org.reactome.r3.util.FileUtility;

import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MechismoOutputLoader {
    private final static Logger logger = Logger.getLogger(MechismoOutputLoader.class);
    private final double mechScoreThresh = 2.0;
    private final double eThresh = 1e-10;
    private final int primaryIdA1Idx = 1;
    private final int posA1Idx = 3;
    private final int resA1Idx = 4;
    private final int mutA1Idx = 5;
    private final int userInputIdx = 6;
    private final int mismatchIdx = 8;
    private final int mechScoreIdx = 17;
    private final int primaryIdB1Idx = 19;
    private final int rxnKlassIdx = 21;
    private final int confIdx = 26;
    private final int eaIdx = 35;
    private final int posB1Idx = 37;
    private final int resB1Idx = 38;
    private final int ebIdx = 43;
    private final String rxnKlass = "hetero";
    private final String sampleIDPatternString = "(TCGA-[0-9A-Z-]+)";
    private String mechismoOuputFilePath = null;
    private Set<String> fis = null;
    private Map<String,Set<String>> samples2fis = null;
    private Map<String,Set<String>> fis2Samples = null;
    //TODO: decide whether or not to add protein B pos/res info to muts?
    private Map<String,Map<String,Set<String>>> samples2fis2muts = null;
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
        this.samples2fis2muts = new HashMap<>();
        FileUtility fileUtility = new FileUtility();
        try {
            fileUtility.setInput(mechismoOuputFilePath);
            String line = null;
            String[] tokens = null;
            while ((line = fileUtility.readLine()) != null) {
                tokens = line.split("\t");
                if(tokens.length > ebIdx &&
                        rxnKlass.equals(tokens[rxnKlassIdx]) &&
                        new Double(tokens[eaIdx]) < eThresh &&
                        new Double(tokens[ebIdx]) < eThresh){
                    String protA = tokens[primaryIdA1Idx];
                    String protB = tokens[primaryIdB1Idx];
                    Double mechScore = Math.abs(new Double(tokens[mechScoreIdx]));
                    if(mechScore > mechScoreThresh) {
                        String forwardFI = String.format("%s\t%s",
                                protA,
                                protB);
                        String backwardFI = String.format("%s\t%s",
                                protB,
                                protA);
                        String userInput = tokens[userInputIdx];
                        String mutString = String.format("%s:pos%s-res%s-mut%s",
                                tokens[primaryIdA1Idx],
                                tokens[posA1Idx],
                                tokens[resA1Idx],
                                tokens[mutA1Idx]);
                        this.sampleIDMatcher = this.sampleIDPattern.matcher(userInput);
                        while (this.sampleIDMatcher.find()) {
                            String sampleID = this.sampleIDMatcher.group();

                            //update samples2fis
                            Set<String> sampleFIs;
                            if (samples2fis.containsKey(sampleID)) {
                                sampleFIs = samples2fis.get(sampleID);
                            } else {
                                sampleFIs = new HashSet<>();
                            }
                            sampleFIs.add(forwardFI);
                            sampleFIs.add(backwardFI);
                            samples2fis.put(sampleID, sampleFIs);

                            //update fis2samples
                            Set<String> fiSamples;
                            if (fis2Samples.containsKey(forwardFI) ||
                                    fis2Samples.containsKey(backwardFI)) {
                                if (!fis2Samples.get(forwardFI)
                                        .equals(fis2Samples.get(backwardFI))) {
                                    throw new IllegalStateException("These should always match");
                                }
                                fiSamples = fis2Samples.get(forwardFI);
                            } else {
                                fiSamples = new HashSet<>();
                            }
                            fiSamples.add(sampleID);
                            fis2Samples.put(forwardFI, fiSamples);
                            fis2Samples.put(backwardFI, fiSamples);

                            //update samples2fis2muts
                            Map<String,Set<String>> fis2muts;
                            if(samples2fis2muts.containsKey(sampleID)){
                                fis2muts = samples2fis2muts.get(sampleID);
                            }else{
                                fis2muts = new HashMap<>();
                            }
                            Set<String> muts;
                            if(fis2muts.containsKey(forwardFI)){
                                muts = fis2muts.get(forwardFI);
                            }else{
                                muts = new HashSet<>();
                            }
                            muts.add(mutString);
                            fis2muts.put(forwardFI,muts);
                            fis2muts.put(backwardFI,muts);
                            samples2fis2muts.put(sampleID,fis2muts);
                        }
                        fis.add(forwardFI);
                        fis.add(backwardFI);
                    }
                    else{
                        //System.out.println(String.format("Ignoring FI %s-%s with mechScore = %.3f",
                        //        protA,protB,mechScore));
                    }
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

    public Map<String,Map<String,Set<String>>> ExtractSamples2FIs2Muts(){
        return this.ExtractSamples2FIs2Muts(this.mechismoOuputFilePath);
    }

    public Map<String,Map<String,Set<String>>> ExtractSamples2FIs2Muts(String mechismoOuputFilePath){
        if(this.samples2fis2muts == null){
            this.ParseMechismoOutputFile(mechismoOuputFilePath);
        }
        return this.samples2fis2muts;
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
