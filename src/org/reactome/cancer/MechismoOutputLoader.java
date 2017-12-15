package org.reactome.cancer;

import org.apache.log4j.Logger;
import org.reactome.r3.util.FileUtility;

import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MechismoOutputLoader {
    private final static Logger logger = Logger.getLogger(MechismoOutputLoader.class);
    private final int nameA1Idx = 0;
    private final int primaryIdA1Idx = 1;
    private final int posA1Idx = 3;
    private final int resA1Idx = 4;
    private final int mutA1Idx = 5;
    private final int userInputIdx = 6;
    private final int mismatchIdx = 8;
    private final int mechScoreIdx = 17;
    private final int nameB1Idx = 18;
    private final int primaryIdB1Idx = 19;
    private final int rxnKlassIdx = 21;
    private final int confIdx = 26;
    private final int eaIdx = 35;
    private final int posB1Idx = 37;
    private final int resB1Idx = 38;
    private final int ebIdx = 43;
    private final String rxnKlass = "hetero";
    private final String junkString = "{ECO:0000312";
    private final String patientPatternString = "([A-Z]+:TCGA-[0-9A-Z-]+)";
    private double mechScoreLowerBoundInclusive;
    private double eThresh;
    private String tcgaCancerType;
    private String mechismoOuputFilePath;
    private String mechismoFIFilterFilePath;

    private Set<FI> fis;
    private Map<Patient, Set<FI>> patientsToFIs;
    private Map<FI, Set<Patient>> fisToPatients;
    private Map<Patient, Map<FI, Set<Mutation>>> patientToFIsToMutations;

    private Map<Gene, Integer> gene2SampleCount;
    private Map<Mutation, Integer> mut2SampleCount;
    private Set<String> fiFilter;
    private Pattern sampleIDPattern;
    private Matcher sampleIDMatcher;

    public MechismoOutputLoader(String mechismoOuputFilePath){
        this(mechismoOuputFilePath,null,null,0.0d,1.0d);
    }

    public MechismoOutputLoader(String mechismoOutputFilePath,
                                String mechismoFIFilterFilePath,
                                String tcgaCancerType,
                                double mechScoreLowerBoundInclusive,
                                double eThresh) {
        this.mechismoOuputFilePath = mechismoOutputFilePath;
        this.mechismoFIFilterFilePath = mechismoFIFilterFilePath;
        this.tcgaCancerType = tcgaCancerType;
        this.mechScoreLowerBoundInclusive = mechScoreLowerBoundInclusive;
        this.eThresh = eThresh;
        this.sampleIDPattern = Pattern.compile(this.patientPatternString);
        ParseMechismoOutputFile(mechismoOuputFilePath);
    }

    public Set<FI> getFIs(Patient patient){
        return this.patientToFIsToMutations.get(patient).keySet();
    }

    public Set<FI> getFIs(){
        return this.fis;
    }

    public Set<Patient> getPatients(FI fi){
        return this.fisToPatients.get(fi);
    }

    public Set<Patient> getPatients(){
       return this.patientsToFIs.keySet();
    }

    public Boolean hasMutations(Patient patient){
        return this.patientToFIsToMutations.containsKey(patient);
    }

    public Set<Mutation> getMutations(Patient patient, FI fi){
        return this.patientToFIsToMutations.get(patient).get(fi);
    }

    private void ParseMechismoFIFilterFile() {
        if (this.mechismoFIFilterFilePath != null) {
            this.fiFilter = new HashSet<>();
            FileUtility fileUtility = new FileUtility();
            try {
                fileUtility.setInput(mechismoFIFilterFilePath);
                String line;
                String[] tokens;
                fileUtility.readLine(); //skip header line
                while ((line = fileUtility.readLine()) != null) {
                    tokens = line.split("\t");
                    if (Integer.parseInt(tokens[2]) > Integer.parseInt(tokens[3]) &&
                            Double.parseDouble(tokens[10]) < 0.01) {
                        this.fiFilter.add(String.format("%s\t%s",
                                tokens[0],
                                tokens[1]));
                        this.fiFilter.add(String.format("%s\t%s",
                                tokens[1],
                                tokens[0]));
                    }
                }
                fileUtility.close();
            } catch (IOException ioe) {
                logger.error(String.format("Couldn't use %s, %s: %s",
                        mechismoOuputFilePath,
                        ioe.getMessage(),
                        Arrays.toString(ioe.getStackTrace())));
            }
        }
    }

    private void ParseMechismoOutputFile(String mechismoOuputFilePath) {
        ParseMechismoFIFilterFile();
        this.fis = new HashSet<>();
        this.patientsToFIs = new HashMap<>();
        this.fisToPatients = new HashMap<>();
        this.mut2SampleCount = new HashMap<>();
        this.gene2SampleCount = new HashMap<>();
        this.patientToFIsToMutations = new HashMap<>();
        FileUtility fileUtility = new FileUtility();
        try {
            fileUtility.setInput(mechismoOuputFilePath);
            String line;
            String[] tokens;
            while ((line = fileUtility.readLine()) != null) {
                tokens = line.split("\t");
                if (tokens.length > ebIdx &&
                        rxnKlass.equals(tokens[rxnKlassIdx]) &&
                        new Double(tokens[eaIdx]) < eThresh &&
                        new Double(tokens[ebIdx]) < eThresh) {
                    //remove junk strings
                    for(int i = 0; i < tokens.length; i++){
                        if(tokens[i].contains(junkString)){
                            tokens[i] = tokens[i].split(" ")[0];
                        }
                    }
                    String hgncNameA = tokens[nameA1Idx];
                    String hgncNameB = tokens[nameB1Idx];
                    if (this.fiFilter == null ||
                            this.fiFilter.contains(String.format("%s\t%s",
                                    hgncNameA,
                                    hgncNameB))) {
                        String userInput = tokens[userInputIdx];
                        if (!userInput.toUpperCase().contains("SILENT")) {
                            Double mechScore = new Double(tokens[mechScoreIdx]);
                            if (Math.abs(mechScore) >= mechScoreLowerBoundInclusive) {
                                String uniprotIDA = tokens[primaryIdA1Idx];
                                String uniprotIDB = tokens[primaryIdB1Idx];
                                Integer position = Integer.parseInt(tokens[posA1Idx]);
                                Character normalResidue = tokens[resA1Idx].charAt(0);
                                Character mutationResidue = tokens[mutA1Idx].charAt(0);
                                Gene gene1 = new Gene(hgncNameA, uniprotIDA);
                                Gene gene2 = new Gene(hgncNameB, uniprotIDB);
                                FI fi = new FI(gene1, gene2);

                                //mutated gene is always the first one listed (partner A)
                                Mutation mutation = new Mutation(
                                        gene1,
                                        position,
                                        normalResidue,
                                        mutationResidue);
                                this.sampleIDMatcher = this.sampleIDPattern.matcher(userInput);
                                while (this.sampleIDMatcher.find()) {
                                    String[] patientData = this.sampleIDMatcher.group().split(":");
                                    if (this.tcgaCancerType == null ||
                                            patientData[0].toUpperCase().contains(this.tcgaCancerType)) {
                                        Integer geneCount = 1;
                                        if (this.gene2SampleCount.containsKey(gene1)) {
                                            geneCount += this.gene2SampleCount.get(gene1);
                                        }
                                        this.gene2SampleCount.put(gene1, geneCount);

                                        Integer mutationCount = 1;
                                        if (this.mut2SampleCount.containsKey(mutation)) {
                                            mutationCount += this.mut2SampleCount.get(mutation);
                                        }
                                        this.mut2SampleCount.put(mutation, mutationCount);

                                        //update patientsToFIs
                                        Patient patient = new Patient(
                                                patientData[1].trim().toUpperCase(),
                                                patientData[0].trim().toUpperCase());

                                        Set<FI> sampleFIs;
                                        if (patientsToFIs.containsKey(patient)) {
                                            sampleFIs = patientsToFIs.get(patient);
                                        } else {
                                            sampleFIs = new HashSet<>();
                                        }
                                        sampleFIs.add(fi);
                                        patientsToFIs.put(patient, sampleFIs);

                                        //update fis2samples
                                        Set<Patient> fiSamples;
                                        if (fisToPatients.containsKey(fi)) {
                                            fiSamples = fisToPatients.get(fi);
                                        } else {
                                            fiSamples = new HashSet<>();
                                        }
                                        fiSamples.add(patient);
                                        fisToPatients.put(fi, fiSamples);

                                        //update patientToFIsToMutations
                                        Map<FI, Set<Mutation>> fis2muts;
                                        if (patientToFIsToMutations.containsKey(patient)) {
                                            fis2muts = patientToFIsToMutations.get(patient);
                                        } else {
                                            fis2muts = new HashMap<>();
                                        }
                                        Set<Mutation> muts;
                                        if (fis2muts.containsKey(fi)) {
                                            muts = fis2muts.get(fi);
                                        } else {
                                            muts = new HashSet<>();
                                        }
                                        muts.add(mutation);
                                        fis2muts.put(fi, muts);
                                        patientToFIsToMutations.put(patient, fis2muts);
                                        fis.add(fi);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            fileUtility.close();
        } catch (IOException ioe) {
            logger.error(String.format("Couldn't use %s, %s: %s",
                    mechismoOuputFilePath,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }
}
