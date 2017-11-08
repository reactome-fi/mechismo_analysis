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
    private final String patientPatternString = "([A-Z]+:TCGA-[0-9A-Z-]+)";
    private double mechScoreLowerBoundInclusive;
    private double eThresh;
    private String mechismoOuputFilePath = null;
    private String mechismoFIFilterFilePath = null;
    private Set<FI> fis = null;
    private Set<String> fiFilter = null;
    private Map<Mutation, Integer> mut2SampleCount = null;
    private Map<Gene, Integer> gene2SampleCount = null;
    private Map<Patient, Set<FI>> samples2fis = null;
    private Map<FI, Set<Patient>> fis2Samples = null;
    //TODO: decide whether or not to add protein B pos/res info to muts?
    private Map<Patient, Map<FI, Set<Mutation>>> samples2fis2muts = null;
    private Pattern sampleIDPattern = null;
    private Matcher sampleIDMatcher = null;

    public MechismoOutputLoader(String mechismoOutputFilePath,
                                String mechismoFIFilterFilePath,
                                double mechScoreLowerBoundInclusive,
                                double eThresh) {
        this.mechismoOuputFilePath = mechismoOutputFilePath;
        this.mechismoFIFilterFilePath = mechismoFIFilterFilePath;
        this.mechScoreLowerBoundInclusive = mechScoreLowerBoundInclusive;
        this.eThresh = eThresh;
        this.sampleIDPattern = Pattern.compile(this.patientPatternString);
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
        this.samples2fis = new HashMap<>();
        this.fis2Samples = new HashMap<>();
        this.mut2SampleCount = new HashMap<>();
        this.gene2SampleCount = new HashMap<>();
        this.samples2fis2muts = new HashMap<>();
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
                    String hgncNameA = tokens[nameA1Idx];
                    String hgncNameB = tokens[nameB1Idx];
                    if (this.fiFilter == null ||
                            this.fiFilter.contains(String.format("%s\t%s",
                            hgncNameA,
                            hgncNameB))) {
                        Double mechScore = new Double(tokens[mechScoreIdx]);
                        if (Math.abs(mechScore) >= mechScoreLowerBoundInclusive) {
                            String uniprotIDA = tokens[primaryIdA1Idx];
                            String uniprotIDB = tokens[primaryIdB1Idx];
                            Integer position = Integer.parseInt(tokens[posA1Idx]);
                            Character normalResidue = new Character(tokens[resA1Idx].charAt(0));
                            Character mutationResidue = new Character(tokens[mutA1Idx].charAt(0));
                            Gene gene1 = new Gene(hgncNameA,uniprotIDA);
                            Gene gene2 = new Gene(hgncNameB,uniprotIDB);
                            FI fi = new FI(gene1,gene2);
                            String userInput = tokens[userInputIdx];
                            //mutated gene is always the first one listed (partner A)
                            Mutation mutation = new Mutation(
                                    gene1,
                                    position,
                                    normalResidue,
                                    mutationResidue);
                            this.sampleIDMatcher = this.sampleIDPattern.matcher(userInput);
                            while (this.sampleIDMatcher.find()) {
                                String[] patientData = this.sampleIDMatcher.group().split(":");

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

                                //update samples2fis
                                Patient patient = new Patient(patientData[1],patientData[0]);
                                Set<FI> sampleFIs;
                                if (samples2fis.containsKey(patient)) {
                                    sampleFIs = samples2fis.get(patient);
                                } else {
                                    sampleFIs = new HashSet<>();
                                }
                                sampleFIs.add(fi);
                                samples2fis.put(patient, sampleFIs);

                                //update fis2samples
                                Set<Patient> fiSamples;
                                if (fis2Samples.containsKey(fi)){
                                    fiSamples = fis2Samples.get(fi);
                                } else {
                                    fiSamples = new HashSet<>();
                                }
                                fiSamples.add(patient);
                                fis2Samples.put(fi, fiSamples);

                                //update samples2fis2muts
                                Map<FI, Set<Mutation>> fis2muts;
                                if (samples2fis2muts.containsKey(patient)) {
                                    fis2muts = samples2fis2muts.get(patient);
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
                                samples2fis2muts.put(patient, fis2muts);
                            }
                            fis.add(fi);
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

    public Map<Patient, Map<FI, Set<Mutation>>> ExtractSamples2FIs2Muts() {
        return this.ExtractSamples2FIs2Muts(this.mechismoOuputFilePath);
    }

    public Map<Patient, Map<FI, Set<Mutation>>> ExtractSamples2FIs2Muts(String mechismoOuputFilePath) {
        if (this.samples2fis2muts == null) {
            this.ParseMechismoOutputFile(mechismoOuputFilePath);
        }
        return this.samples2fis2muts;
    }

    public Map<FI, Set<Patient>> ExtractFIs2Samples() {
        return this.ExtractFIs2Samples(this.mechismoOuputFilePath);
    }

    public Map<FI, Set<Patient>> ExtractFIs2Samples(String mechismoOuputFilePath) {
        if (this.fis2Samples == null) {
            this.ParseMechismoOutputFile(mechismoOuputFilePath);
        }
        return this.fis2Samples;
    }

    public Map<Patient, Set<FI>> ExtractSamples2FIs() {
        return this.ExtractSamples2FIs(this.mechismoOuputFilePath);
    }

    public Map<Patient, Set<FI>> ExtractSamples2FIs(String mechismoOuputFilePath) {
        if (this.samples2fis == null) {
            this.ParseMechismoOutputFile(mechismoOuputFilePath);
        }
        return this.samples2fis;
    }

    public Set<FI> ExtractMechismoFIs() {
        return this.ExtractMechismoFIs(this.mechismoOuputFilePath);
    }

    public Set<FI> ExtractMechismoFIs(String mechismoOuputFilePath) {
        if (this.fis == null) {
            this.ParseMechismoOutputFile(mechismoOuputFilePath);
        }
        return this.fis;
    }

    public Map<Mutation, Integer> getMut2SampleCount() {
        return mut2SampleCount;
    }

    public Map<Gene, Integer> getGene2SampleCount() {
        return gene2SampleCount;
    }
}
