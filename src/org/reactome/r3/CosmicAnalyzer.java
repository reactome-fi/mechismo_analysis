/*
 * Created on Sep 6, 2016
 *
 */
package org.reactome.r3;

import org.junit.Test;
import org.reactome.cancer.MutationObservation;
import org.reactome.r3.util.FileUtility;

import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * This class is used to perform analysis on the cosmic data.
 *
 * @author gwu
 */
public class CosmicAnalyzer {
    //    private final String DIR_NAME = "datasets/COSMIC/v78/";
    // As of February 27, 2017, use v80
    private final String DIR_NAME = "datasets/COSMIC/v80/";

    /**
     * Default constructor.
     */
    public CosmicAnalyzer() {

    }

    @Test
    public void checkMutationFile() throws IOException {
        String fileName = DIR_NAME + "CosmicMutantExport.tsv";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        System.out.println(line);
        fu.close();
    }

    /**
     * Filter entries to a list of coordinates in the peptide.
     *
     * @param entries
     * @param coordinates
     * @return
     */
    public List<MutationObservation> filterEntries(List<MutationObservation> entries,
                                                   Collection<Integer> coordinates) {
        List<MutationObservation> rtn = new ArrayList<MutationObservation>();
        for (MutationObservation entry : entries) {
            if (coordinates.contains(entry.getCoordinate()))
                rtn.add(entry);
//            else if (coordinates.contains(entry.coordinate + 1))
//                rtn.add(entry);
//            else if (coordinates.contains(entry.coordinate - 1))
//                rtn.add(entry);
        }
        return rtn;
    }

    /**
     * Load mutations with missense and nonsense, pathogenic only impact, and needZygosity.
     *
     * @param genes
     * @return
     * @throws IOException
     */
    public Map<String, List<MutationObservation>> loadMutations(Set<String> genes) throws IOException {
        Set<String> mutationTypes = new HashSet<String>();
        mutationTypes.add("Substitution - Missense");
//        mutationTypes.add("Substitution - Nonsense");
        mutationTypes.add("Insertion - In frame");
        mutationTypes.add("Deletion - In frame");
        return loadMutations(genes, true, true, mutationTypes);
    }

    @Test
    public void checkMutationTypes() throws IOException {
        String fileName = DIR_NAME + "CosmicMutantExport.tsv";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        // Header line
        String line = fu.readLine();
        String[] headers = line.split("\t");
        for (int i = 0; i < headers.length; i++)
            System.out.println(i + "\t" + headers[i]);
        Set<String> mutationTypes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            mutationTypes.add(tokens[19]); // header 19 should be "Mutation Description"
        }
        fu.close();
        System.out.println("Total mutation types: " + mutationTypes.size());
        List<String> list = new ArrayList<String>(mutationTypes);
        Collections.sort(list);
        for (String type : list)
            System.out.println(type);
    }

    /**
     * Load entries in COSMIC for a set of genes.
     *
     * @param genes
     * @param pathogenicOnly
     * @param needZygosity
     * @param mutationTypes
     * @return
     * @throws IOException
     */
    public Map<String, List<MutationObservation>> loadMutations(Set<String> genes,
                                                                boolean pathogenicOnly,
                                                                boolean needZygosity,
                                                                Set<String> mutationTypes) throws IOException {
        List<MutationObservation> list = new ArrayList<MutationObservation>();
        String fileName = DIR_NAME + "CosmicMutantExport.tsv";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // Check gene
            if (!genes.contains(tokens[0]))
                continue;
            // Check functional impact
            if (pathogenicOnly && !tokens[27].equals("PATHOGENIC"))
                continue;
            // Check Mutation zygosity
            if (!needZygosity && tokens[20].equals("het"))
                continue;
            // Check mutation type
            if (!mutationTypes.contains(tokens[19]))
                continue;
            MutationObservation entry = new MutationObservation();
            entry.setGene(tokens[0]);
            entry.setFathmmType(tokens[27]);
            entry.setMutationType(tokens[19]);
            entry.setSample(tokens[4]);
            entry.setMutationZygosity(tokens[20]);
            // Need to parse the mutation
            entry.setMutation(tokens[18]);
            entry.setCoordinate(parseCoordinate(entry.getMutation()));
            list.add(entry);
        }
        fu.close();
        // Do a sort into a map
        Map<String, List<MutationObservation>> geneToEntries = new HashMap<String, List<MutationObservation>>();
        for (MutationObservation entry : list) {
            List<MutationObservation> entries = geneToEntries.get(entry.getGene());
            if (entries == null) {
                entries = new ArrayList<MutationObservation>();
                geneToEntries.put(entry.getGene(), entries);
            }
            entries.add(entry);
        }
        Map<String, List<MutationObservation>> rtn = new HashMap<String, List<MutationObservation>>();
        for (String gene : geneToEntries.keySet()) {
            List<MutationObservation> entries = mergeEntries(geneToEntries.get(gene));
            rtn.put(gene, entries);
        }
        return rtn;
    }

    private List<MutationObservation> mergeEntries(List<MutationObservation> entries) {
        if (entries.size() < 2)
            return entries;
        // Duplications found in the downloaded file (e.g. BRAF, L505H),
        // use the following to remove duplications
        Map<String, MutationObservation> keyToEntry = new HashMap<String, MutationObservation>();
        for (MutationObservation entry : entries) {
            String key = entry.getGene() + ";" + entry.getMutation() + ";" + entry.getSample();
            if (keyToEntry.containsKey(key))
                continue;
            keyToEntry.put(key, entry);
        }
        return new ArrayList<MutationObservation>(keyToEntry.values());
    }

    @Test
    public void testLoadMutations() throws IOException {
        String gene = "CDK4";
        gene = "BRAF";
        Set<String> genes = new HashSet<String>();
        genes.add(gene);
        Set<String> mutationTypes = new HashSet<String>();
        mutationTypes.add("Substitution - Missense");
        mutationTypes.add("Substitution - Nonsense");
        Map<String, List<MutationObservation>> geneToEntries = loadMutations(genes,
                true,
                true,
                mutationTypes);
        List<MutationObservation> entries = geneToEntries.get(gene);
        System.out.println("Total entries for " + gene + ": " + entries.size());
        Collections.sort(entries, new Comparator<MutationObservation>() {
            public int compare(MutationObservation entry1, MutationObservation entry2) {
                return entry1.getCoordinate().compareTo(entry2.getCoordinate());
            }
        });
        for (MutationObservation entry : entries)
            System.out.println(entry);
    }

    private Integer parseCoordinate(String mutation) {
        // Something like this: p.H132D
        Pattern pattern = Pattern.compile("\\d+"); // Just want to search for the number
        Matcher matcher = pattern.matcher(mutation);
        if (matcher.find()) {
            return new Integer(matcher.group());
        }
        return null;
    }
}
