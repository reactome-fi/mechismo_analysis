/*
 * Created on Sep 6, 2016
 *
 */
package org.reactome.r3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;

/**
 * This class is used to perform analysis on the cosmic data.
 * @author gwu
 *
 */
public class CosmicAnalyzer {
    private final String DIR_NAME = "datasets/COSMIC/v78/";
    
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
     * @param entries
     * @param coordinates
     * @return
     */
    public List<CosmicEntry> filterEntries(List<CosmicEntry> entries,
                                           Collection<Integer> coordinates) {
        List<CosmicEntry> rtn = new ArrayList<CosmicAnalyzer.CosmicEntry>();
        for (CosmicEntry entry : entries) {
            if (coordinates.contains(entry.coordinate))
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
     * @param genes
     * @return
     * @throws IOException
     */
    public Map<String, List<CosmicEntry>> loadMutations(Set<String> genes) throws IOException {
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
        String line = fu.readLine();
        Set<String> mutationTypes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            mutationTypes.add(tokens[19]);
        }
        fu.close();
        System.out.println("Total mutation types: " + mutationTypes.size());
        for (String type : mutationTypes)
            System.out.println(type);
    }
    
    /**
     * Load entries in COSMIC for a set of genes.
     * @param genes
     * @param pathogenicOnly
     * @param needZygosity
     * @param mutationTypes
     * @return
     * @throws IOException
     */
    public Map<String, List<CosmicEntry>> loadMutations(Set<String> genes,
                                                        boolean pathogenicOnly,
                                                        boolean needZygosity,
                                                        Set<String> mutationTypes) throws IOException {
        List<CosmicEntry> list = new ArrayList<CosmicAnalyzer.CosmicEntry>();
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
            CosmicEntry entry = new CosmicEntry();
            entry.gene = tokens[0];
            entry.fathmmType = tokens[27];
            entry.mutationType = tokens[19];
            entry.sample = tokens[4];
            entry.mutationZygosity = tokens[20];
            // Need to parse the mutation
            entry.mutation = tokens[18];
            entry.coordinate = parseCoordinate(entry.mutation);
            list.add(entry);
        }
        fu.close();
        // Do a sort into a map
        Map<String, List<CosmicEntry>> geneToEntries = new HashMap<String, List<CosmicEntry>>();
        for (CosmicEntry entry : list) {
            List<CosmicEntry> entries = geneToEntries.get(entry.getGene());
            if (entries == null) {
                entries = new ArrayList<CosmicAnalyzer.CosmicEntry>();
                geneToEntries.put(entry.gene, entries);
            }
            entries.add(entry);
        }
        Map<String, List<CosmicEntry>> rtn = new HashMap<String, List<CosmicEntry>>();
        for (String gene : geneToEntries.keySet()) {
            List<CosmicEntry> entries = mergeEntries(geneToEntries.get(gene));
            rtn.put(gene, entries);
        }
        return rtn;
    }
    
    private List<CosmicEntry> mergeEntries(List<CosmicEntry> entries) {
        if (entries.size() < 2)
            return entries;
        // Duplications found in the downloaded file (e.g. BRAF, L505H),
        // use the following to remove duplications
        Map<String, CosmicEntry> keyToEntry = new HashMap<String, CosmicAnalyzer.CosmicEntry>();
        for (CosmicEntry entry : entries) {
            String key = entry.getGene() + ";" + entry.getMutation() + ";" + entry.getSample();
            if (keyToEntry.containsKey(key))
                continue;
            keyToEntry.put(key, entry);
        }
        return new ArrayList<CosmicEntry>(keyToEntry.values());
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
        Map<String, List<CosmicEntry>> geneToEntries = loadMutations(genes,
                                                                     true,
                                                                     true,
                                                                     mutationTypes);
        List<CosmicEntry> entries = geneToEntries.get(gene);
        System.out.println("Total entries for " + gene + ": " + entries.size());
        Collections.sort(entries, new Comparator<CosmicEntry>() {
            public int compare(CosmicEntry entry1, CosmicEntry entry2) {
                return entry1.coordinate.compareTo(entry2.coordinate);
            }
        });
        for (CosmicEntry entry : entries)
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
    
    public static class CosmicEntry { 
        
        private String gene;
        private Integer coordinate;
        private String mutation;
        private String mutationType;
        private String fathmmType;
        private String mutationZygosity;
        private String sample;
        
        public CosmicEntry() {
        }

        public String getSample() {
            return sample;
        }

        public void setSample(String sample) {
            this.sample = sample;
        }

        public String getMutationZygosity() {
            return mutationZygosity;
        }

        public void setMutationZygosity(String mutationZygosity) {
            this.mutationZygosity = mutationZygosity;
        }

        public String getGene() {
            return gene;
        }

        public void setGene(String gene) {
            this.gene = gene;
        }

        public Integer getCoordinate() {
            return coordinate;
        }

        public void setCoordinate(Integer coordinate) {
            this.coordinate = coordinate;
        }

        public String getMutation() {
            return mutation;
        }

        public void setMutation(String mutation) {
            this.mutation = mutation;
        }

        public String getMutationType() {
            return mutationType;
        }

        public void setMutationType(String mutationType) {
            this.mutationType = mutationType;
        }

        public String getFathmmType() {
            return fathmmType;
        }

        public void setFathmmType(String fathmmType) {
            this.fathmmType = fathmmType;
        }
        
        public String toString() {
            return gene + " " + mutation + " " + fathmmType + " " + sample + " " + mutationZygosity;
        }
        
    }
    
}
