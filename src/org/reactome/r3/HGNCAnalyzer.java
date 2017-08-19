/*
 * Created on May 24, 2017
 *
 */
package org.reactome.r3;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.junit.Test;

/**
 * @author gwu
 *
 */
public class HGNCAnalyzer {
    private String dirName = "datasets/HGNC/";
    private String geneFileName = dirName + "gene_with_protein_product_052417.txt";
    
    /**
     * Default constructor.
     */
    public HGNCAnalyzer() {
    }
    
    public Set<String> loadHumanGenes() throws IOException {
        Set<String> genes = Files.lines(Paths.get(geneFileName))
             .skip(1)
             .map(line -> line.split("\t")[1])
             .collect(Collectors.toSet());
        return genes;
    }
    
    @Test
    public void testLoadHumanGeneToUniProtId() throws IOException {
        Map<String, String> geneToId = loadHumanGeneToUniProtId();
        System.out.println("Total genes: " + geneToId.size());
    }
    
    @Test
    public void testLoadHumanUniProtToGene() throws IOException {
        Map<String, String> idToGene = loadHumanGeneToUniProtId();
        System.out.println("Total UniProt ids: " + idToGene.size());
    }
    
    public Map<String, String> loadHumanUniProtToGene() throws IOException {
        Map<String, String> idToGene = new HashMap<>();
        Files.lines(Paths.get(geneFileName))
             .skip(1)
             .map(line -> line.split("\t"))
             .filter(tokens -> tokens.length >= 26 && tokens[25].length() > 0)
             .forEach(tokens -> {
                 String gene = tokens[1];
                 String id = tokens[25];
                 if (idToGene.containsKey(id))
                     throw new IllegalStateException(id + " has more than one gene!");
                 idToGene.put(id, gene);
             });
        return idToGene;
    }
    
    public Map<String, String> loadHumanGeneToUniProtId() throws IOException {
        Map<String, String> geneToId = new HashMap<>();
        Files.lines(Paths.get(geneFileName))
        .skip(1)
        .map(line -> line.split("\t"))
        .forEach(tokens -> {
            if (tokens.length < 26)
                return;
            String gene = tokens[1];
            String id = tokens[25];
            if (id.length() == 0)
                return;
            if (geneToId.containsKey(gene))
                System.err.println("Error: " + gene + " is duplicated!");
            geneToId.put(gene, id);
        });
        return geneToId;
    }
    
    @Test
    public void testLoadHumanGenes() throws IOException {
        Set<String> genes = loadHumanGenes();
        System.out.println("This number should be 19065: " + genes.size());
    }
    
}
