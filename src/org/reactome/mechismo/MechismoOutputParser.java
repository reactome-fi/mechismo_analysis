package org.reactome.mechismo;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.fi.util.InteractionUtilities;

/**
 * The actual parser used to parse output downloaded from mechsimo.
 * @author wug
 *
 */
public class MechismoOutputParser {
    
    public MechismoOutputParser() {
        
    }
    
    @Test
    public void checkOutput() throws IOException {
        String fileName = "datasets/Mechismo/TCGA/TCGA_mech_output.tsv";
        // The following is used to check the used PDB structures
        Map<String, Set<String>> interactionToStructures = new HashMap<>();
        Files.lines(Paths.get(fileName))
             .skip(1)
             .forEach(line -> {
                 String[] tokens = line.split("\t");
                 if (tokens.length < 33)
                     return;
                 String interactions = InteractionUtilities.generateFIFromGene(tokens[0], tokens[18]);
                 interactionToStructures.compute(interactions, (key, set) -> {
                     if (set == null)
                         set = new HashSet<>();
                     set.add(tokens[32]);
                     return set;
                 });
             });
        
        System.out.println("Total interactions in interactionsToStructures: " + interactionToStructures.size());
        
        interactionToStructures.forEach((key, set) -> {
            if (set.size() > 1)
                System.out.println(key + "\t" + set);
        });
        
//        // The following is used to check how many interactions are covered by the output file
//        Set<String> interactions = Files.lines(Paths.get(fileName))
//                                        .skip(1)
//                                        .map(line -> line.split("\t"))
//                                        .filter(tokens -> tokens.length > 18)
//                                        .map(tokens -> tokens[0] + "\t" + tokens[18])
//                                        .collect(Collectors.toSet());
//        System.out.println("Total interactions: " + interactions.size());
//        // Want to do a sorting
//        Set<String> sorted = interactions.stream()
//                                         .map(interaction -> {
//                                             String[] tokens = interaction.split("\t");
//                                             String newFI = InteractionUtilities.generateFIFromGene(tokens[0], tokens[1]);
//                                             return newFI;
//                                         })
//                                         .collect(Collectors.toSet());
//        System.out.println("Sorted interactions: " + sorted.size());
    }

}
