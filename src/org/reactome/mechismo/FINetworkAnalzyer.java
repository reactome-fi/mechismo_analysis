package org.reactome.mechismo;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.junit.Test;
import org.reactome.fi.util.FileUtility;
import org.reactome.mechismo.model.AnalysisResult;
import org.reactome.mechismo.model.CancerType;
import org.reactome.mechismo.model.Interaction;
import org.reactome.mechismows.MechismowsReader;
import org.reactome.r3.util.InteractionUtilities;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

/**
 * This class is used to analyze the mechismo results in the perspective of the 
 * Reactome FI network, composed FIs extracted from Reactome only.
 * @author wug
 *
 */
public class FINetworkAnalzyer {

    public FINetworkAnalzyer() {
    }

    private Set<String> loadReactomeFIs() throws IOException {
        String fileName = "results/ProteinFIsInReactions_032017.txt";
        try (Stream<String> lines = Files.lines(Paths.get(fileName))) {
            Set<String> fis = lines.map(line -> {
                String[] tokens = line.split("\t");
                String fi = InteractionUtilities.generateFIFromGene(tokens[2], tokens[3]);
                return fi;
            })
                    .collect(Collectors.toSet());
            return fis;
        }
    }
    
    /**
     * Generate a tab delimited file having FDR values for individual cancer types.
     * @throws Exception
     */
    @Test
    public void sortInteractions() throws Exception {
        Set<String> reactomeFIs = loadReactomeFIs();
        // Output all FIs in mechismo output file
        String fileName = "results/FI_Cancer_FDR_050918.txt";
        // Output FIs having analysis results attached (no protein/chemicals)
        fileName = "results/FI_Cancer_FDR_Filtered_050918.txt";
        fileName = "results/FI_Cancer_FDR_Filtered_051018.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        MechismowsReader reader = new MechismowsReader();
        AnnotationConfigApplicationContext context = reader.getContext();
        List<CancerType> cancerTypes = reader.loadCancerTypes(context);
        List<String> sortedCancerNames = sortCancerTypes(cancerTypes);
        List<Interaction> interactions = reader.loadInteractions(context);
        double smallestFDR = getSmallestNonZeroFDR(interactions);
        double largestPScore = -Math.log10(smallestFDR / 10.0d); // Treat as the largest p-score
        // Perform filtering
        interactions = interactions.stream()
                                   .filter(i -> i.getAnalysisResults() != null) // Make sure analysis results there
                                   .filter(i -> !i.getName().contains("[CHEM:")) // Make sure they are protein FIs
                                   .collect(Collectors.toList());
        fu.printLine(generateHeader(sortedCancerNames));
        interactions.forEach(interaction -> {
            String text = outputResultsAsString(interaction, sortedCancerNames, reactomeFIs, largestPScore);
            if (text == null || text.length() == 0)
                return;
            try {
                fu.printLine(text); 
            }
            catch(IOException e) {
                e.printStackTrace();
            }
        });
        context.close();
        fu.close();
    }
    
    private Double getSmallestNonZeroFDR(List<Interaction> interactions) {
        Double smallest = Double.MAX_VALUE;
        for (Interaction interaction : interactions) {
            Set<AnalysisResult> results = interaction.getAnalysisResults();
            for (AnalysisResult result : results) {
                if (result.getFdr() > 0.0d && result.getFdr() < smallest)
                    smallest = result.getFdr();
            }
        }
        return smallest;
    }
    
    private List<String> sortCancerTypes(List<CancerType> types) {
        List<String> rtn = types.stream()
                                .map(type -> type.getAbbreviation())
                                .sorted((t1, t2) -> {
                                    if (t1.equals("PANCAN"))
                                        return -1;
                                    if (t2.equals("PANCAN"))
                                        return 1;
                                    return t1.compareTo(t2);
                                })
                                .collect(Collectors.toList());
        return rtn;
    }
    
    private String generateHeader(List<String> cancerTypes) {
        StringBuilder builder = new StringBuilder();
        builder.append("Partner1\tPartner2\tInReactome\tSignificantCancers\tRankScore\tPanCancer_PScore");
        cancerTypes.forEach(type -> builder.append("\t").append(type));
        return builder.toString();
    }
    
    /**
     * A rank score is calculated as sum of negative log fdr value. If fdr == 0.0,
     * the largest score will be used
     * @param results
     * @return
     */
    private double calculateRankScore(Set<AnalysisResult> results, double largestScore) {
        double score = 0.0d;
        for (AnalysisResult result : results) {
            Double fdr = result.getFdr();
            if (fdr.equals(0.0d))
                score += largestScore;
            else
                score -= Math.log(fdr);
        }
        return score;
    }
    
    private int calculateSignificantCancers(Set<AnalysisResult> results) {
        int counter = 0;
        for (AnalysisResult result : results) {
            if (result.getFdr() <= 0.05d)
                counter ++;
        }
        return counter;
    }
    
    private double calculatePancanPScore(Map<String, AnalysisResult> typeToResult, double largestPScore) {
        AnalysisResult result = typeToResult.get("PANCAN");
        if (result == null)
            return 0.0d;
        Double fdr = result.getFdr();
        if (fdr.equals(0.0d))
            return largestPScore;
        return -Math.log10(fdr);
    }
    
    private String outputResultsAsString(Interaction interaction,
                                         List<String> cancerTypes,
                                         Set<String> reactomeFIs,
                                         double largestPScore) {
        Set<AnalysisResult> results = interaction.getAnalysisResults();
        if (results == null || results.size() == 0)
            return null;
        StringBuilder builder = new StringBuilder();
        builder.append(interaction.getName());
        builder.append("\t").append(reactomeFIs.contains(interaction.getName()) ? 1 : 0);
        int significantCancers = calculateSignificantCancers(results);
        builder.append("\t").append(significantCancers);
        double rankScore = calculateRankScore(results, largestPScore);
        builder.append("\t").append(rankScore);
        // CancerType objects may be loaded from different hibernate session,
        // therefore abbreviations are used as keys.
        Map<String, AnalysisResult> typeToResult = new HashMap<>();
        results.forEach(result -> typeToResult.put(result.getCancerType().getAbbreviation(),
                                                   result));
        double pancanPScore = calculatePancanPScore(typeToResult, largestPScore);
        builder.append("\t").append(pancanPScore);
        cancerTypes.forEach(type -> {
            AnalysisResult result = typeToResult.get(type);
            builder.append("\t");
            if (result == null)
                builder.append("");
            else
                builder.append(result.getFdr());
        });
        return builder.toString();
    }
    
}
