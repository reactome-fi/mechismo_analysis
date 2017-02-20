/*
 * Created on Aug 4, 2016
 *
 */
package org.reactome.r3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;

import com.sun.jndi.ldap.ext.StartTlsResponseImpl;

/**
 * This class is used to handle data process related to dbNSFP.
 * @author gwu
 *
 */
public class DbNSFPAnalyzer {
    
    /**
     * Default constructor.
     */
    public DbNSFPAnalyzer() {
    }
    
    @Test
    public void annotateMAFFileWithdbNSFP(String mafFileName,
                                          String annotatedMafFileName,
                                          String scoreFileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(scoreFileName);
        String line = fu.readLine();
        String[] headers = line.split("\t");
//        for (int i = 0; i < headers.length; i++)
//            System.out.println(i + "\t" + headers[i]);
//        fu.close();
//        if (true)
//            return;
        Map<String, Double> keyToScore = new HashMap<String, Double>();
        int scoreIndex = 63; // Use MetaLR_rankscore
        while ((line = fu.readLine()) != null) {
//            System.out.println(line);
            String[] tokens = line.split("\t");
//            2   ref
//            3   alt
//            7   hg19_chr
//            8   hg19_pos(1-based)
//            62  MetaLR_score            
//            63  MetaLR_rankscore
            // It is possible no score is available if a new stop codon is introduced
            if (tokens[scoreIndex].equals("."))
                continue;
            String key = tokens[7] + ":" + tokens[8] + ":" + tokens[2] + ":" + tokens[3];
            Double score = new Double(tokens[scoreIndex]);
            keyToScore.put(key, score);
        }
        fu.close();
        System.out.println("Size of map: " + keyToScore.size());
        
        // Start annotation
        fu.setInput(mafFileName);
        // Output
        fu.setOutput(annotatedMafFileName);
        line = fu.readLine(); // Header
        fu.printLine(line + "\t" + headers[scoreIndex]);
        
        // Get the needed indices
        IndexGroup indexGroup = parseMAFHeaderLine(line);
        
        int counter = 0;
        final Map<String, Integer> keyToCounter = new HashMap<String, Integer>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Double score = null;
            if (tokens[indexGroup.startIndex].equals(tokens[indexGroup.endIndex])) {
                String key = tokens[indexGroup.chrIndex] + ":" + tokens[indexGroup.startIndex] + ":" + tokens[indexGroup.refIndex] + ":" + tokens[indexGroup.altIndex];
                score = keyToScore.get(key);
                if (score != null) {
                    Integer counts = keyToCounter.get(key);
                    if (counts == null)
                        keyToCounter.put(key, 1);
                    else
                        keyToCounter.put(key, ++counts);
                }
            }
            fu.printLine(line + "\t" + (score == null ? "" : score));
            if (score != null)
                counter ++;
        }
        fu.close();
        System.out.println("Total annotated lines: " + counter);
        // Print out top used keys: for debugging
        List<String> keyList = new ArrayList<String>(keyToCounter.keySet());
        Collections.sort(keyList, new Comparator<String>() {
            public int compare(String key1, String key2) {
                Integer count1 = keyToCounter.get(key1);
                Integer count2 = keyToCounter.get(key2);
                if (count1 == count2)
                    return key1.compareTo(key2);
                else
                    return count2.compareTo(count1);
            }
        });
        for (int i = 0; i < keyList.size(); i++) {
            String key = keyList.get(i);
            counter = keyToCounter.get(key);
            if (counter == 1)
                break;
            System.out.println((i + 1) + "\t" + key + "\t" + keyToCounter.get(key));
        }
    }
    
    /**
     * Use this method to generate an input file for the Java class provided by dbNSFP, search_dbNSFP32a. An example on
     * how to use this method is as following:
     * java -Xmx8G search_dbNSFP32a -i ../../ICGC/2016_04/Barcelona_consensus.filter.genes.dbNSFP.input.txt 
     *                              -o Barcelona_consensus.filter.genes.dbNSFP.output.txt 
     *                              -v hg19
     * @param mafFileName
     * @param outFileName
     * @throws IOException
     */
    @Test
    public void generateInputFordbNSFPFromMAF(String mafFileName,
                                              String outFileName) throws IOException {
        FileUtility fu = new FileUtility();
        // Generate a file format that can be recognized by dbNSFP Java class
        fu.setInput(mafFileName);
        String outputFileName = outFileName;
        fu.setOutput(outputFileName);
        String line = fu.readLine();
        // Get the needed indices
        IndexGroup indexGroup = parseMAFHeaderLine(line);
        // We don't need the title
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (!tokens[indexGroup.startIndex].equals(tokens[indexGroup.endIndex])) {
                // Most likely they are deletion or insertion.
                // Escape them!
                //System.out.println(line);
                continue;
            }
            fu.printLine(tokens[indexGroup.chrIndex] + "\t" + // Chromosome
                         tokens[indexGroup.startIndex] + "\t" + // Position
                         tokens[indexGroup.refIndex] + "\t" + // Reference
                         tokens[indexGroup.altIndex] + "\t"); // alternative
        }
        fu.close();
    }
    
    private IndexGroup parseMAFHeaderLine(String line) {
     // Get the needed indices
        String[] headers = line.split("\t");
        int chrIndex = -1, startIndex = -1, refIndex = -1, altIndex = -1, endIndex = -1;
        int tumorSeqIndex = -1;
        for (int i = 0; i < headers.length; i++) {
            if (headers[i].equals("Chromosome"))
                chrIndex = i;
            else if (headers[i].equals("Start_Position"))
                startIndex = i;
            else if (headers[i].equals("End_Position"))
                endIndex = i;
            else if (headers[i].equals("Reference_Allele"))
                refIndex = i;
            else if (headers[i].equals("Alternate_Allele"))
                altIndex = i;
            else if (headers[i].equals("Tumor_Seq_Allele2")) // Used by TCGA MAF
                tumorSeqIndex = i;
        }
        if (altIndex < 0)
            altIndex = tumorSeqIndex;
        IndexGroup indexGroup = new IndexGroup();
        indexGroup.chrIndex = chrIndex;
        indexGroup.startIndex = startIndex;
        indexGroup.endIndex = endIndex;
        indexGroup.refIndex = refIndex;
        indexGroup.altIndex = altIndex;
        indexGroup.endIndex = endIndex;
        return indexGroup;
    }
    
    private class IndexGroup {
        int chrIndex = -1, startIndex = -1, refIndex = -1, altIndex = -1, endIndex = -1;
    }
    
}
