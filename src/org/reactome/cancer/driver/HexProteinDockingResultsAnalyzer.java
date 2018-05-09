/*
 * Created on Oct 26, 2016
 *
 */
package org.reactome.cancer.driver;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Vector3d;

import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.RotationConvention;
import org.apache.commons.math3.geometry.euclidean.threed.RotationOrder;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.io.FileConvert;
import org.junit.Test;
import org.reactome.r3.Interactome3dAnalyzer;
import org.reactome.r3.util.FileUtility;

/**
 * This class is used to process the results from Hex docking software.
 * @author gwu
 *
 */
public class HexProteinDockingResultsAnalyzer {
    private final String RESULT_DIR = "results/DriverGenes/Drivers_0816/HexResults/";
    private FileUtility fu = new FileUtility();
    
    /**
     * Default constructor.
     */
    public HexProteinDockingResultsAnalyzer() {
    }
    
    private String[] getSwissProtIds(String fileName) {
        File file = new File(fileName);
        String shortName = file.getName();
        String[] tokens = shortName.split("_");
        return new String[]{tokens[0], tokens[1]};
    }
    
    public void analyzeTransformResult(String fileName) throws IOException, StructureException {
        String[] swissProtIds = getSwissProtIds(fileName);
        fu.setInput(fileName);
        String line = null;
        Structure receptor = null;
        Structure ligand = null;
        Interactome3dAnalyzer helper = new Interactome3dAnalyzer();
        helper.setUpAtomCache();
        Map<Structure, Double> ligandToEnergy = new HashMap<Structure, Double>();
        System.out.println("Solution\tETotal\tBoltzmann\t" + swissProtIds[0] + "\t" + swissProtIds[1] + "\tTotal");
        int maxiumContact = 0;
        Structure maxiumContactStructure = null;
        while ((line = fu.readLine()) != null) {
            if (line.trim().length() == 0)
                continue;
            else if (line.startsWith("# ReceptorSourceFile")) {
                String receptorSourceFile = getFileName(line);
                receptor = StructureIO.getStructure(receptorSourceFile);
                receptor.getChains().get(0).setSwissprotId(swissProtIds[0]);
            }
            else if (line.startsWith("# LigandSourceFile")) {
                String ligandSourceFile = getFileName(line);
                ligand = StructureIO.getStructure(ligandSourceFile);
                ligand.getChains().get(0).setSwissprotId(swissProtIds[1]);
                outputContacts("Original", 
                               null,
                               swissProtIds,
                               analyzeContact(receptor,
                                              ligand,
                                              helper));
            }
            else if (line.startsWith("#"))
                continue; // Other annotation lines we don't care
            else {
                // One line for one ligand solution
                String[] tokens = line.split("( )+");
                // Format as the following:
                // 0       1        2   3      4 5 6 7     8    9     10             11
                // Cluster Solution RMS Etotal x y z alpha beta gamma receptor-model ligand-model
                Double etotal = null;
                if (!tokens[3].equals("nan"))
                    etotal = new Double(tokens[3]);
                double[][] rotation = createRotationMatrix(Double.parseDouble(tokens[7]),
                                                           Double.parseDouble(tokens[8]),
                                                           Double.parseDouble(tokens[9]));
                Structure solution = ligand.clone();
                Calc.rotate(solution, 
                            rotation);
                Vector3d translation = new Vector3d(new double[] {
                        Double.parseDouble(tokens[4]),
                        Double.parseDouble(tokens[5]),
                        Double.parseDouble(tokens[6])
                });
                Calc.translate(solution, translation);
                
                // Get the final structure
                int totalContact = outputContacts(tokens[1],
                                                  etotal,
                                                  swissProtIds,
                                                  analyzeContact(receptor, solution, helper));
                
                Structure complex = createComplexStructure(receptor, solution);
//                viewStructure(complex);
                
                if (totalContact > maxiumContact) {
                    maxiumContact = totalContact;
                    maxiumContactStructure = complex;
                }
//                
//                saveStructure(complex, fileName, tokens[1]);
//                break;
            }
        }
        fu.close();
        
//        BiojavaJmol jmolPane = new BiojavaJmol();
//        jmolPane.setStructure(maxiumContactStructure);
//        saveStructure(maxiumContactStructure, fileName, "2000");
    }
    
    
    
    private void viewStructure(Structure complex) {
//        BiojavaJmol jmolPane = new BiojavaJmol();
//        jmolPane.setStructure(complex);
    }
    
    private void saveStructure(Structure complex,
                               String dockFileName,
                               String solution) throws IOException {
        FileConvert convert = new FileConvert(complex);
        String pdb = convert.toPDB();
        
        // Get the file name
        File file = new File(dockFileName);
        String shortFileName = file.getName();
        int index = shortFileName.indexOf(".");
        String pdbFileName = shortFileName.substring(0, index) + "_" + solution + ".pdb";
        File pdbFile = new File(file.getParent(), pdbFileName);
        
        FileUtility fu = new FileUtility();
        fu.setOutput(pdbFile.getAbsolutePath());
        fu.printLine(pdb);
        fu.close();
    }
                               

    private Structure createComplexStructure(Structure receptor, Structure ligand) {
        Structure complex = new StructureImpl();
        Chain receptorChain = receptor.getChains().get(0);
        // Give chains different ids. Otherwise,
        // they will be treated as one single chain.
        receptorChain.setChainID("A");
        complex.addChain(receptorChain);
        Chain ligandChain = ligand.getChains().get(0);
        ligandChain.setChainID("B");
        complex.addChain(ligandChain);
        return complex;
    }
    
    private Map<Chain, List<Integer>> analyzeContact(Structure receptor,
                                                     Structure ligand,
                                                     Interactome3dAnalyzer helper) throws StructureException {
        // Assume there is only one chain in each structure
        Chain receptorChain = receptor.getChains().get(0);
        Chain ligandChain = ligand.getChains().get(0);
        Map<Chain, List<Integer>> chainToContact = helper.extractContacts(receptorChain, ligandChain);
        return chainToContact;
    }
    
    private double calculateBoltzmannDistribtion(double etotal) {
        // The unit for etotal should be kj/mol. We need to convert it into
        // j for one molecule
        // Boltzmann distribution: https://en.wikipedia.org/wiki/Boltzmann_distribution
        // which is defined as exp(-E/kT). Here k = R/Na, R is gas constant, 8.31 JK-1mol-1.
        // Therefore, Boltzman can be calculate as: exp(etotal * 1000 / (R * T)). We use human
        // temperature T = 273 + 37
        double bp = Math.exp(-etotal * 1000.0 / (8.31 * (273 + 37))); 
        return bp;
    }
    
    private int outputContacts(String title,
                               Double etotal,
                               String[] swissProts,
                               Map<Chain, List<Integer>> chainToContacts) {
        // Switch from chain to uniprot to make sure the order is correct in the output
        Map<String, Integer> proteinToSize = new HashMap<String, Integer>();
        StringBuilder builder = new StringBuilder();
        builder.append(title);
        builder.append("\t").append(etotal);
        // Calculate Boltamann distribution
        if (etotal == null)
            builder.append("\tNA");
        else {
            double bp = calculateBoltzmannDistribtion(etotal);
            builder.append("\t").append(bp);
        }
        for (Chain chain : chainToContacts.keySet()) {
            List<Integer> contacts = chainToContacts.get(chain);
            proteinToSize.put(chain.getSwissprotId(), contacts.size());
        }
        int total = 0;
        for (String protein : swissProts) {
            Integer contacts = proteinToSize.get(protein);
            if (contacts == null)
                contacts = 0;
            builder.append("\t").append(contacts);
            total += contacts;
        }
        builder.append("\t").append(total);
        System.out.println(builder.toString());
        return total;
    }
    
    private String getFileName(String line) {
        int index = line.indexOf(":");
        String fileName = line.substring(index + 1).trim();
        return fileName;
    }
    
    /**
     * Use Rotation class in apache math to create a rotation matrix that can
     * be used by biojava structure API.
     * @param alpha
     * @param beta
     * @param gamma
     * @return
     */
    private double[][] createRotationMatrix(double alpha,
                                            double beta,
                                            double gamma) {
        // Based on the document and the comparison results from hex saving between
        // matrix and transform, RotationOrder.ZYZ should be used, which means
        // rotate around axis z first by alpha, then around y by beta, and around
        // z again for gamma.
        Rotation rotation = new Rotation(RotationOrder.ZYZ,
                                         RotationConvention.VECTOR_OPERATOR,
                                         alpha,
                                         beta,
                                         gamma);
        return rotation.getMatrix();
    }
    
    @Test
    public void testAnalyzeTransformResults() throws Exception {
        String fileName = RESULT_DIR + "Reaction_5672965_102616/P01116_Q05397_1.txt";
        fileName = RESULT_DIR + "Reaction_5672965_102616/P01112_Q9Y4G8_0.txt";
//        String fileName = "/Users/gwu/Documents/temp/docking_results_transform_1.hex";
        analyzeTransformResult(fileName);
    }
    
    public static void main(String[] args) throws Exception {
        HexProteinDockingResultsAnalyzer analyzer = new HexProteinDockingResultsAnalyzer();
        analyzer.testAnalyzeTransformResults();
    }
    
    @Test
    public void testRotation() {
        double[] angles = new double[] {
//                3.799728, 0.935227, 3.129773
                3.884926, 0.902009, 2.972368
        };
        double[][] matrix = createRotationMatrix(angles[0],
                                                 angles[1],
                                                 angles[2]);
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++)
                System.out.print(matrix[i][j] + "    ");
            System.out.println();
        }
    }
}
