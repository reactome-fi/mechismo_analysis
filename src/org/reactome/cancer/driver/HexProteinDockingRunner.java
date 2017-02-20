/*
 * Created on Oct 24, 2016
 *
 */
package org.reactome.cancer.driver;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.r3.Interactome3dAnalyzer;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.ProcessRunner;

/**
 * This class is used to run Hex docking program. The implementation here should be ported into a
 * HPC clustering environment for running.
 * @author gwu
 *
 */
public class HexProteinDockingRunner {
    // Directory holding hex related stuff that can be configured
    private String hexResultDir = "results/DriverGenes/Drivers_0816/HexResults/";
    // A shell script to run hex
    private String hexRunCommand = "/Users/gwu/ProgramFiles/Hex/bin/hex"; 
    // For logging
    private String logOutFile = hexResultDir + "hex_out.txt";
    private String logErrFile = hexResultDir + "hex_err.txt";
    private int ncpu = 4;
    
    /**
     * Default constructor.
     */
    public HexProteinDockingRunner() {
    }
    
    public int getNcpu() {
        return ncpu;
    }

    public void setNcpu(int ncpu) {
        this.ncpu = ncpu;
    }

    public String getHexResultDir() {
        return hexResultDir;
    }


    public void setHexResultDir(String hexResultDir) {
        this.hexResultDir = hexResultDir;
    }


    public String getHexRunCommand() {
        return hexRunCommand;
    }

    public void setHexRunCommand(String hexRunCommand) {
        this.hexRunCommand = hexRunCommand;
    }

    public void runHex(String pdbFile1,
                       String pdbFile2,
                       String resultFile) throws IOException {
        String macroFile = generateMacroFile(pdbFile1, pdbFile2, resultFile);
        String[] commands = new String[] {
                hexRunCommand,
                "-ncpu " + ncpu,
                "-nogui",
                "-e " + macroFile 
        };
        
        ProcessRunner runner = new ProcessRunner();
        String[] rtn = runner.runScript(commands);
        
        output(rtn[0], logOutFile);
        output(rtn[1], logErrFile);
    }
    
    private void output(String out,
                        String fileName) throws IOException {
        if (out == null || out.length() == 0)
            return;
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName, true);
        fu.printLine(out);
        fu.close();
    }
    
    private String generateMacroFile(String pdb1,
                                     String pdb2,
                                     String resultFileName) throws IOException {
        File pdbFile1 = new File(pdb1);
        File pdbFile2 = new File(pdb2);
        File resultFile = new File(resultFileName);
        String template = hexResultDir + "hex_docking_template.mac";
        String macroFile = hexResultDir + "hex_docking.mac";
        FileUtility fu = new FileUtility();
        fu.setInput(template);
        fu.setOutput(macroFile);
        String line = null;
        while ((line = fu.readLine()) != null) {
            line = line.replace("$receptor_pdb", pdbFile1.getAbsolutePath());
            line = line.replace("$ligand_pdb", pdbFile2.getAbsolutePath());
            line = line.replace("$result_file", resultFile.getAbsolutePath());
            fu.printLine(line);
        }
        fu.close();
        return macroFile;
    }
    
    /**
     * Batch run docking for one specific reaction.
     * @throws IOException
     */
    @Test
    public void runDockingForReaction() throws IOException {
        String fiFile = "results/DriverGenes/Drivers_0816/FIsInReaction5672965.txt";
        Set<String> fisWithPPIs = new CancerDriverReactomeAnalyzer().loadFIsWithPPIFeature(fiFile);
        System.out.println("Total fis with PPIs: " + fisWithPPIs.size());
        
        String srcDir = "results/DriverGenes/Drivers_0816/interactome3D/reaction_5672965/representative/";
        Map<String, Set<String>> proteinToPDBs = new Interactome3dAnalyzer().loadProteinToPDBFiles(srcDir);
        
        ncpu = 6;
        
        for (String fi : fisWithPPIs) {
            long time1 = System.currentTimeMillis();
            System.out.println("Docking " + fi + "...");
            String[] proteins = fi.split("\t");
            Set<String> pdbs1 = proteinToPDBs.get(proteins[0]);
            Set<String> pdbs2 = proteinToPDBs.get(proteins[1]);
            if (pdbs1 == null || pdbs2 == null) {
                System.out.println("No structure for invovled protein(s)!");
                continue;
            }
            for (String pdb1 : pdbs1) {
                for (String pdb2 : pdbs2) {
                    String resultFile = getResultFile(proteins);
                    System.out.println(resultFile);
                    runHex(pdb1, 
                           pdb2,
                           resultFile);
                }
            }
            long time2 = System.currentTimeMillis();
            System.out.println("Time: " + (time2 - time1) / 60.0d + " seconds.");
//            break;
        }
    }
    
    private String getResultFile(String[] proteins) {
        int count = 0;
        while (true) {
            String fileName = proteins[0] + "_" + proteins[1] + "_" + count + ".txt";
            File file = new File(hexResultDir + fileName);
            if (!file.exists())
                return file.getAbsolutePath();
            count ++;
        }
    }
    
    @Test
    public void testRun() {
        try {
            String dir = "results/DriverGenes/Drivers_0816/interactome3D/reaction_5672965/representative/proteins_1/";
            String pdb1 = dir + "O14511-MDL-O14511.2-1fhg_A.pdb";
            String pdb2 = dir + "P01111-EXP-3con_A.pdb";
            String resultFile = hexResultDir + "O14511_P01111.txt";
            runHex(pdb1, pdb2, resultFile);
        }
        catch(IOException e) {
            e.printStackTrace();
        }
    }
    
}
