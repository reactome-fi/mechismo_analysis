package org.reactome.main;

import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.cancer.driver.CancerDriverReactomeAnalyzer;
import org.reactome.cancer.driver.Interactome3dDriverAnalyzer;
import org.reactome.cancer.driver.MechismoAnalyzer;
import org.reactome.r3.util.Configuration;

import java.util.Scanner;

/**
 * Created by burkhart on 3/1/17.
 */
public class ProcessDatasetsDriver {

    public ProcessDatasetsDriver() {
    }

    private MySQLAdaptor getReactomeDBA() throws Exception {
        Configuration configuration = Configuration.getConfiguration();
        MySQLAdaptor dba = configuration.getReactomeDBA();
        if (dba != null)
            return dba;
        // Try to use information from System.in if we cannot get a pre-configured one
        System.out.print("MySQL Username: ");
        Scanner scanner = new Scanner(System.in);
        String un = scanner.next();

        System.out.print("MySQL Password: ");
        String pw = scanner.next();
        scanner.close();

        dba = new MySQLAdaptor("localhost",
                "reactome_59_plus_i",
                un,
                pw);
        return dba;
    }

    @Test
    public void testFindInteractionsWithMutatedInterfaces() {
        run("1");
    }


    void run() {
        System.out.print("Interface Enrichment 0\n" +
                "Mechismo/Reactome Overlay 1\n" +
                "Mechismo/Reactome Interface Enrichment 2\n" +
                "Prepare Heatmap Data 3\n" +
                "Compare Known Drivers 4\n" +
                "Map Mechismo Reactome Reactions 5\n" +
                "Analyze Mechismo Interface Co-occurrence 6\n" +
                "Cancel <enter>\n");
        Scanner scanner = new Scanner(System.in);
        String ex = scanner.nextLine();
        scanner.close();
        run(ex);
    }

    private void run(String ex) {
        Interactome3dDriverAnalyzer interactome3dDriverAnalyzer = new Interactome3dDriverAnalyzer();
        try {
            MySQLAdaptor dba = getReactomeDBA();

            CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer = new CancerDriverReactomeAnalyzer();
            cancerDriverReactomeAnalyzer.setDBA(dba);
            if (ex.equals("")) {
                System.out.println("Execution Cancelled.");
            } else if (Integer.parseInt(ex) == 0) {
                interactome3dDriverAnalyzer.findInteractionsWithMutatedInterfaces(cancerDriverReactomeAnalyzer,
                        "datasets/firehose_all/",
                        // The MAF filenames look like:
                        // TCGA-AY-4070-01.hg19.oncotator.hugo_entrez_remapped.maf.txt
                        // TCGA-AY-4071-01.hg19.oncotator.hugo_entrez_remapped.maf.txt
                        "^.+\\.maf\\.txt$",
                        "datasets/interactome3d/pdb_structures/",
                        "results/interactions_with_mutated_interfaces.csv");
            } else if (Integer.parseInt(ex) == 1) {
                interactome3dDriverAnalyzer.findMechismoInteractionsInReactome(cancerDriverReactomeAnalyzer,
                        "datasets/Mechismo/COSMICv74_somatic_noSNPs_GWS_mechismo_output.tsv",
                        "^[^\\s]+\\s+(?<prot1>[^\\s]+)\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+(?<mech>[^\\s]+)\\s+[0-9a-zA-Z-]+\\s+(?<prot2>[0-9a-zA-Z-]+).*$",
                        "datasets/Mechismo/cancer_types/large_intestine+colon+carcinoma+adenocarcinoma/sample_id.txt",
                        "^(?<patientID>[^\\s]+)\\s+[^\\s]+\\s+(?<prot1>[^\\s]+)\\s+[^\\s]+\\s+[^\\s\\[\\]]+\\s+(?<prot2>[^\\s]+)\\s+(?<mech>[^\\s]+).*",
                        "results/interactions_with_mechismo_enrichment.csv");
            } else if (Integer.parseInt(ex) == 2) {
                interactome3dDriverAnalyzer.calculateMechismoInteractionCorrelationWithInterfaceMutationEnrichment(cancerDriverReactomeAnalyzer,
                        "datasets/Mechismo/COSMICv74_somatic_noSNPs_GWS_mechismo_output.tsv",
                        "^[^\\s]+\\s+(?<prot1>[^\\s]+)\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+(?<mech>[^\\s]+)\\s+[0-9a-zA-Z-]+\\s+(?<prot2>[0-9a-zA-Z-]+).*$",
                        "datasets/Mechismo/cancer_types/large_intestine+colon+carcinoma+adenocarcinoma/sample_id.txt",
                        "^(?<patientID>[^\\s]+)\\s+[^\\s]+\\s+(?<prot1>[^\\s]+)\\s+[^\\s]+\\s+[^\\s\\[\\]]+\\s+(?<prot2>[^\\s]+)\\s+(?<mech>[^\\s]+).*",
                        "results/interactions_with_mechismo_enrichment.csv",
                        "datasets/firehose_all/",
                        // The MAF filenames look like:
                        // TCGA-AY-4070-01.hg19.oncotator.hugo_entrez_remapped.maf.txt
                        // TCGA-AY-4071-01.hg19.oncotator.hugo_entrez_remapped.maf.txt
                        "^.+\\.maf\\.txt$",
                        "datasets/interactome3d/pdb_structures/",
                        "results/interactions_with_mutated_interfaces.csv",
                        "results/interactions_with_mechismo_enrichment_and_mutated_interfaces.csv");
            } else if (Integer.parseInt(ex) == 3) {
                interactome3dDriverAnalyzer.prepareHeatmapData(cancerDriverReactomeAnalyzer,
                        "datasets/interactome3d/pdb_structures/",
                        "results/heatmapData/",
                        "datasets/firehose_data/all_oncotated_calls",
                        "datasets/Mechismo/cancer_types/");
            } else if (Integer.parseInt(ex) == 4) {
                interactome3dDriverAnalyzer.compareKnownDrivers("results/knownDriverData/",
                        "datasets/firehose_data/all_oncotated_calls",
                        "datasets/guanming_known_drivers.txt");
            } else if (Integer.parseInt(ex) == 5) {
                String[] cancerTypes = new String[]{
                        "BRCA",
                        "KI",
                        "GBM",
                        "OV",
                        "LUAD",
                        "UCEC",
                        "KIRC",
                        "HNSC",
                        "LGG",
                        "THCA",
                        "LUSC",
                        "PRAD",
                        "COAD",
                        "STAD",
                        "LIHC",
                        "BLCA",
                        "SKCM",
                        "CESC",
                        "KIRP",
                        "SARC",
                        "LAML",
                        "PAAD"
                };

                Integer[] onePcts = new Integer[]{
                        12,
                        10,
                        7,
                        7,
                        6,
                        6,
                        6,
                        6,
                        6,
                        6,
                        6,
                        6,
                        5,
                        5,
                        5,
                        5,
                        4,
                        4,
                        3,
                        3,
                        3,
                        2
                };

                for(int i = 0; i < cancerTypes.length; i++) {
                    new MechismoAnalyzer().mapReactomeReactions(
                            cancerDriverReactomeAnalyzer,
                            "/home/burkhart/Software/Ogmios/datasets/Mechismo/TCGA_mech_output.tsv",
                            "/home/burkhart/Software/Ogmios/datasets/ReactionNetwork_070517.txt",
                            "/home/burkhart/Software/Ogmios/datasets/Mechismo/MechismoSamplesToReactions_103017.txt",
                            null,//"/home/burkhart/Software/Ogmios/datasets/Mechismo/tcga_mechismo_stat_pancancer.tsv",
                            "/home/burkhart/Software/Ogmios/results/Mechismo/",
                            cancerTypes[i], //output filename prefix
                            cancerTypes[i], //set cancer type = null for pancancer
                            1,
                            100,
                            0.0,
                            1.0,
                            true, // Ignore Dependent
                            false,
                            false,
                            false,
                            false,
                            "No", // Rxn Filter
                            onePcts[i]
                    );
                }
            }else if(Integer.parseInt(ex) == 6){
                new MechismoAnalyzer().analyzeInterfaceCooccurrence(
                        cancerDriverReactomeAnalyzer,
                        "/home/burkhart/Software/Ogmios/datasets/Mechismo/TCGA_mech_output.tsv",
                        null,
                        "SKCM",
                        0.0d,
                        1.0d,
                        "/home/burkhart/Software/Ogmios/datasets/FIsInGene_031516_with_annotations.txt",
                        "/home/burkhart/Software/Ogmios/datasets/ReactionNetwork_070517.txt",
                        "/home/burkhart/Software/Ogmios/results/Mechismo/"
                );
            } else {
                System.out.println("An error occurred.");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
