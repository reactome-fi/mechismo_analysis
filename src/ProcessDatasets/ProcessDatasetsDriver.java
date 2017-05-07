package ProcessDatasets;

import org.reactome.cancer.driver.CancerDriverReactomeAnalyzer;
import org.reactome.cancer.driver.Interactome3dDriverAnalyzer;

/**
 * Created by burkhart on 3/1/17.
 */
public class ProcessDatasetsDriver {
    public ProcessDatasetsDriver(String[] args){

    }

    public void run(){
        Interactome3dDriverAnalyzer interactome3dDriverAnalyzer = new Interactome3dDriverAnalyzer();
        try {
            System.out.print("MySQL Username: ");
            String un = new java.util.Scanner(System.in).next();

            System.out.print("MySQL Password: ");
            String pw = new java.util.Scanner(System.in).next();

            CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer = new CancerDriverReactomeAnalyzer();
            cancerDriverReactomeAnalyzer.setMySqlCredentials(un,pw);

            System.out.print("Interface Enrichment 0\n" +
                    "Mechismo/Reactome Overlay 1\n" +
                    "Mechismo/Reactome Interface Enrichment 2\n" +
                    "Cancel <enter>\n");
            String ex = new java.util.Scanner(System.in).next();

            if(Integer.parseInt(ex) == 0){
                interactome3dDriverAnalyzer.findInteractionsWithMutatedInterfaces(cancerDriverReactomeAnalyzer,
                        "datasets/firehose_all/",
                        // The MAF filenames look like:
                        // TCGA-AY-4070-01.hg19.oncotator.hugo_entrez_remapped.maf.txt
                        // TCGA-AY-4071-01.hg19.oncotator.hugo_entrez_remapped.maf.txt
                        "^.+\\.maf\\.txt$",
                        "datasets/interactome3d/2016_06/prebuilt/representative/",
                        "results/interactions_with_mutated_interfaces.csv");
            }else if(Integer.parseInt(ex) == 1){
                interactome3dDriverAnalyzer.findMechismoInteractionsInReactome(cancerDriverReactomeAnalyzer,
                        "datasets/Mechismo/COSMICv74_somatic_noSNPs_GWS_mechismo_output.tsv",
                        "^[^\\s]+\\s+(?<prot1>[^\\s]+)\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+(?<mech>[^\\s]+)\\s+[0-9a-zA-Z-]+\\s+(?<prot2>[0-9a-zA-Z-]+).*$",
                        "datasets/Mechismo/cancer_types/large_intestine+colon+carcinoma+adenocarcinoma/sample_id.txt",
                        "^(?<patientID>[^\\s]+)\\s+[^\\s]+\\s+(?<prot1>[^\\s]+)\\s+[^\\s]+\\s+[^\\s\\[\\]]+\\s+(?<prot2>[^\\s]+)\\s+(?<mech>[^\\s]+).*",
                        "results/interactions_with_mechismo_enrichment.csv");
            }else if(Integer.parseInt(ex) == 2){
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
                        "datasets/interactome3d/2016_06/prebuilt/representative/",
                        "results/interactions_with_mutated_interfaces.csv",
                        "results/interactions_with_mechismo_enrichment_and_mutated_interfaces.csv");
            } else{
                System.out.println("Execution Cancelled.");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
