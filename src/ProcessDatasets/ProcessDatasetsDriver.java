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

            interactome3dDriverAnalyzer.findInteractionsWithMutatedInterfaces(cancerDriverReactomeAnalyzer,
                    "datasets/gdac.broadinstitute.org_COAD.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/",
                    // The MAF filenames look like:
                    // TCGA-AY-4070-01.hg19.oncotator.hugo_entrez_remapped.maf.txt
                    // TCGA-AY-4071-01.hg19.oncotator.hugo_entrez_remapped.maf.txt
                    "^.+\\.maf\\.txt$",
                    "datasets/interactome3d/2016_06/prebuilt/representative/",
                    "results/interactions_with_mutated_interfaces.csv");


        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
