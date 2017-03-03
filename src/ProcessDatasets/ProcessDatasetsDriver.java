package ProcessDatasets;

import org.reactome.cancer.MAFFileLoader;
import org.reactome.cancer.driver.CancerDriverReactomeAnalyzer;
import org.reactome.cancer.driver.Interactome3dDriverAnalyzer;

import java.io.IOException;

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

            CancerDriverReactomeAnalyzer reactomeAnalyzer = new CancerDriverReactomeAnalyzer();
            reactomeAnalyzer.SetMySqlCredentials(un,pw);

            interactome3dDriverAnalyzer.checkAllHumanReactions(reactomeAnalyzer);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
