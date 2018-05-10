package org.reactome.mechismows;

import java.util.List;

import org.junit.Test;
import org.reactome.mechismo.config.AppConfig;
import org.reactome.mechismo.model.CancerType;
import org.reactome.mechismo.model.Interaction;
import org.reactome.mechismo.service.InteractionService;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

public class MechismowsReader {
    
    public MechismowsReader() {
        
    }
    
    @Test
    public void testLoadInteractions() throws Exception {
        AnnotationConfigApplicationContext context = getContext();
        List<Interaction> interactions = loadInteractions(context);
        System.out.println("Total interactions: " + interactions.size());
        context.close();
    }
    
    public AnnotationConfigApplicationContext getContext() {
        return new AnnotationConfigApplicationContext(AppConfig.class);
    }
    
    public List<Interaction> loadInteractions(AnnotationConfigApplicationContext context) {
        InteractionService service = context.getBean(InteractionService.class);
        List<Interaction> interactions = service.list(Interaction.class);
        return interactions;
    }
    
    public List<CancerType> loadCancerTypes(AnnotationConfigApplicationContext context) {
        InteractionService service = context.getBean(InteractionService.class);
        return service.list(CancerType.class);
    }

}
