/*
 * Created on Jun 23, 2010
 *
 */
package org.reactome.annotate;

import java.util.List;
import java.util.Set;

import javax.xml.bind.annotation.XmlRootElement;

/**
 * This class wrap a GeneSetAnnotation for a network module
 * @author wgm
 *
 */
@XmlRootElement
public class ModuleGeneSetAnnotation {
    private List<GeneSetAnnotation> annotations;
    private int module;
    private Set<String> ids;
    
    public ModuleGeneSetAnnotation() {
    }

    public List<GeneSetAnnotation> getAnnotations() {
        return annotations;
    }

    public void setAnnotations(List<GeneSetAnnotation> annotations) {
        this.annotations = annotations;
    }

    public int getModule() {
        return module;
    }

    public void setModule(int module) {
        this.module = module;
    }

    public Set<String> getIds() {
        return ids;
    }

    public void setIds(Set<String> ids) {
        this.ids = ids;
    }
    
}
