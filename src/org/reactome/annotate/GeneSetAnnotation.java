/*
 * Created on Jun 22, 2010
 *
 */
package org.reactome.annotate;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

/**
 * A data structure to store annotations from pathways, GO terms or other sources.
 * @author wgm
 *
 */
@XmlRootElement
public class GeneSetAnnotation {
    // Annotated terms (e.g. pathway or go terms)
    private String topic;
    private Integer numberInTopic;
    private Double ratioOfTopic;
    // ids hit by the topic
    private Integer hitNumber;
    private Double pValue;
    // Use String to indicate a case like <0.001
    private String fdr;
    // Hit ids
    private List<String> hitIds;
    
    public GeneSetAnnotation() {
    }
    
    public Integer getNumberInTopic() {
        return numberInTopic;
    }

    public void setNumberInTopic(Integer numberInTopic) {
        this.numberInTopic = numberInTopic;
    }

    public Double getRatioOfTopic() {
        return ratioOfTopic;
    }

    public void setRatioOfTopic(Double ratioOfTopic) {
        this.ratioOfTopic = ratioOfTopic;
    }

    public String getTopic() {
        return topic;
    }

    public void setTopic(String topic) {
        this.topic = topic;
    }

    public Integer getHitNumber() {
        return hitNumber;
    }

    public void setHitNumber(Integer hitNumber) {
        this.hitNumber = hitNumber;
    }

    public Double getPValue() {
        return pValue;
    }

    public void setPValue(Double pValue) {
        this.pValue = pValue;
    }

    public String getFdr() {
        return fdr;
    }

    public void setFdr(String fdr) {
        this.fdr = fdr;
    }

    public List<String> getHitIds() {
        return hitIds;
    }

    public void setHitIds(List<String> hitIds) {
        this.hitIds = hitIds;
    }
    
    public void addHitIds(String hitId) {
        if (hitIds == null)
            hitIds = new ArrayList<String>();
        hitIds.add(hitId);
    }
    
}
