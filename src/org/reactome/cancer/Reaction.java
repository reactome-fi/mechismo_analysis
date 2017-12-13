package org.reactome.cancer;

import java.util.Objects;

public final class Reaction{
    private final Long reactionID;
    private final String reactionName;
    public Reaction(Long reactionID,
                    String reactionName){
        this.reactionID = reactionID;
        this.reactionName = reactionName;
    }

    public Long getReactionID() {
        return reactionID;
    }

    public String getReactionName() {
        return reactionName;
    }

    @Override
    public int hashCode(){
        return Objects.hash(this.reactionID,this.reactionName);
    }

    @Override
    public boolean equals(Object o){
        return o == null
        ? false
        : o instanceof Reaction && o.hashCode() == this.hashCode();
    }

    @Override
    public String toString(){
        return String.format("%d:%s",this.reactionID,this.reactionName).replace(
                ",", "~");
    }
}
