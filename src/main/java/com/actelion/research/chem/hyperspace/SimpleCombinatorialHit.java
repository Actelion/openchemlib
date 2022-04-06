package com.actelion.research.chem.hyperspace;

import java.io.Serializable;

public class SimpleCombinatorialHit implements Serializable {
    private static final long serialVersionUID = 1970031554026082875L;
    public final String rxnId;
    public final SimpleSynthon[][] synthons;
    public SimpleCombinatorialHit(String rxnId, SimpleSynthon[][] synthons) {
        this.rxnId = rxnId;
        this.synthons = synthons;
    }
}