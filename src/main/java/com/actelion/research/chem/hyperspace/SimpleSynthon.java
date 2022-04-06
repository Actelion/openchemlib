package com.actelion.research.chem.hyperspace;

import java.io.Serializable;

public class SimpleSynthon implements Serializable {
    private static final long serialVersionUID = 4625857679359578042L;
    public final String idcode;
    public final String synthonId;
    public final String rxnId;
    public final String synthonSet;
    public SimpleSynthon(String idcode, String id, String rxn_id, String synthonSet) {
        this.idcode     = idcode;
        this.synthonId  = id;
        this.rxnId      = rxn_id;
        this.synthonSet = synthonSet;
    }
}

