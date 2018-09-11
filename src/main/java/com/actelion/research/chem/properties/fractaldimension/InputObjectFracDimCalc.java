package com.actelion.research.chem.properties.fractaldimension;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.datamodel.IdentifiedObject;

/**
 * InputObjectFracDimCalc
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 28.08.18.
 */
class InputObjectFracDimCalc extends IdentifiedObject<StereoMolecule> {

    private String smiles;

    public InputObjectFracDimCalc(StereoMolecule data, long id, String smiles) {
        super(data, id);
        this.smiles = smiles;
    }

    public InputObjectFracDimCalc(InputObjectFracDimCalc inputObjectFracDimCalc) {
        super(inputObjectFracDimCalc.getData(), inputObjectFracDimCalc.getId());
        this.smiles = inputObjectFracDimCalc.smiles;
    }

    public String getSmiles() {
        return smiles;
    }
}
