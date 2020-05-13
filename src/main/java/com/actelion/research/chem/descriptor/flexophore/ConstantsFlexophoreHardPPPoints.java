package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.phesa.pharmacophore.PharmacophoreCalculator;

public class ConstantsFlexophoreHardPPPoints {

    public static final String ATTR_DONOR = "d";
    public static final String ATTR_ACCEPTOR = "a";
    public static final String ATTR_NEGATIVE_CHARGE = "-";
    public static final String ATTR_POSITIVE_CHARGE = "+";
    public static final String ATTR_AROMATIC = "r";
    public static final String ATTR_LIPO = "l";


    public static String toStringPPPoints(int type) {

        String s = "";

        switch (type){
            case PharmacophoreCalculator.ACCEPTOR_ID:
                s=ATTR_ACCEPTOR;
                break;
            case PharmacophoreCalculator.DONOR_ID:
                s=ATTR_DONOR;
                break;
            case PharmacophoreCalculator.CHARGE_NEG_ID:
                s=ATTR_NEGATIVE_CHARGE;
                break;
            case PharmacophoreCalculator.CHARGE_POS_ID:
                s=ATTR_POSITIVE_CHARGE;
                break;
            case PharmacophoreCalculator.AROM_ID:
                s=ATTR_AROMATIC;
                break;
            case PharmacophoreCalculator.LIPO_ID:
                s=ATTR_LIPO;
                break;
            default:
                throw new RuntimeException("Unknown pharmacophore point type: " + type + "!");
        }
        return s;
    }
}
