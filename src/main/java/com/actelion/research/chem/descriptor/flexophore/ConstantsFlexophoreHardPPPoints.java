package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.phesa.pharmacophore.PharmacophoreCalculator;
import com.actelion.research.chem.phesa.pharmacophore.pp.IPharmacophorePoint;

public class ConstantsFlexophoreHardPPPoints {

    public static final String ATTR_DONOR = "d";
    public static final String ATTR_ACCEPTOR = "a";
    public static final String ATTR_NEGATIVE_CHARGE = "-";
    public static final String ATTR_POSITIVE_CHARGE = "+";
    public static final String ATTR_AROMATIC = "r";
    public static final String ATTR_ALIPHATIC = "l";


    public static String toStringPPPoints(int type) {

        String s = "";

        if(type==IPharmacophorePoint.Functionality.ACCEPTOR.getIndex()){
            s=ATTR_ACCEPTOR;
        } else if(type==IPharmacophorePoint.Functionality.DONOR.getIndex()){
            s=ATTR_DONOR;
        } else if(type==IPharmacophorePoint.Functionality.NEG_CHARGE.getIndex()){
            s=ATTR_NEGATIVE_CHARGE;
        } else if(type==IPharmacophorePoint.Functionality.POS_CHARGE.getIndex()){
            s=ATTR_POSITIVE_CHARGE;
        } else if(type==IPharmacophorePoint.Functionality.AROM_RING.getIndex()){
            s=ATTR_AROMATIC;
        } else if(type== PharmacophoreCalculator.LIPO_ID){
            s=ATTR_ALIPHATIC;
        } else {
            throw new RuntimeException("Unknown pharmacophore point type: " + type + "!");
        }

        return s;
    }
}
