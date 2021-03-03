package com.actelion.research.chem.phesa.pharmacophore;

import java.util.Base64;
import java.util.Base64.Decoder;
import java.util.Base64.Encoder;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.PeriodicTable;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.phesa.EncodeFunctions;
import com.actelion.research.chem.phesa.Gaussian3D;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.pharmacophore.IPharmacophorePoint.Functionality;

/**
 * Differs from the PPGaussian by lacking a directionality vector and is not associated with a molecule
 * Can be used to represent interaction sites in protein binding sites
 * @author joel
 *
 */

public class SimplePPGaussian extends Gaussian3D {
	
	private IPharmacophorePoint.Functionality functionality;
	
	public SimplePPGaussian(String encodedGaussian) {
		decode(encodedGaussian);
	}
	
	public SimplePPGaussian(Coordinates center, IPharmacophorePoint.Functionality functionality) {
		atomicNo = 6;
		atomId = -1;
		this.center = center;
		this.functionality = functionality;
	}

	@Override
	public double calculateHeight() {
		return MolecularVolume.p;
	}

	@Override
	public double calculateWidth() {
		double vdwR = PeriodicTable.getElement(atomicNo).getVDWRadius();
		return MolecularVolume.alpha_pref/(vdwR*vdwR);
	}

	@Override
	public String encode() {
		Encoder encoder = Base64.getEncoder();
		StringBuilder molVolString = new StringBuilder();
		molVolString.append(Integer.toString(atomicNo));
		molVolString.append(" ");
		double[] coords = new double[] {center.x,center.y,center.z};
		molVolString.append(encoder.encodeToString(EncodeFunctions.doubleArrayToByteArray(coords)));
        molVolString.append(" ");
		molVolString.append(encoder.encodeToString(EncodeFunctions.doubleToByteArray(weight)));
        molVolString.append(" ");
        molVolString.append(functionality.getIndex());
		return molVolString.toString();
	}
	
	private void decode(String string64)  {
		Decoder decoder = Base64.getDecoder();
		String[] strings = string64.split(" ");
		if(strings.length==1) { // no pharmacophore information encoded
			return;
		}
		atomId = -1;
		atomicNo = Integer.decode(strings[0]);
		double[] coords = EncodeFunctions.byteArrayToDoubleArray(decoder.decode(strings[1].getBytes()));
		center = new Coordinates(coords[0],coords[1],coords[2]);
		weight = EncodeFunctions.byteArrayToDouble(decoder.decode(strings[2].getBytes()));
		int functionalityID = Integer.valueOf(strings[3]);
        for (IPharmacophorePoint.Functionality f : IPharmacophorePoint.Functionality.values()) {
            if (f.getIndex()==functionalityID) {
                functionality = f;
                break;
            }
        }

		
		alpha = calculateWidth(); //the width of the Gaussian depends on the atomic radius of the atom
		volume = calculateVolume();
		coeff = calculateHeight();
	}
	
	public Functionality getFunctionality() {
		return functionality;
	}

	
	

}
