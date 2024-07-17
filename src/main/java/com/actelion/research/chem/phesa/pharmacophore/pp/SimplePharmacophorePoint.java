package com.actelion.research.chem.phesa.pharmacophore.pp;

import java.util.Base64;
import java.util.Base64.Decoder;
import java.util.Base64.Encoder;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.phesa.EncodeFunctions;

/**
 * doesn't possess directionality terms (only a dummy vector) and is used to represent receptor pharmacophores (not dependent on a molecular
 * conformation, as opposed to ligand-based pharmacophores)
 * @author wahljo1
 *
 */

public class SimplePharmacophorePoint implements IPharmacophorePoint {
	
	private Coordinates center;
	private int atomID;
	private IPharmacophorePoint.Functionality functionality;
	private static final Coordinates directionality = new Coordinates(1.0,0.0,0.0);
	
	
	public SimplePharmacophorePoint(int atomID, Coordinates center, IPharmacophorePoint.Functionality functionality) {
		this.atomID = atomID;
		this.center = center;
		this.functionality = functionality;
	}
	
	public SimplePharmacophorePoint(SimplePharmacophorePoint point) {
		this.atomID = point.atomID;
		this.center = new Coordinates(point.center);
		this.functionality = point.functionality;
	}
	
	private SimplePharmacophorePoint(String ppString, StereoMolecule mol) {
		decode(ppString);
	}

	
	public static SimplePharmacophorePoint fromString(String ppString, StereoMolecule mol) {
		return new SimplePharmacophorePoint(ppString,mol);
	}

	@Override
	public Coordinates getCenter() {
		return center;
	}
	
	

	@Override
	public void updateCoordinates(Coordinates[] coords) {
		center = new Coordinates(coords[atomID]);
		
	}

	@Override
	public Coordinates getDirectionality() {
		return directionality;
	}

	@Override
	public String encode() {
		Encoder encoder = Base64.getEncoder();
		StringBuilder molVolString = new StringBuilder();
		molVolString.append("s");
		molVolString.append(" ");
		molVolString.append(Integer.toString(atomID));
		molVolString.append(" ");
		double[] coords = new double[] {center.x,center.y,center.z};
		molVolString.append(encoder.encodeToString(EncodeFunctions.doubleArrayToByteArray(coords)));
        molVolString.append(" ");
        molVolString.append(functionality.getIndex());
		return molVolString.toString();
	}
	
	public static SimplePharmacophorePoint decode(String encoded) {
		Decoder decoder = Base64.getDecoder();
		String[] strings = encoded.split(" ");
		int id = Integer.decode(strings[1]);
		double[] coords = EncodeFunctions.byteArrayToDoubleArray(decoder.decode(strings[2].getBytes()));
		Coordinates c = new Coordinates(coords[0],coords[1],coords[2]);
		Functionality func = null;
		int functionalityID = Integer.valueOf(strings[3]);
        for (IPharmacophorePoint.Functionality f : IPharmacophorePoint.Functionality.values()) {
            if (f.getIndex()==functionalityID) {
                func = f;
                break;
            }
        }
        SimplePharmacophorePoint sPP = new SimplePharmacophorePoint(id,c,func);
        return sPP;
	}
	

	@Override
	public double getSimilarity(IPharmacophorePoint pp) {
		double sim = 0.0;
		if(pp.getFunctionalityIndex()==functionality.getIndex())
			sim = 1.0;
		return sim;
	}

	@Override
	public int getCenterID() {
		return atomID;
	}
	
	@Override
	public void setCenterID(int centerID) {
		atomID = centerID;
	}

	@Override
	public void setDirectionality(Coordinates directionality) {
		return;
		
	}

	@Override
	public void updateAtomIndices(int[] map) {
		atomID = map[atomID];
	}

	@Override
	public int[] getAtomIndices() {
		int [] a = {atomID};
		return a;
	}

	@Override
	public IPharmacophorePoint copyPharmacophorePoint() {
			return new SimplePharmacophorePoint(this);
	}

	@Override
	public void getDirectionalityDerivativeCartesian(double[] grad, double[] v, Coordinates di, double sim) {
		return;
		
	}
	
	@Override
	public Coordinates getRotatedDirectionality(double[][] rotMatrix,double scaleFactor) {
		Coordinates directMod = new Coordinates(directionality);
		return directMod;
	}

	@Override
	public int getFunctionalityIndex() {
		return functionality.getIndex();
	}
	
	@Override
	public double getVectorSimilarity(IPharmacophorePoint pp2,Coordinates directionalityMod) {
		return 1.0;
	}

}
