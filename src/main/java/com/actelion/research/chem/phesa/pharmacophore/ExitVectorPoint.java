package com.actelion.research.chem.phesa.pharmacophore;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.interactionstatistics.InteractionSimilarityTable;

public class ExitVectorPoint implements IPharmacophorePoint {
	private int coreAtom;
	private int exitAtom;
	private Coordinates directionality;
	private Coordinates center;

	
	public ExitVectorPoint(StereoMolecule mol, int c, int e) {
		coreAtom = c;
		exitAtom = e;
		updateCoordinates(mol);
	}
	
	private ExitVectorPoint(String ppString, StereoMolecule mol) {
		decode(ppString,mol);
	}
	
	public ExitVectorPoint(ExitVectorPoint eP) {
		coreAtom = eP.coreAtom;
		exitAtom = eP.exitAtom;
		directionality = new Coordinates(eP.directionality);
		center = new Coordinates(eP.center);
	}
	
	public static ExitVectorPoint fromString(String ppString, StereoMolecule mol) {
		return new ExitVectorPoint(ppString,mol);
	}

	@Override
	public void updateCoordinates(StereoMolecule mol) {
		center = new Coordinates(mol.getAtomX(coreAtom),mol.getAtomY(coreAtom),mol.getAtomZ(coreAtom));
		directionality = mol.getCoordinates(exitAtom).subC(mol.getCoordinates(coreAtom));
		directionality.scale(1.0/directionality.getLength());
		
	}
	
	@Override
	public void updateCoordinates(Conformer conf) {
		center = new Coordinates(conf.getX(coreAtom),conf.getY(coreAtom),conf.getZ(coreAtom));
		directionality = conf.getCoordinates(exitAtom).subC(conf.getCoordinates(coreAtom));
		directionality.scale(1.0/directionality.getLength());
		
	}
	


	@Override
	public Coordinates getCenter() {
		return center;
	}

	@Override
	public Coordinates getDirectionality() {
		return directionality;
	}

	@Override
	public String encode() {
		StringBuilder molVolString = new StringBuilder();
		molVolString.append("e");
		molVolString.append(" ");
		molVolString.append(Integer.toString(coreAtom));
		molVolString.append(" ");
		molVolString.append(Integer.toString(exitAtom));
		return molVolString.toString();
	}
	

	private void decode(String ppString, StereoMolecule mol) {
		String[] strings = ppString.split(" ");
		coreAtom = Integer.decode(strings[1]);
		exitAtom = Integer.decode(strings[2]);
		updateCoordinates(mol);
	}

	@Override
	public double getSimilarity(IPharmacophorePoint pp) {
		if(pp instanceof ExitVectorPoint) {
			return 1.0;
		}
		return 0.0;
	}




	@Override
	public int getCenterID() {
		return coreAtom;
	}

	@Override
	public void setDirectionality(Coordinates directionality) {
		this.directionality = directionality;
		
	}
	
	@Override
	public void updateAtomIndeces(int[] map) {
		coreAtom = map[coreAtom];
		exitAtom = map[exitAtom];
		
	}

	@Override
	public IPharmacophorePoint copyPharmacophorePoint() {
		// TODO Auto-generated method stub
		return new ExitVectorPoint(this);
	}

	@Override
	public void getDirectionalityDerivativeCartesian(double[] grad, double[] v, Coordinates di, double sim) {
		grad[3*exitAtom] = sim*di.x/3.0;
		grad[3*exitAtom+1] = sim*di.y/3.0;
		grad[3*exitAtom+2] = sim*di.z/3.0;
		grad[3*coreAtom] = sim*-di.x/3.0;
		grad[3*coreAtom+1] = sim*-di.y/3.0;
		grad[3*coreAtom+2] = sim*-di.z/3.0;
		
	}
	
	@Override
	public int getFunctionalityIndex() {
		return IPharmacophorePoint.Functionality.EXIT_VECTOR.getIndex();
	}
	

	
	
		
	
}
