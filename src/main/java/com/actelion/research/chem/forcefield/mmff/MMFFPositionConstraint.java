package com.actelion.research.chem.forcefield.mmff;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.DoubleFormat;

public class MMFFPositionConstraint implements EnergyTerm {
	private final StereoMolecule mol;
	private final double[] refPos;
	private final boolean[] constrained;
	private final double k;
	private final double d;
	
	
	public MMFFPositionConstraint (StereoMolecule mol, double k, double d) {
		this.mol = mol;
		refPos = new double[3*mol.getAllAtoms()];
		constrained = new boolean[mol.getAllAtoms()];
		for(int a=0;a<mol.getAllAtoms();a++) {
			if(mol.getAtomicNo(a)==1)
				constrained[a]=false;
			else 
				constrained[a]=true;
			refPos[3*a] = mol.getAtomX(a);
			refPos[3*a+1] = mol.getAtomY(a);
			refPos[3*a+2] = mol.getAtomZ(a);
		}
		this.k = k;
		this.d = d;
	}

	@Override
	public double getEnergy(double[] pos) {
		return getEnergy(pos, null, null, false);
	}

	@Override
	public double getEnergy(double[] pos, StringBuilder detail, String detailID, boolean skipHydrogen) {
		double energy = 0.0;
		for(int a=0;a<pos.length;a+=3) {
			int atomId = a/3;
			if(!constrained[atomId] || (skipHydrogen && mol.getAtomicNo(atomId) == 1))
				continue;
			double dx = pos[a]-refPos[a];
			double dy = pos[a+1]-refPos[a+1];
			double dz = pos[a+2]-refPos[a+2];
			double dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
			double prefactor = 0.0;
			if(dist>d)
				prefactor = dist-d;
			else 
				prefactor = 0.0;
			energy+=0.5*k*prefactor*prefactor;
		}

		if (detail != null)
			detail.append(detailID+"\tposConstraint\tNaN\tNaN\tall\t"+DoubleFormat.toString(energy)+"\n");

		return energy;
	}

	@Override
	public void getGradient(double[] pos, double[] grad) {
		for(int a=0;a<pos.length;a+=3) {
			int atomId = a/3;
			if(!constrained[atomId])
				continue;
			double dx = pos[a]-refPos[a];
			double dy = pos[a+1]-refPos[a+1];
			double dz = pos[a+2]-refPos[a+2];
			double dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
			double prefactor = 0.0;
			if(dist>d)
				prefactor = dist-d;
			else 
				prefactor = 0.0;
			grad[a] += prefactor*dx/Math.max(dist, 1E-08);
			grad[a+1] += prefactor*dy/Math.max(dist, 1E-08);
			grad[a+2] += prefactor*dz/Math.max(dist, 1E-08);
		}
	}
}
