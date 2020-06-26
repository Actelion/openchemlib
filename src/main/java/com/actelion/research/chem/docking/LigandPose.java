package com.actelion.research.chem.docking;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;

import com.actelion.research.calc.Matrix;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.BondRotationHelper;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.TorsionDB;
import com.actelion.research.chem.conf.TorsionDetail;
import com.actelion.research.chem.docking.InteractionTerm;
import com.actelion.research.chem.conf.torsionstrain.StatisticalTorsionPotential;
import com.actelion.research.chem.conf.torsionstrain.StatisticalTorsionTerm;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import com.actelion.research.chem.interactionstatistics.InteractionDistanceStatistics;
import com.actelion.research.chem.interactionstatistics.SplineFunction;
import com.actelion.research.chem.phesa.Evaluable;
import com.actelion.research.chem.potentialenergy.AngleConstraint;
import com.actelion.research.chem.potentialenergy.BondConstraint;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;
import com.actelion.research.chem.potentialenergy.TorsionConstraint;

import java.util.AbstractMap.SimpleEntry;

public class LigandPose implements Evaluable{
	
	private BondRotationHelper torsionHelper;
	private List<int[]> ligAtomPairs; //separated by more than 3 bonds for internal strain
	private double[] torsionValues;
	private double[] rotation; //quaternion
	private double[] translation;
	private double[] state; //coordinates
	private double[] gradient; //gradient with respect to coordinates
	private Conformer ligConf;
	private StereoMolecule mol;
	private List<PotentialEnergyTerm> ligStrain;
	private int[] ligAtomTypes;
	
	public LigandPose(Conformer ligConf) {
		this.ligConf = ligConf;
		init();
	}
	
	private void init() {
		InteractionDistanceStatistics.getInstance().initialize();
		StatisticalTorsionPotential.getInstance().initialize();
		mol = ligConf.getMolecule();
		state = new double[3*mol.getAllAtoms()];
		gradient = new double[3*mol.getAllAtoms()];
		ligStrain = new ArrayList<PotentialEnergyTerm>();
		List<Integer> ligAtomTypesList = new ArrayList<>();
		for(int a=0;a<mol.getAtoms();a++) {
			ligAtomTypesList.add(InteractionAtomTypeCalculator.getAtomType(mol, a));
		}
		ligAtomTypes = new int[ligAtomTypesList.size()];
		IntStream.range(0, ligAtomTypes.length).forEach(e -> ligAtomTypes[e] = ligAtomTypesList.get(e));
		Coordinates com = DockingUtils.getCOM(ligConf);
		ligConf.translate(-com.x, -com.y, -com.z); // center-of-mass is at the origin of the coordinate system
		DockingUtils.createInitialOrientation(ligConf); // principal moments of inertia are aligned to axis of coordinate system
		updateState();
		translation = new double[] {0.0,0.0,0.0};
		rotation = new double[] {1.0,0.0,0.0,0.0};
		torsionHelper = new BondRotationHelper(mol);
		ligAtomPairs = new ArrayList<int[]>();
		createStrainFunction();
		
		
	}
	
	private void createStrainFunction() {
		findLigAtomPairs();
		for(int[] pair : ligAtomPairs) {
			int at1 = pair[0];
			int at2 = pair[1];
			InteractionTerm term = InteractionTerm.create(ligConf,ligConf,at1,at2,ligAtomTypes,ligAtomTypes);
			if(term==null)
				continue;
			ligStrain.add(term);
		}
		//add torsions around rotatable bonds to term list and restrain the others to the current value, to prevent distortions
		for(int b=0;b<mol.getBonds();b++) {
			boolean isRotBond = false;
			for(int rotBondIndex=0;rotBondIndex<torsionHelper.getRotatableBonds().length;rotBondIndex++) {
				int rotBond = torsionHelper.getRotatableBonds()[rotBondIndex];
				if(b==rotBond) {
					isRotBond = true;
					int[] torsionAtoms = torsionHelper.getTorsionAtoms()[rotBondIndex];
					String torsionID = torsionHelper.getTorsionIDs()[rotBondIndex];
					StatisticalTorsionTerm term = StatisticalTorsionTerm.create(ligConf,torsionAtoms,torsionID);
					if(term==null)
						continue;
					ligStrain.add(term);
					break;
				}
			}
			if(!isRotBond) {
				int at1 = mol.getBondAtom(0, b);
				int at2 = mol.getBondAtom(0, 1);
				if(mol.getConnAtoms(at1)==1 || mol.getConnAtoms(at2)==1) //terminal bond
					continue;
				int at3 = mol.getConnAtom(at1, 0) == at2 ? mol.getConnAtom(at1, 1) : mol.getConnAtom(at1, 0);
				int at4 = mol.getConnAtom(at2, 0) == at1 ? mol.getConnAtom(at2, 1) : mol.getConnAtom(at2, 0);
				int[] torsionAtoms = new int[] {at1,at2,at3,at4};
				Coordinates c1 = ligConf.getCoordinates(at1);
				Coordinates c2 = ligConf.getCoordinates(at2);
				Coordinates c3 = ligConf.getCoordinates(at3);
				Coordinates c4 = ligConf.getCoordinates(at4);
				double dihedral = 360.0*Coordinates.getDihedral(c1, c2, c3, c4)/(2*Math.PI);
				TorsionConstraint constraint = new TorsionConstraint(ligConf,torsionAtoms,dihedral,5.0);
				ligStrain.add(constraint);
			}
			int at1 = mol.getBondAtom(0, b);
			int at2 = mol.getBondAtom(0, 1);
			Coordinates c1 = mol.getCoordinates(at1);
			Coordinates c2 = mol.getCoordinates(at2);
			double dist = c1.distance(c2);
			BondConstraint constraint = new BondConstraint(ligConf,new int[] {at1,at2},dist);
			ligStrain.add(constraint);
		}
		// create angle constraints
        for (int at1=0; at1<mol.getAtoms(); at1++) {
            if (mol.getConnAtoms(at1) > 1) {
                for (int i=0; i<mol.getConnAtoms(at1); i++) {
                    int at2 = mol.getConnAtom(at1, i);
                    for (int k=i+1; k<mol.getConnAtoms(at1); k++) {
                        int at3 = mol.getConnAtom(at1, k);
                        Coordinates c1 = mol.getCoordinates(at1);
                        Coordinates c2 = mol.getCoordinates(at2);
                        Coordinates c3 = mol.getCoordinates(at3);
                        Coordinates v1 = c2.subC(c1);
                        Coordinates v2 = c3.subC(c1);
                        v1.unit();
                        v2.unit();
                        double alpha = Math.acos(v1.dot(v2));
                        AngleConstraint constraint = new AngleConstraint(ligConf,new int[] {at2,at1,at3},alpha);
                        ligStrain.add(constraint);
                    }
                }
            }
        }

	}
	
	private void findLigAtomPairs() {
		Set<SimpleEntry<Integer,Integer>> invalidPairs  = new HashSet<SimpleEntry<Integer,Integer>>();
		StereoMolecule mol = ligConf.getMolecule();
		for(int a=0;a<mol.getAtoms();a++) {
			for(int i=0;i<mol.getConnAtoms(a);i++) {
				int aa = mol.getConnAtom(a, i);
				SimpleEntry<Integer,Integer> entry = a<aa ? new SimpleEntry<>(a,aa) : new SimpleEntry<>(a,aa);
				invalidPairs.add(entry);
				for(int j=0;j<mol.getConnAtoms(aa);j++) {
					int aaa = mol.getConnAtom(aa, j);
					entry = a<aaa ? new SimpleEntry<>(a,aaa) : new SimpleEntry<>(a,aaa);
					invalidPairs.add(entry);
					for(int k=0;k<mol.getConnAtoms(aaa);k++) {
						int aaaa = mol.getConnAtom(aaa, k);
						entry = a<aaaa ? new SimpleEntry<>(a,aaaa) : new SimpleEntry<>(a,aaaa);
						invalidPairs.add(entry);
					}
				}
			}
		}
		for(int i=0;i<mol.getAtoms();i++) {
			for(int j=i+1;j<mol.getAtoms();j++) {
				SimpleEntry<Integer,Integer> entry =new SimpleEntry<>(i,j);
				if(invalidPairs.contains(entry))
					continue;
				else 
					ligAtomPairs.add(new int[] {i,j});
					
			}
		}
	}
		
		
	//constrain bonds that are not rotatable, constrain bond lengths and angles
	public double getFGValue(double gradient[]) {
		double energy = 0.0;
		for(PotentialEnergyTerm term : ligStrain) {
			energy+=term.getFGValue(gradient);
		}
		return energy;
			
	}
		
	
	public void updateState() {
		for(int a=0;a<mol.getAllAtoms();a++) {
			Coordinates c = ligConf.getCoordinates(a);
			state[3*a] = c.x;
			state[3*a+1] = c.y;
			state[3*a+2] = c.z;
		}
	}
	

	@Override
	public void setState(double[] state){
		assert this.state.length==state.length;
		for(int i=0;i<state.length;i++) {
			this.state[i] = state[i];
		}
		for(int a=0;a<mol.getAllAtoms();a++) {
			Coordinates c = new Coordinates(state[3*a],state[3*a+1],state[3*a+2]);
			ligConf.setCoordinates(a, c);
		}
	}
	
	public double[] getState(double[] v){
		for(int i=0;i<this.state.length;i++) {
			v[i] = state[i];
			
		}
		return v;
	}
	
	public double[] getState() {
		return this.getState(new double[state.length]);
	}
		
}
	
	
	


