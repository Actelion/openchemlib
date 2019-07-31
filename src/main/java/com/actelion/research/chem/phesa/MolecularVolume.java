package com.actelion.research.chem.phesa;

import com.actelion.research.chem.AtomFunctionAnalyzer;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import java.util.ArrayList;
import com.actelion.research.util.EncoderFloatingPointNumbers;



/**
 * @version: 1.0, February 2018
 * Author: J. Wahl
 * class to approximate the volume of a molecule as a sum of atom-centered Gaussians, as introduced by Grant and Pickup, J. Phys. Chem. 1995, 99, 3503-3510
 * no higher order terms (atom-atom overlaps) are calculated to reduce computational costs
*/



public class MolecularVolume {
	static public final double p = 2.82842712475; // height of Atomic Gaussian, 2*sqrt(2), commonly used in the literature: Haque and Pande, DOI 10.1002/jcc.11307 
	static public final double alpha_pref = 2.41798793102; // taken from DOI 10.1002/jcc.11307
	static public final double pPP = 4.242640687125; // height of Pharmacophore-Gaussian, 1.5 times height of atomic gaussian
	static public final double alpha_prefPP = 3.11495261034551579; // chosen in a way that the Volume of a Pharmacophore Point coincides with Volume of a Carbon Atom
	
	

	private double volume;
	private Coordinates com;
	private ArrayList<AtomicGaussian> atomicGaussians;
	private ArrayList<PPGaussian> ppGaussians;
	private ArrayList<Coordinates> hydrogens;

	
	public MolecularVolume(ArrayList<AtomicGaussian> atomicGaussiansInp,ArrayList<PPGaussian> ppGaussiansInp, ArrayList<Coordinates> hydrogenCoords) {
		this.volume = 0.0;
		this.atomicGaussians = new ArrayList<AtomicGaussian>();
		for(AtomicGaussian ag : atomicGaussiansInp) {
			atomicGaussians.add(new AtomicGaussian(ag));
		}
		this.ppGaussians = new ArrayList<PPGaussian>();
		for(PPGaussian pg : ppGaussiansInp) {
			ppGaussians.add(new PPGaussian(pg));
		}
		this.hydrogens = new ArrayList<Coordinates>();
		for(Coordinates hydrogen : hydrogenCoords) {
			this.hydrogens.add(hydrogen);
			
		}

		this.calcCOM();

		
	}
	

	


	
	public MolecularVolume(StereoMolecule mol) {
		this.volume = 0.0;
		this.hydrogens = new ArrayList<Coordinates>();
		this.calc(mol);
		this.calcPPVolume(mol);
		this.calcCOM();

		
	}
	
	
	
	
	public MolecularVolume(MolecularVolume original) {
		this.volume = new Double(original.volume);
		this.atomicGaussians = new ArrayList<AtomicGaussian>();
		this.ppGaussians = new ArrayList<PPGaussian>();
		for(AtomicGaussian ag : original.getAtomicGaussians()) {
			this.atomicGaussians.add(new AtomicGaussian(ag));
		}
		for(PPGaussian pg : original.getPPGaussians()) {
			this.ppGaussians.add(new PPGaussian(pg));
		}

		this.hydrogens = new ArrayList<Coordinates>();
		for(Coordinates hydrogen : original.hydrogens) {
			this.hydrogens.add(hydrogen);
			
		}
		this.com = original.com;
		
	}
	

	

	
	/**
	 * calculates the molecular Volume for a StereoMolecule with 3D coordinates
	 * @param mol
	 */
	
	private void calc(StereoMolecule mol) { 
		this.atomicGaussians = new ArrayList<AtomicGaussian>();
		int nrOfAtoms = mol.getAllAtoms();
		for (int i=0;i<nrOfAtoms;i++) {
			
			if(mol.getAtomicNo(i)==1){ //hydrogens don't contribute to the molecular volume
				this.hydrogens.add(mol.getCoordinates(i));
				continue;
			}
			Coordinates coords = new Coordinates(mol.getCoordinates(i));
			AtomicGaussian atomicGaussian = new AtomicGaussian(i,mol.getAtomicNo(i),coords);
			this.atomicGaussians.add(atomicGaussian);
		}
	}
	
	/**
	 * calculates the pharmacophore points and corresponding volumes of a 3d molecule
	 * @param mol
	 */
	
	private void calcPPVolume(StereoMolecule mol) {
		this.ppGaussians = new ArrayList<PPGaussian>();


		for(int i=0;i<mol.getAllAtoms();i++) {
				if (mol.getAtomicNo(i)==1) {
					if(isDonorHydrogen(mol,i)) {
						int d = mol.getConnAtom(i,0);
						int interactionClass = InteractionAtomTypeCalculator.getAtomType(mol, d);
						if(interactionClass<0) {
							continue;
						}
						DonorPoint dp = new DonorPoint(mol,d,i,interactionClass);
						this.ppGaussians.add(new PPGaussian(6,dp));	
				}
			}
				else if (mol.getAtomicNo(i)==7 || mol.getAtomicNo(i)==8) {
					if(isAcceptor(mol,i)) {
						int neighbours = mol.getAllConnAtoms(i);
						ArrayList<Integer> neighbourList= new ArrayList<Integer>();
						for(int j=0;j<neighbours;j++) 
							neighbourList.add(mol.getConnAtom(i,j));
						
						int interactionClass = InteractionAtomTypeCalculator.getAtomType(mol, i);
						if(interactionClass<0) {
							continue;
						}
						if(mol.getAtomicNo(i)==8 && neighbours==1 && (mol.getConnBondOrder(i, 0)==2 || AtomFunctionAnalyzer.isAcidicOxygen(mol, i) )) {
							int aa1 = mol.getConnAtom(mol.getConnAtom(i,0),0);
							if(aa1==i) 
								aa1 = mol.getConnAtom(mol.getConnAtom(i,0),1);
							neighbourList.add(aa1);
							AcceptorPoint ap = new AcceptorPoint(mol,i,neighbourList,interactionClass,1);
							this.ppGaussians.add(new PPGaussian(6,ap));
							ap = new AcceptorPoint(mol,i,neighbourList,interactionClass,2);
							this.ppGaussians.add(new PPGaussian(6,ap));
						}
						else {
						AcceptorPoint ap = new AcceptorPoint(mol,i,neighbourList,interactionClass);
						this.ppGaussians.add(new PPGaussian(6,ap));	
						}
				}
			}

		
		
	}
	}


	
		
	/**
	 * calculates volume weighted center of mass of the molecular Volume
	 */
	
	private void calcCOM(){ 
		double volume = 0.0;
		double comX = 0.0;
		double comY = 0.0;
		double comZ = 0.0;
		for(AtomicGaussian atGauss : this.atomicGaussians){
			volume += atGauss.getVolume();
			comX += atGauss.getCenter().x*atGauss.getVolume();
			comY += atGauss.getCenter().y*atGauss.getVolume();
			comZ += atGauss.getCenter().z*atGauss.getVolume();
		}

		comX = comX/volume;
		comY = comY/volume;
		comZ = comZ/volume;
		this.volume = volume;
		this.com = new Coordinates(comX,comY,comZ);

	}
	
	private boolean isAcceptor(StereoMolecule mol, int a) {
		if (mol.getAtomCharge(a)<=0) { //charge is not positive -> no acceptor
			if (mol.isAromaticAtom(a)) { 
				if (mol.getAllConnAtoms(a)<3) {
					return true;
				}
				else {
					return false; //is in aromatic ring and has at least 3 bonded atoms -> no acceptor 
				}
			}
			else if (mol.getAtomicNo(a)==7){ // atom is not aromatic
				if (mol.isFlatNitrogen(a)) 
					return false;
				for(int i=0;i<mol.getAllConnAtoms(a);i++) {
					int aa = mol.getConnAtom(a, i);
					if (mol.getAtomicNo(aa)==6) {
						for(int j=0;j<mol.getAllConnAtoms(aa);j++) {
							int aaa = mol.getConnAtom(aa,j);
							if(a==aaa) continue;
							if (mol.getBondOrder(mol.getBond(aa,aaa))==2) { 
								if (mol.getAtomicNo(aaa)==7) return false; //amide structure
								if (mol.getAtomicNo(aaa)==8) return false;
								if (mol.getAtomicNo(aaa)==16) return false;
							}
						}
					}

				}
			}

		return true;
		}
		else return false;
	}
	
	private boolean isDonorHydrogen(StereoMolecule mol, int h) {
		int dh = mol.getConnAtom(h, 0);
		if (mol.getAtomCharge(dh)>=0 && (mol.getAtomicNo(dh)==7 || mol.getAtomicNo(dh)==8) ) { //charge is not positive -> no acceptor
			return true;
		}
		else return false;
	}

	public Coordinates getCOM() {
		return this.com;
	}

	public ArrayList<AtomicGaussian> getAtomicGaussians() {
		return this.atomicGaussians;
	}
	
	public ArrayList<PPGaussian> getPPGaussians() {
		return this.ppGaussians;
	}
	


	public ArrayList<Coordinates> getHydrogens() {
		return this.hydrogens;
	}
	
	public StereoMolecule getConformer(StereoMolecule mol) {
		int nrOfAtoms = mol.getAllAtoms();
		StereoMolecule conformer = new StereoMolecule(mol);
		int hydrogenCounter = 0;
		ArrayList<Coordinates> hydrogens = this.getHydrogens();
		for(int i=0;i<nrOfAtoms;i++) {
			if(mol.getAtomicNo(i)==1){
				conformer.getCoordinates(i).set(hydrogens.get(hydrogenCounter));
				hydrogenCounter+=1;
			} 
			//else {
			//	conformer.getCoordinates(i).set(molVol.getAtomicGaussians().get(i).getCenter());
			//}
		for(int j=0;j<this.getAtomicGaussians().size();j++) {
			int atomId =this.getAtomicGaussians().get(j).getAtomId();
			conformer.getCoordinates(atomId).set(this.getAtomicGaussians().get(j).getCenter());
		}
		}

		return conformer;
	}

	


	public String encodeFull() {
		StringBuilder molVolString = new StringBuilder();
		molVolString.append(Integer.toString(atomicGaussians.size()));
		molVolString.append("  ");
		for(AtomicGaussian ag : atomicGaussians) {
			molVolString.append(ag.encode());
			molVolString.append("  ");

		}

		molVolString.append(Integer.toString(ppGaussians.size()));
		molVolString.append("  ");
		for(PPGaussian pg : ppGaussians) {
			molVolString.append(pg.encode().trim());
			molVolString.append("  ");

		}
		

		double[] hydrogenCoords = new double[3*hydrogens.size()];
		for(int i=0;i<hydrogens.size();i++) {
			hydrogenCoords[3*i] = hydrogens.get(i).x;
			hydrogenCoords[3*i+1] = hydrogens.get(i).y;
			hydrogenCoords[3*i+2] = hydrogens.get(i).z;
		}

		molVolString.append(EncoderFloatingPointNumbers.encode(hydrogenCoords,13));

		
		return molVolString.toString();

	}
	
	public String encodeCoordsOnly() {
		StringBuilder molVolString = new StringBuilder();
		//molVolString.append(Integer.toString(atomicGaussians.size()));
		double[] coords = new double[3*atomicGaussians.size()];
		for(int i=0;i<atomicGaussians.size();i++) {
			coords[3*i]=atomicGaussians.get(i).getCenter().x;
			coords[3*i+1]=atomicGaussians.get(i).getCenter().y;
			coords[3*i+2]=atomicGaussians.get(i).getCenter().z;
		}
		molVolString.append(EncoderFloatingPointNumbers.encode(coords, 13));
		molVolString.append("  ");

		
		coords = new double[3*ppGaussians.size()];
		for(int i=0;i<ppGaussians.size();i++) {
			coords[3*i]=ppGaussians.get(i).getCenter().x;
			coords[3*i+1]=ppGaussians.get(i).getCenter().y;
			coords[3*i+2]=ppGaussians.get(i).getCenter().z;
		}
		
		molVolString.append(EncoderFloatingPointNumbers.encode(coords, 13));
		molVolString.append("  ");
		
		coords = new double[3*ppGaussians.size()];
		for(int i=0;i<ppGaussians.size();i++) {
			coords[3*i]=ppGaussians.get(i).getPharmacophorePoint().getDirectionality().x;
			coords[3*i+1]=ppGaussians.get(i).getPharmacophorePoint().getDirectionality().y;
			coords[3*i+2]=ppGaussians.get(i).getPharmacophorePoint().getDirectionality().z;
		}
		molVolString.append(EncoderFloatingPointNumbers.encode(coords, 13));
		molVolString.append("  ");
		


		double[] hydrogenCoords = new double[3*hydrogens.size()];
		for(int i=0;i<hydrogens.size();i++) {
			hydrogenCoords[3*i] = hydrogens.get(i).x;
			hydrogenCoords[3*i+1] = hydrogens.get(i).y;
			hydrogenCoords[3*i+2] = hydrogens.get(i).z;
		}

		molVolString.append(EncoderFloatingPointNumbers.encode(hydrogenCoords,13));


		return molVolString.toString();

	}

	public static MolecularVolume decodeCoordsOnly(String string, MolecularVolume reference)  {
		ArrayList<AtomicGaussian> referenceAtomicGaussians = reference.getAtomicGaussians(); 
		ArrayList<PPGaussian> referencePPGaussians = reference.getPPGaussians(); 
		
		String[] splitString = string.split("  ");
		double[] atomicGaussiansCoords = EncoderFloatingPointNumbers.decode(splitString[0]);
		double[] ppGaussiansCoords = EncoderFloatingPointNumbers.decode(splitString[1]);
		double[] ppGaussiansDirectionalities = EncoderFloatingPointNumbers.decode(splitString[2]);
		double[] hydrogensCoords = EncoderFloatingPointNumbers.decode(splitString[3]);
		
		ArrayList<AtomicGaussian> atomicGaussians = new ArrayList<AtomicGaussian>();
		ArrayList<PPGaussian> ppGaussians = new ArrayList<PPGaussian>();
		ArrayList<Coordinates> hydrogens = new ArrayList<Coordinates>();
		
		int nrOfAtomicGaussians = atomicGaussiansCoords.length/3;
		int nrOfHydrogens = hydrogensCoords.length/3;
		int nrOfPPGaussians = ppGaussiansCoords.length/3;
		
		for(int i=0;i<nrOfAtomicGaussians;i++) {
			Coordinates coords = new Coordinates(atomicGaussiansCoords[i*3],atomicGaussiansCoords[i*3+1],atomicGaussiansCoords[i*3+2]);
			AtomicGaussian at = new AtomicGaussian(referenceAtomicGaussians.get(i));
			at.setCenter(coords);
			atomicGaussians.add(at);
		}
		
		for(int i=0;i<nrOfPPGaussians;i++) {
			Coordinates coords = new Coordinates(ppGaussiansCoords[i*3],ppGaussiansCoords[i*3+1],ppGaussiansCoords[i*3+2]);
			PPGaussian pp = new PPGaussian(referencePPGaussians.get(i));
			pp.setCenter(new Coordinates(coords.x, coords.y, coords.z));
			Coordinates directionality = new Coordinates(ppGaussiansDirectionalities[i*3],ppGaussiansDirectionalities[i*3+1],ppGaussiansDirectionalities[i*3+2]);
			pp.getPharmacophorePoint().setDirectionality(directionality);
			ppGaussians.add(pp);
		}
		


		for(int i=0;i<nrOfHydrogens;i++) {
			hydrogens.add(new Coordinates(hydrogensCoords[i*3],hydrogensCoords[i*3+1],hydrogensCoords[i*3+2]));
		}



		return new MolecularVolume(atomicGaussians,ppGaussians,hydrogens);
	}
	
	
	public static MolecularVolume decodeFull(String string, StereoMolecule refMol)  {
		String[] splitString = string.split("  ");
		int nrOfAtomicGaussians = Integer.decode(splitString[0].trim());
		int firstIndex = 1;
		int lastIndex = 1+nrOfAtomicGaussians;
		ArrayList<AtomicGaussian> atomicGaussians = new ArrayList<AtomicGaussian>();
		ArrayList<PPGaussian> ppGaussians = new ArrayList<PPGaussian>();
		ArrayList<Coordinates> hydrogens = new ArrayList<Coordinates>();
		
		for(int i=firstIndex;i<lastIndex;i++) {
			atomicGaussians.add(AtomicGaussian.fromString(splitString[i].trim()));
		}
		
		int nrOfPPGaussians = Integer.decode(splitString[lastIndex]);
		firstIndex = lastIndex+1;
		lastIndex = firstIndex + nrOfPPGaussians;
		for(int i=firstIndex;i<lastIndex;i++) {
			ppGaussians.add(PPGaussian.fromString(splitString[i],refMol));
		}
		

		
		double[] coords = EncoderFloatingPointNumbers.decode(splitString[splitString.length-1]);
		int nrOfHydrogens = coords.length/3;
		for(int i=0;i<nrOfHydrogens;i++) {
			hydrogens.add(new Coordinates(coords[i*3],coords[i*3+1],coords[i*3+2]));
			
		}


		return new MolecularVolume(atomicGaussians,ppGaussians,hydrogens);
	}
}

