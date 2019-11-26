package com.actelion.research.chem.phesa;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.phesa.pharmacophore.IPharmacophorePoint;
import com.actelion.research.chem.phesa.pharmacophore.IonizableGroupDetector;
import com.actelion.research.chem.phesa.pharmacophore.PPGaussian;
import com.actelion.research.chem.phesa.pharmacophore.PharmacophoreCalculator;

import java.util.ArrayList;
import java.util.List;

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

	
	

	private double volume;
	private Coordinates com;
	private ArrayList<AtomicGaussian> atomicGaussians;
	private ArrayList<PPGaussian> ppGaussians;
	private ArrayList<ExclusionGaussian> exclusionGaussians;
	private ArrayList<Coordinates> hydrogens;

	
	public MolecularVolume(ArrayList<AtomicGaussian> atomicGaussiansInp,ArrayList<PPGaussian> ppGaussiansInp, ArrayList<ExclusionGaussian> exclusionGaussians, ArrayList<Coordinates> hydrogenCoords) {
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
		
		this.exclusionGaussians = new ArrayList<ExclusionGaussian>();
		for(ExclusionGaussian eg : exclusionGaussians) {
			this.exclusionGaussians.add(new ExclusionGaussian(eg));
		}
		

		this.calcCOM();

		
	}
	

	public void updateCOM() {
		this.calcCOM();
	}
	
	private void updateAtomIndeces(List<? extends Gaussian3D> gaussians,int[] map) {
		for(Gaussian3D gaussian:gaussians)
			gaussian.updateAtomIndeces(map);
	}
	
	public void updateAtomIndeces(int[] map) {
		updateAtomIndeces(ppGaussians,map);
		updateAtomIndeces(atomicGaussians,map);
		updateAtomIndeces(exclusionGaussians,map);
	}
	


	
	public MolecularVolume(StereoMolecule mol) {
		this.volume = 0.0;
		this.hydrogens = new ArrayList<Coordinates>();
		this.exclusionGaussians = new ArrayList<ExclusionGaussian>();
		this.calc(mol);
		this.calcPPVolume(mol);
		this.calcCOM();

	}
	
	public MolecularVolume(MolecularVolume original, Conformer conf) {
		this(original);
		update(conf);

	}
	
	
	
	
	
	public MolecularVolume(MolecularVolume original) {
		this.volume = new Double(original.volume);
		this.atomicGaussians = new ArrayList<AtomicGaussian>();
		this.ppGaussians = new ArrayList<PPGaussian>();
		this.exclusionGaussians = new ArrayList<ExclusionGaussian>();
		for(AtomicGaussian ag : original.getAtomicGaussians()) {
			this.atomicGaussians.add(new AtomicGaussian(ag));
		}
		for(PPGaussian pg : original.getPPGaussians()) {
			this.ppGaussians.add(new PPGaussian(pg));
		}

		this.hydrogens = new ArrayList<Coordinates>();
		for(Coordinates hydrogen : original.hydrogens) {
			this.hydrogens.add(new Coordinates(hydrogen));
			
		}
		
		for(ExclusionGaussian eg : original.exclusionGaussians) {
			this.exclusionGaussians.add(new ExclusionGaussian(eg));
		}
		
		this.com = new Coordinates(original.com);
		
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
				this.hydrogens.add(new Coordinates(mol.getCoordinates(i)));
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
		ppGaussians = new ArrayList<PPGaussian>();
		List<IPharmacophorePoint> ppPoints = new ArrayList<IPharmacophorePoint>();
		IonizableGroupDetector detector = new IonizableGroupDetector(mol);
		ppPoints.addAll(detector.detect());
		ppPoints.addAll(PharmacophoreCalculator.getPharmacophorePoints(mol));
		for(IPharmacophorePoint ppPoint : ppPoints )
			ppGaussians.add(new PPGaussian(6,ppPoint));
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
	
	
	public Coordinates getCOM() {
		this.calcCOM();
		return this.com;
	}

	public ArrayList<AtomicGaussian> getAtomicGaussians() {
		return this.atomicGaussians;
	}
	
	public ArrayList<PPGaussian> getPPGaussians() {
		return this.ppGaussians;
	}
	
	public ArrayList<ExclusionGaussian> getExclusionGaussians() {
		return this.exclusionGaussians;
	}
	


	public ArrayList<Coordinates> getHydrogens() {
		return this.hydrogens;
	}
	
	private void updateHydrogens(StereoMolecule mol) {
		int h = 0;
		for(int i = mol.getAtoms();i<mol.getAllAtoms();i++) {
			hydrogens.get(h).set(new Coordinates(mol.getCoordinates(i)));
			h++;
		}
			
	}
	
	private void updateHydrogens(Conformer conf) {
		int h = 0;
		for(int i = conf.getMolecule().getAtoms();i<conf.getMolecule().getAllAtoms();i++) {
			hydrogens.get(h).set(new Coordinates(conf.getCoordinates(i)));
			h++;
		}
			
	}
	

	
	public void update(StereoMolecule mol) {
		updateCoordinates(getAtomicGaussians(),mol);
		updateCoordinates(getPPGaussians(),mol);
		updateCoordinates(getExclusionGaussians(),mol);
		updateHydrogens(mol);
	}
	
	public void update(Conformer conf) {
		updateCoordinates(getAtomicGaussians(),conf);
		updateCoordinates(getPPGaussians(),conf);
		updateCoordinates(getExclusionGaussians(),conf);
		updateHydrogens(conf);
	}
	
	private void updateCoordinates(ArrayList<? extends Gaussian3D> gaussians, StereoMolecule mol) {
		for(Gaussian3D gaussian : gaussians) {
			gaussian.updateCoordinates(mol);
		}
		
	}
	
	private void updateCoordinates(ArrayList<? extends Gaussian3D> gaussians, Conformer conf) {
		for(Gaussian3D gaussian : gaussians) {
			gaussian.updateCoordinates(conf);
		}
		
	}
	
	public void translateToCOM(Coordinates com) {


		for (AtomicGaussian ag : getAtomicGaussians()){
			ag.getCenter().sub(com);  //translate atomicGaussians. Moves center of mass to the origin.
		}

		
		for (PPGaussian pg : getPPGaussians()){
			pg.getCenter().sub(com);  //translate atomicGaussians. Moves center of mass to the origin.
		}
		
		for (ExclusionGaussian eg : getExclusionGaussians()){
			eg.getCenter().sub(com);  //translate atomicGaussians. Moves center of mass to the origin.
		}
		
		for (Coordinates hydrogen : getHydrogens()){
			hydrogen.sub(com);  //translate atomicGaussians. Moves center of mass to the origin.
		}
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
		
		molVolString.append(Integer.toString(exclusionGaussians.size()));
		
		molVolString.append("  ");
		
		for(ExclusionGaussian eg : exclusionGaussians) {
			molVolString.append(eg.encode());
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
		
		coords = new double[3*exclusionGaussians.size()];
		for(int i=0;i<exclusionGaussians.size();i++) {
			coords[3*i]=exclusionGaussians.get(i).getReferenceVector().x;
			coords[3*i+1]=exclusionGaussians.get(i).getReferenceVector().y;
			coords[3*i+2]=exclusionGaussians.get(i).getReferenceVector().z;
		}
		molVolString.append(EncoderFloatingPointNumbers.encode(coords, 13));
		molVolString.append("  ");
		
		coords = new double[3*exclusionGaussians.size()];
		for(int i=0;i<exclusionGaussians.size();i++) {
			coords[3*i]=exclusionGaussians.get(i).getShiftVector().x;
			coords[3*i+1]=exclusionGaussians.get(i).getShiftVector().y;
			coords[3*i+2]=exclusionGaussians.get(i).getShiftVector().z;
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
		ArrayList<ExclusionGaussian> referenceExclusionGaussians = reference.getExclusionGaussians();
		
		String[] splitString = string.split("  ");
		double[] atomicGaussiansCoords = EncoderFloatingPointNumbers.decode(splitString[0]);
		double[] ppGaussiansCoords = EncoderFloatingPointNumbers.decode(splitString[1]);
		double[] ppGaussiansDirectionalities = EncoderFloatingPointNumbers.decode(splitString[2]);
		double[] exclusionGaussiansRefCoords = EncoderFloatingPointNumbers.decode(splitString[3]);
		double[] exclusionGaussiansShiftCoords = EncoderFloatingPointNumbers.decode(splitString[4]);
		double[] hydrogensCoords = EncoderFloatingPointNumbers.decode(splitString[5]);
		
		ArrayList<AtomicGaussian> atomicGaussians = new ArrayList<AtomicGaussian>();
		ArrayList<PPGaussian> ppGaussians = new ArrayList<PPGaussian>();
		ArrayList<ExclusionGaussian> exclusionGaussians = new ArrayList<ExclusionGaussian>();
		ArrayList<Coordinates> hydrogens = new ArrayList<Coordinates>();
		
		int nrOfAtomicGaussians = atomicGaussiansCoords.length/3;
		int nrOfHydrogens = hydrogensCoords.length/3;
		int nrOfPPGaussians = ppGaussiansCoords.length/3;
		int nrOfExclusionGaussians = exclusionGaussiansRefCoords.length/3;
		
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
		
		for(int i=0;i<nrOfExclusionGaussians;i++) {
			Coordinates coords = new Coordinates(exclusionGaussiansRefCoords[i*3],exclusionGaussiansRefCoords[i*3+1],exclusionGaussiansRefCoords[i*3+2]);
			ExclusionGaussian eg = new ExclusionGaussian(referenceExclusionGaussians.get(i));
			eg.setReferenceVector(new Coordinates(coords.x, coords.y, coords.z));
			Coordinates shift = new Coordinates(exclusionGaussiansShiftCoords[i*3],exclusionGaussiansShiftCoords[i*3+1],exclusionGaussiansShiftCoords[i*3+2]);
			eg.setShiftVector(shift);
			exclusionGaussians.add(eg);
		}
		


		for(int i=0;i<nrOfHydrogens;i++) {
			hydrogens.add(new Coordinates(hydrogensCoords[i*3],hydrogensCoords[i*3+1],hydrogensCoords[i*3+2]));
		}



		return new MolecularVolume(atomicGaussians,ppGaussians, exclusionGaussians,hydrogens);
	}
	
	
	public static MolecularVolume decodeFull(String string, StereoMolecule refMol)  {
		String[] splitString = string.split("  ");
		int nrOfAtomicGaussians = Integer.decode(splitString[0].trim());
		int firstIndex = 1;
		int lastIndex = 1+nrOfAtomicGaussians;
		ArrayList<AtomicGaussian> atomicGaussians = new ArrayList<AtomicGaussian>();
		ArrayList<PPGaussian> ppGaussians = new ArrayList<PPGaussian>();
		ArrayList<ExclusionGaussian> exclusionGaussians = new ArrayList<ExclusionGaussian>();
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
		
		int nrOfExclusionGaussians = Integer.decode(splitString[lastIndex]);
		firstIndex = lastIndex+1;
		lastIndex = firstIndex + nrOfExclusionGaussians;
		for(int i=firstIndex;i<lastIndex;i++) {
			exclusionGaussians.add(ExclusionGaussian.fromString(splitString[i]));
		}
		

		
		double[] coords = EncoderFloatingPointNumbers.decode(splitString[splitString.length-1]);
		int nrOfHydrogens = coords.length/3;
		for(int i=0;i<nrOfHydrogens;i++) {
			hydrogens.add(new Coordinates(coords[i*3],coords[i*3+1],coords[i*3+2]));
			
		}

		return new MolecularVolume(atomicGaussians,ppGaussians,exclusionGaussians,hydrogens);
	}
}

