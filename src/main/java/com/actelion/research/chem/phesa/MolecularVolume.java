package com.actelion.research.chem.phesa;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.phesa.pharmacophore.ExitVectorPoint;
import com.actelion.research.chem.phesa.pharmacophore.IPharmacophorePoint;
import com.actelion.research.chem.phesa.pharmacophore.IonizableGroupDetector;
import com.actelion.research.chem.phesa.pharmacophore.PPGaussian;
import com.actelion.research.chem.phesa.pharmacophore.PharmacophoreCalculator;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

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
	private ArrayList<VolumeGaussian> volumeGaussians; //exclusion and inclusion spheres
	private ArrayList<Coordinates> hydrogens;

	
	public MolecularVolume(List<AtomicGaussian> atomicGaussiansInp,List<PPGaussian> ppGaussiansInp,List<VolumeGaussian> volGaussians, 
			List<Coordinates> hydrogenCoords) {
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
		
		this.volumeGaussians = new ArrayList<VolumeGaussian>();
		for(VolumeGaussian eg : volGaussians) {
			this.volumeGaussians.add(new VolumeGaussian(eg));
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
		updateAtomIndeces(volumeGaussians,map);
	}
	


	
	public MolecularVolume(StereoMolecule mol) {
		this.volume = 0.0;
		this.hydrogens = new ArrayList<Coordinates>();
		this.volumeGaussians = new ArrayList<VolumeGaussian>();
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
		this.volumeGaussians = new ArrayList<VolumeGaussian>();
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
		
		for(VolumeGaussian eg : original.volumeGaussians) {
			this.volumeGaussians.add(new VolumeGaussian(eg));
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
			else if(mol.getAtomicNo(i)==0) {
				Coordinates coords = new Coordinates(mol.getCoordinates(i));
				AtomicGaussian atomicGaussian = new AtomicGaussian(i,6,coords);
				atomicGaussian.setWeight(0.0);
				this.atomicGaussians.add(atomicGaussian);
			}
			else {
				Coordinates coords = new Coordinates(mol.getCoordinates(i));
				AtomicGaussian atomicGaussian = new AtomicGaussian(i,mol.getAtomicNo(i),coords);
				this.atomicGaussians.add(atomicGaussian);
			}
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
		for(VolumeGaussian volGauss : this.volumeGaussians){
			volume += volGauss.getRole()*volGauss.getVolume();
			comX += volGauss.getRole()*volGauss.getCenter().x*volGauss.getVolume();
			comY += volGauss.getRole()*volGauss.getCenter().y*volGauss.getVolume();
			comZ += volGauss.getRole()*volGauss.getCenter().z*volGauss.getVolume();
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
	
	public ArrayList<VolumeGaussian> getVolumeGaussians() {
		return this.volumeGaussians;
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
		updateCoordinates(getVolumeGaussians(),mol);
		updateHydrogens(mol);
	}
	
	public void update(Conformer conf) {
		updateCoordinates(getAtomicGaussians(),conf);
		updateCoordinates(getPPGaussians(),conf);
		updateCoordinates(getVolumeGaussians(),conf);
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
		
		for (VolumeGaussian vg : getVolumeGaussians()){
			vg.translateRef(com.scaleC(-1.0));  //translate atomicGaussians. Moves center of mass to the origin.
		}
		

		
		for (Coordinates hydrogen : getHydrogens()){
			hydrogen.sub(com);  //translate atomicGaussians. Moves center of mass to the origin.
		}
		calcCOM();
	}
	
	public List<PPGaussian> getExitVectorGaussians() {
		return ppGaussians.stream().filter(e -> (e.getPharmacophorePoint() instanceof ExitVectorPoint)).collect(Collectors.toList());
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
		
		molVolString.append(Integer.toString(volumeGaussians.size()));
		
		molVolString.append("  ");
		
		for(VolumeGaussian vg : volumeGaussians) {
			molVolString.append(vg.encode());
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
		
		coords = new double[3*volumeGaussians.size()];
		for(int i=0;i<volumeGaussians.size();i++) {
			coords[3*i]=volumeGaussians.get(i).getReferenceVector().x;
			coords[3*i+1]=volumeGaussians.get(i).getReferenceVector().y;
			coords[3*i+2]=volumeGaussians.get(i).getReferenceVector().z;
		}
		molVolString.append(EncoderFloatingPointNumbers.encode(coords, 13));
		molVolString.append("  ");
		
		coords = new double[3*volumeGaussians.size()];
		for(int i=0;i<volumeGaussians.size();i++) {
			coords[3*i]=volumeGaussians.get(i).getShiftVector().x;
			coords[3*i+1]=volumeGaussians.get(i).getShiftVector().y;
			coords[3*i+2]=volumeGaussians.get(i).getShiftVector().z;
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
		ArrayList<VolumeGaussian> referenceVolGaussians = reference.getVolumeGaussians();
		
		
		String[] splitString = string.split("  ");
		double[] atomicGaussiansCoords = EncoderFloatingPointNumbers.decode(splitString[0]);
		double[] ppGaussiansCoords = EncoderFloatingPointNumbers.decode(splitString[1]);
		double[] ppGaussiansDirectionalities = EncoderFloatingPointNumbers.decode(splitString[2]);
		double[] volumeGaussiansRefCoords = EncoderFloatingPointNumbers.decode(splitString[3]);
		double[] volumeGaussiansShiftCoords = EncoderFloatingPointNumbers.decode(splitString[4]);
		double[] hydrogensCoords = EncoderFloatingPointNumbers.decode(splitString[5]);
		
		ArrayList<AtomicGaussian> atomicGaussians = new ArrayList<AtomicGaussian>();
		ArrayList<PPGaussian> ppGaussians = new ArrayList<PPGaussian>();
		ArrayList<VolumeGaussian> volumeGaussians = new ArrayList<VolumeGaussian>();
		ArrayList<Coordinates> hydrogens = new ArrayList<Coordinates>();
		
		int nrOfAtomicGaussians = atomicGaussiansCoords.length/3;
		int nrOfHydrogens = hydrogensCoords.length/3;
		int nrOfPPGaussians = ppGaussiansCoords.length/3;
		int nrOfVolumeGaussians = volumeGaussiansRefCoords.length/3;
		
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
		
		for(int i=0;i<nrOfVolumeGaussians;i++) {
			Coordinates coords = new Coordinates(volumeGaussiansRefCoords[i*3],volumeGaussiansRefCoords[i*3+1],volumeGaussiansRefCoords[i*3+2]);
			VolumeGaussian vg = new VolumeGaussian(referenceVolGaussians.get(i));
			vg.setReferenceVector(new Coordinates(coords.x, coords.y, coords.z));
			Coordinates shift = new Coordinates(volumeGaussiansShiftCoords[i*3],volumeGaussiansShiftCoords[i*3+1],volumeGaussiansShiftCoords[i*3+2]);
			vg.setShiftVector(shift);
			volumeGaussians.add(vg);
		}
		


		for(int i=0;i<nrOfHydrogens;i++) {
			hydrogens.add(new Coordinates(hydrogensCoords[i*3],hydrogensCoords[i*3+1],hydrogensCoords[i*3+2]));
		}



		return new MolecularVolume(atomicGaussians,ppGaussians, volumeGaussians,hydrogens);
	}
	
	
	public static MolecularVolume decodeFull(String string, StereoMolecule refMol)  {
		String[] splitString = string.split("  ");
		int nrOfAtomicGaussians = Integer.decode(splitString[0].trim());
		int firstIndex = 1;
		int lastIndex = 1+nrOfAtomicGaussians;
		List<AtomicGaussian> atomicGaussians = new ArrayList<AtomicGaussian>();
		List<PPGaussian> ppGaussians = new ArrayList<PPGaussian>();
		List<VolumeGaussian> volumeGaussians = new ArrayList<VolumeGaussian>();
		List<Coordinates> hydrogens = new ArrayList<Coordinates>();
		
		for(int i=firstIndex;i<lastIndex;i++) {
			atomicGaussians.add(AtomicGaussian.fromString(splitString[i].trim()));
		}
		int nrOfPPGaussians = Integer.decode(splitString[lastIndex]);
		firstIndex = lastIndex+1;
		lastIndex = firstIndex + nrOfPPGaussians;
		for(int i=firstIndex;i<lastIndex;i++) {
			ppGaussians.add(PPGaussian.fromString(splitString[i],refMol));
		}
		
		int nrOfVolumeGaussians = Integer.decode(splitString[lastIndex]);
		firstIndex = lastIndex+1;
		lastIndex = firstIndex + nrOfVolumeGaussians;
		for(int i=firstIndex;i<lastIndex;i++) {
			volumeGaussians.add(VolumeGaussian.fromString(splitString[i]));
		}


		double[] coords = EncoderFloatingPointNumbers.decode(splitString[splitString.length-1]);
		int nrOfHydrogens = coords.length/3;
		for(int i=0;i<nrOfHydrogens;i++) {
			hydrogens.add(new Coordinates(coords[i*3],coords[i*3+1],coords[i*3+2]));
			
		}
		return new MolecularVolume(atomicGaussians,ppGaussians,volumeGaussians,hydrogens);
	}
}

