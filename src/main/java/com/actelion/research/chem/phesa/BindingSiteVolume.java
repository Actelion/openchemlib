package com.actelion.research.chem.phesa;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.phesa.pharmacophore.SimplePPGaussian;

public class BindingSiteVolume {
	
	private List<SimplePPGaussian> ppGaussians;
	private List<AtomicGaussian> atomicGaussians;
	
	public BindingSiteVolume() {
		ppGaussians = new ArrayList<>();
		atomicGaussians = new ArrayList<>();
		
	}
	
	public void addPharmacophorePoint(SimplePPGaussian ppGaussian) {
		ppGaussians.add(ppGaussian);
	}
	
	public void addAtomVolume(AtomicGaussian atomGaussian) {
		atomicGaussians.add(atomGaussian);
	}
	
	public String encode() {
		StringBuilder sb = new StringBuilder();
		sb.append(Integer.toString(atomicGaussians.size()));
		sb.append("  ");
		atomicGaussians.forEach(e -> {
			sb.append(e.encode());
			sb.append("  ");
		});
		sb.append(Integer.toString(ppGaussians.size()));
		sb.append("  ");
		ppGaussians.forEach(e -> {
			sb.append(e.encode());
			sb.append("  ");
		});
		return sb.toString();
	}
	
	public static BindingSiteVolume decode(String s) {
		String[] splitString = s.split("  ");
		int nrOfAtomicGaussians = Integer.decode(splitString[0].trim());
		int firstIndex = 1;
		int lastIndex = 1+nrOfAtomicGaussians;
		List<AtomicGaussian> atomicGaussians = new ArrayList<AtomicGaussian>();
		List<SimplePPGaussian> ppGaussians = new ArrayList<SimplePPGaussian>();
		
		for(int i=firstIndex;i<lastIndex;i++) {
			atomicGaussians.add(AtomicGaussian.fromString(splitString[i].trim()));
		}
		int nrOfPPGaussians = Integer.decode(splitString[lastIndex]);
		firstIndex = lastIndex+1;
		lastIndex = firstIndex + nrOfPPGaussians;
		for(int i=firstIndex;i<lastIndex;i++) {
			ppGaussians.add(new SimplePPGaussian(splitString[i]));
		}
		BindingSiteVolume receptorVol = new BindingSiteVolume();
		atomicGaussians.forEach(e -> receptorVol.addAtomVolume(e));
		ppGaussians.forEach(e -> receptorVol.addPharmacophorePoint(e));
		
		return receptorVol;
		
	}

	public List<SimplePPGaussian> getPPGaussians() {
		return ppGaussians;
	}

	public void setPPGaussians(List<SimplePPGaussian> ppGaussians) {
		this.ppGaussians = ppGaussians;
	}

	public List<AtomicGaussian> getAtomicGaussians() {
		return atomicGaussians;
	}

	public void setAtomicGaussians(List<AtomicGaussian> atomicGaussians) {
		this.atomicGaussians = atomicGaussians;
	}
	
	
	

}
