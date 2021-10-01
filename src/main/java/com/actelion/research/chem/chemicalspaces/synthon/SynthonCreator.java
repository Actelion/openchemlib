package com.actelion.research.chem.chemicalspaces.synthon;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.chem.io.RXNFileCreator;
import com.actelion.research.chem.reaction.Reaction;

public class SynthonCreator {
	
	
	public static Reaction[] create(Reaction rxn) throws Exception {
		if(rxn.getProducts()>1)
			throw new Exception("only reactions with one product are supported");
		Reaction[] synthonTransformations = new Reaction[rxn.getReactants()];
		Map<Integer,Integer> mappedAtomToReactant = new HashMap<Integer,Integer>(); //stores the information on which mapped atom in the product ist contributed by which reactant
		for(int r=0;r<rxn.getReactants();r++) {
			if(rxn.getReactants()==0)
				throw new Exception("cannot create Synthons for reactions with one reactant");
			StereoMolecule reactant = rxn.getReactant(r);
			for(int a=0;a<reactant.getAtoms();a++) {
				int reactantMapNo = reactant.getAtomMapNo(a);
				if(reactantMapNo==0)
					continue;
				else {
					mappedAtomToReactant.put(reactantMapNo, r);
				}
			}
		}
		//find bonds in products that are formed between a) atoms from different reactants or b) between an unmapped and a mapped atom
		//in a properly mapped reaction, unmapped atoms stem from a reagent or solvent
		for(int p=0;p<rxn.getProducts();p++) {
			StereoMolecule product = rxn.getProduct(p);
			List<Integer> bondsToCut = new ArrayList<>();
			for(int b=0;b<product.getBonds();b++) {
				int bondAtom1 = product.getBondAtom(0, b);
				int bondAtom2 = product.getBondAtom(1, b);
				int mappedBondAtom1 = product.getAtomMapNo(bondAtom1);
				int mappedBondAtom2 = product.getAtomMapNo(bondAtom2);
				if(mappedBondAtom1==0) {
					if(mappedBondAtom2==0)
						continue;
					else 
						bondsToCut.add(b);
						
				}
				else if(mappedBondAtom2==0) {
					bondsToCut.add(b);
					continue;
				}
				else if(mappedAtomToReactant.get(mappedBondAtom1)!=mappedAtomToReactant.get(mappedBondAtom2)) {
					bondsToCut.add(b);
				}
			}

			StereoMolecule decomposedProduct = new StereoMolecule(product); 
			decomposedProduct.ensureHelperArrays(Molecule.cHelperCIP);
			int counter = 0;
			for(int cutBond : bondsToCut) {
				//first check if cutting this bond would result in isolating a fragment that is made of nonmapped atoms --> stem from a reagent
				//we want to assign these atoms to one of the synthons, so such a cut is forbidden
				StereoMolecule tmpMol = new StereoMolecule(decomposedProduct);
				tmpMol.ensureHelperArrays(Molecule.cHelperNeighbours);
				tmpMol.setBondType(cutBond, Molecule.cBondTypeDeleted);
				tmpMol.deleteMarkedAtomsAndBonds();
				StereoMolecule[] tmpFrags = tmpMol.getFragments();
				boolean acceptCut = true;
				for(StereoMolecule tmpFrag : tmpFrags) {
					tmpFrag.ensureHelperArrays(Molecule.cHelperNeighbours);
					int mapSum = 0;
					for(int a=0;a<tmpFrag.getAtoms();a++) {
						mapSum+=tmpFrag.getAtomMapNo(a);
					}
					if(mapSum==0) {
						acceptCut=false;
						break;
					}
				}
				if(acceptCut) {
					int bondType = decomposedProduct.getBondType(cutBond);
					decomposedProduct.setBondType(cutBond, Molecule.cBondTypeDeleted);
					int connector1 = decomposedProduct.addAtom(92+counter);
					int connector2 = decomposedProduct.addAtom(92+counter);
					int newBond1 = decomposedProduct.addBond(decomposedProduct.getBondAtom(0, cutBond), connector1);
					decomposedProduct.setBondType(newBond1, bondType);
					int newBond2 = decomposedProduct.addBond(connector2,decomposedProduct.getBondAtom(1, cutBond));
					decomposedProduct.setBondType(newBond2, bondType);
					counter++;
				}
				
			}
			decomposedProduct.deleteMarkedAtomsAndBonds();
			decomposedProduct.ensureHelperArrays(Molecule.cHelperCIP);
			StereoMolecule[] frags = decomposedProduct.getFragments();
			for(StereoMolecule frag : frags) {
				CoordinateInventor inventor = new CoordinateInventor();
				inventor.setRandomSeed(0x1234567890L);  // create reproducible coordinates
				inventor.invent(frag);
				Reaction reaction = new Reaction();
				int reactantID = -1;
				for(int a=0;a<frag.getAtoms();a++) {
					if(frag.getAtomMapNo(a)!=0) {
						reactantID = mappedAtomToReactant.get(frag.getAtomMapNo(a));
					}
						
				}
				if(reactantID==-1)
					throw new Exception("could not process reaction");
				frag.ensureHelperArrays(Molecule.cHelperCIP);
				reaction.addReactant(rxn.getReactant(reactantID));
				reaction.addProduct(frag);
				synthonTransformations[reactantID] = reaction;
				//Writer writer = new BufferedWriter(new FileWriter(rxn.getName() + "_" + i + ".rxn"));
				//new RXNFileCreator(reaction).writeRXNfile(writer);
				//writer.close();

				

				
			}
		
		}
		return synthonTransformations;
		
		
		
		
	
	}
	/**
	 * does not allow cuts through ring bonds
	 * @param rxn
	 * @return
	 * @throws Exception
	 */
	public static Reaction[] createStrict(Reaction rxn) throws Exception {
		if(rxn.getProducts()>1)
			throw new Exception("only reactions with one product are supported");
		Reaction[] synthonTransformations = new Reaction[rxn.getReactants()];
		Map<Integer,Integer> mappedAtomToReactant = new HashMap<Integer,Integer>(); //stores the information on which mapped atom in the product is contributed by which reactant
		for(int r=0;r<rxn.getReactants();r++) {
			if(rxn.getReactants()==0)
				throw new Exception("cannot create Synthons for reactions with one reactant");
			StereoMolecule reactant = rxn.getReactant(r);
			for(int a=0;a<reactant.getAtoms();a++) {
				int reactantMapNo = reactant.getAtomMapNo(a);
				if(reactantMapNo==0)
					continue;
				else {
					mappedAtomToReactant.put(reactantMapNo, r);
				}
			}
		}
		//find bonds in products that are formed between a) atoms from different reactants or b) between an unmapped and a mapped atom
		//in a properly mapped reaction, unmapped atoms stem from a reagent or solvent
		for(int p=0;p<rxn.getProducts();p++) {
			StereoMolecule product = rxn.getProduct(p);
			List<Integer> bondsToCut = new ArrayList<>();
			for(int b=0;b<product.getBonds();b++) {
				int bondAtom1 = product.getBondAtom(0, b);
				int bondAtom2 = product.getBondAtom(1, b);
				int mappedBondAtom1 = product.getAtomMapNo(bondAtom1);
				int mappedBondAtom2 = product.getAtomMapNo(bondAtom2);
				if(mappedBondAtom1==0) {
					if(mappedBondAtom2==0)
						continue;
					else 
						bondsToCut.add(b);
						
				}
				else if(mappedBondAtom2==0) {
					bondsToCut.add(b);
					continue;
				}
				else if(mappedAtomToReactant.get(mappedBondAtom1)!=mappedAtomToReactant.get(mappedBondAtom2)) {
					bondsToCut.add(b);
				}
			}

			StereoMolecule decomposedProduct = new StereoMolecule(product); 
			decomposedProduct.ensureHelperArrays(Molecule.cHelperCIP);
			int counter = 0;
			for(int cutBond : bondsToCut) {
				//first check if cutting this bond would result in isolating a fragment that is made of nonmapped atoms --> stem from a reagent
				//we want to assign these atoms to one of the synthons, so such a cut is forbidden
				StereoMolecule tmpMol = new StereoMolecule(decomposedProduct);
				tmpMol.ensureHelperArrays(Molecule.cHelperNeighbours);
				tmpMol.setBondType(cutBond, Molecule.cBondTypeDeleted);
				tmpMol.deleteMarkedAtomsAndBonds();
				StereoMolecule[] tmpFrags = tmpMol.getFragments();
				boolean acceptCut = true;
				for(StereoMolecule tmpFrag : tmpFrags) {
					tmpFrag.ensureHelperArrays(Molecule.cHelperNeighbours);
					int mapSum = 0;
					for(int a=0;a<tmpFrag.getAtoms();a++) {
						mapSum+=tmpFrag.getAtomMapNo(a);
					}
					if(mapSum==0) {
						acceptCut=false;
						break;
					}
				}
				if(acceptCut) {
					int bondType = decomposedProduct.getBondType(cutBond);
					decomposedProduct.setBondType(cutBond, Molecule.cBondTypeDeleted);
					int connector1 = decomposedProduct.addAtom(92+counter);
					int connector2 = decomposedProduct.addAtom(92+counter);
					int newBond1 = decomposedProduct.addBond(decomposedProduct.getBondAtom(0, cutBond), connector1);
					decomposedProduct.setBondType(newBond1, bondType);
					int newBond2 = decomposedProduct.addBond(connector2,decomposedProduct.getBondAtom(1, cutBond));
					decomposedProduct.setBondType(newBond2, bondType);
					counter++;
				}
				
			}
			decomposedProduct.deleteMarkedAtomsAndBonds();
			decomposedProduct.ensureHelperArrays(Molecule.cHelperCIP);
			StereoMolecule[] frags = decomposedProduct.getFragments();
			for(StereoMolecule frag : frags) {
				CoordinateInventor inventor = new CoordinateInventor();
				inventor.setRandomSeed(0x1234567890L);  // create reproducible coordinates
				inventor.invent(frag);
				Reaction reaction = new Reaction();
				int reactantID = -1;
				for(int a=0;a<frag.getAtoms();a++) {
					if(frag.getAtomMapNo(a)!=0) {
						reactantID = mappedAtomToReactant.get(frag.getAtomMapNo(a));
					}
						
				}
				if(reactantID==-1)
					throw new Exception("could not process reaction");
				frag.ensureHelperArrays(Molecule.cHelperCIP);
				reaction.addReactant(rxn.getReactant(reactantID));
				reaction.addProduct(frag);
				synthonTransformations[reactantID] = reaction;
				//Writer writer = new BufferedWriter(new FileWriter(rxn.getName() + "_" + i + ".rxn"));
				//new RXNFileCreator(reaction).writeRXNfile(writer);
				//writer.close();

				

				
			}
		
		}
		return synthonTransformations;
		
		
		
		
	
	}
	
	
	/*
	public static void main(String[] args) {
		List<Reaction> htmcReactions = ReactionList.INSTANCE.getReactionsHTMC();
		for(Reaction reaction : htmcReactions) {
			System.out.println(reaction.getName());
			try {
				create(reaction, false);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	*/
	
	


}
