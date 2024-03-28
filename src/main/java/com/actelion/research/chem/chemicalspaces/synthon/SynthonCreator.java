package com.actelion.research.chem.chemicalspaces.synthon;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.chem.io.RXNFileCreator;
import com.actelion.research.chem.reaction.Reaction;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.util.*;

public class SynthonCreator {

	/**
	 * only works for reactions with a single product
	 * @param rxn
	 * @return
	 * @throws Exception
	 */
	public static Reaction[] create(Reaction rxn) throws Exception {
		if(rxn.getProducts()>1)
			throw new Exception("only reactions with one product are supported");
		Reaction[] synthonTransformations = new Reaction[rxn.getReactants()];
		Map<Integer,Integer> mappedAtomToReactant = new HashMap<>(); //stores the information on which mapped atom in the product ist contributed by which reactant
		for(int r=0;r<rxn.getReactants();r++) {
			if(rxn.getReactants()==0)
				throw new Exception("cannot create Synthons for reactions with one reactant");
			StereoMolecule reactant = rxn.getReactant(r);
			for(int a=0;a<reactant.getAtoms();a++) {
				int reactantMapNo = reactant.getAtomMapNo(a);
				if(reactantMapNo!=0)
					mappedAtomToReactant.put(reactantMapNo, r);
			}
		}
		//find bonds in products that are formed between a) atoms from different reactants or b) between an unmapped and a mapped atom
		//in a properly mapped reaction, unmapped atoms stem from a reagent or solvent
		StereoMolecule product = rxn.getProduct(0);
		//detect homogeneous rings: homogeneous rings are rings that either consist of only one reactant, or 
		//one reactant and unmapped (reagent) atoms. Bonds in homogeneous rings should not be cut!
		RingCollection rc = product.getRingSet();
		List<Integer> homogeneousRingBonds = new ArrayList<>();
		for(int ringNo=0;ringNo<rc.getSize();ringNo++) {
			Set<Integer> ringReactant = new HashSet<Integer>();
			int[] atoms = rc.getRingAtoms(ringNo);
			for(int a : atoms) {
				int mappedA = product.getAtomMapNo(a);
				if(mappedA!=0)
					ringReactant.add(mappedAtomToReactant.get(mappedA));
			}
			if(ringReactant.size()==1) {//homogeneous ring!
				int[] ringBonds = rc.getRingBonds(ringNo);
				Arrays.stream(ringBonds).forEach(e -> homogeneousRingBonds.add(e));
			}
			
		}
		List<Integer> bondsToCut = new ArrayList<>();
		for(int b=0;b<product.getBonds();b++) {
			if(homogeneousRingBonds.contains(b))
				continue;
			int bondAtom1 = product.getBondAtom(0, b);
			int bondAtom2 = product.getBondAtom(1, b);
			int mappedBondAtom1 = product.getAtomMapNo(bondAtom1);
			int mappedBondAtom2 = product.getAtomMapNo(bondAtom2);
			if(mappedBondAtom1==0) {
				if(mappedBondAtom2!=0) //bond between a mapped and an unmapped atom
					bondsToCut.add(b);
			}
			else if(mappedBondAtom2==0) { //bond between a mapped and an umapped atom
				bondsToCut.add(b);
			}
			else if(!mappedAtomToReactant.get(mappedBondAtom1).equals(mappedAtomToReactant.get(mappedBondAtom2))) { //atoms from two different reactants
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
			Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(rxn.getName() + "_" + reactantID + ".rxn"), StandardCharsets.UTF_8));
			new RXNFileCreator(reaction).writeRXNfile(writer);
			writer.close();

		
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
