package com.actelion.research.chem.chemicalspaces.ptree;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.SmilesCreator;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.chemicalspaces.synthon.SynthonReactor;

/**
 * reactor class that can handle reactions that are defined by synthons
 * @author joel
 *
 */


public class PharmTreeSynthonReactionHelper  {
	
	private List<StereoMolecule> synthons;
	private Map<Integer,Set<Integer>> ringWithLinks;
	private Map<Integer,Set<Integer>> ringReactantIndeces;
	private Map<Integer,List<Integer>> reactantsWithLinkers;
	private StereoMolecule genericProduct;
	private List<int[]> synthonsToGenericProductMap;
	private List<int[]> deletionMap;
	private List<Set<Integer>> productRings;
	private int[] bonds; // bond orders formed by linkers
	
	/**
	 * a reaction can be fully characterized by a list of synthons defining the bonds that are formed during the reaction
	 * special attention is given to ring systems that are formed ruing the reactions
	 * @param synthons
	 */
	public PharmTreeSynthonReactionHelper(List<StereoMolecule> synthons) {
		this.synthons = synthons;
		ringWithLinks = new HashMap<Integer,Set<Integer>>();
		ringReactantIndeces = new HashMap<Integer,Set<Integer>>();
		bonds = new int[4];
		classifyReaction();
		
	}
	
	
	private void classifyReaction() {
		reactantsWithLinkers = new HashMap<Integer,List<Integer>>();
		for(StereoMolecule reactant : synthons) {

			for(int a=0;a<reactant.getAtoms();a++) {
				if(reactant.getAtomLabel(a).equals("U")) {
					int bondOrder = reactant.getConnBondOrder(a, 0);
					bonds[0] = bondOrder;
				}
				else if(reactant.getAtomLabel(a).equals("Np")) {
					int bondOrder = reactant.getConnBondOrder(a, 0);
					bonds[1] = bondOrder;
				}
				else if(reactant.getAtomLabel(a).equals("Pu")) {
					int bondOrder = reactant.getConnBondOrder(a, 0);
					bonds[2] = bondOrder;
				}
				else if(reactant.getAtomLabel(a).equals("Am")) {
					int bondOrder = reactant.getConnBondOrder(a, 0);
					bonds[3] = bondOrder;
				}
			}
		}
		for(int r=0;r<synthons.size();r++) {
			reactantsWithLinkers.putIfAbsent(r, new ArrayList<Integer>());
			StereoMolecule reactant = synthons.get(r);
			for(int a=0;a<reactant.getAtoms();a++) {
				if(reactant.getAtomLabel(a).equals("U")) {
					reactantsWithLinkers.get(r).add(1);
				}
				else if(reactant.getAtomLabel(a).equals("Np")) {
					reactantsWithLinkers.get(r).add(2);
				}
				else if(reactant.getAtomLabel(a).equals("Pu")) {
					reactantsWithLinkers.get(r).add(3);
				}
				else if(reactant.getAtomLabel(a).equals("Am")) {
					reactantsWithLinkers.get(r).add(4);
				}
			}
		}
		Map<Integer,Integer> bondLabels = new HashMap<Integer,Integer>();
		synthonsToGenericProductMap = new ArrayList<int[]>();
		deletionMap = new ArrayList<int[]>();
		genericProduct = synthonReaction(synthons, synthonsToGenericProductMap, deletionMap,bondLabels);

		/*
		 * to do: get formed rings as StereoMolecules, so that they later can be represented as separate building block
		 * or try first: merge all functionalities of the rings into one node attached to one BB 
		 */
				
		
		RingCollection rc = genericProduct.getRingSet();
		productRings = new ArrayList<Set<Integer>>();
		for(int r=0;r<rc.getSize();r++) {
			productRings.add(Arrays.stream(rc.getRingAtoms(r)).boxed().collect(Collectors.toSet()));
		}
		int eductRings = 0;
		for(int r=0;r<synthons.size();r++) 
			eductRings+=synthons.get(r).getRingSet().getSize();
		if(eductRings==productRings.size())
			productRings = new ArrayList<Set<Integer>>();
		ringReactantIndeces = new HashMap<Integer,Set<Integer>>(); //information on which formed ring contains which reactants  
		for(int r=0;r<synthons.size();r++) {
			StereoMolecule reactant = synthons.get(r);
			for(int a=0;a<reactant.getAtoms();a++) {
				if(reactant.getAtomicNo(a)>=SynthonReactor.CONNECTOR_OFFSET)
					continue;
				for(int fr=0;fr<productRings.size();fr++) {
					if(productRings.get(fr).contains(mapSynthonToGenericProductIndex(r,a))) {
						ringReactantIndeces.putIfAbsent(fr, new HashSet<Integer>());
						ringReactantIndeces.get(fr).add(r);
					}
				}

					
			}
		}

		for(int r=0;r<rc.getSize();r++) {
			int[] ringBonds = rc.getRingBonds(r);
			for(int b : ringBonds) {
				if(bondLabels.containsKey(b)) {
					ringWithLinks.putIfAbsent(r, new HashSet<Integer>());
					ringWithLinks.get(r).add(bondLabels.get(b));
				}
			}
		}
		
	}
	
	/**
	 * forms the product from the given building blocks and the generic reactants that fully
	 * define the reaction, also calculates that mapping from the reactant indeces to the product indeces 
	 * the bond that is formed by connecting the two atoms adjacent to R1 is labeled 1, 2 for R2, etc.
	 * @param buildingBlocks
	 * @param genericReactants
	 * @param atomMap
	 * @return
	 */
	private  StereoMolecule synthonReaction(List<StereoMolecule> bbs,  List<int[]> atomMap, List<int[]> deletionMap,
			Map<Integer,Integer> bondLabels) {
		List<StereoMolecule> buildingBlocks = new ArrayList<StereoMolecule>();
		bbs.stream().forEach(e -> buildingBlocks.add(new StereoMolecule(e)));
		buildingBlocks.forEach(e -> e.ensureHelperArrays(Molecule.cHelperCIP));
		HashMap<Integer,List<int[]>> rgrps = new HashMap<Integer,List<int[]>>();
		for(int m=0;m<buildingBlocks.size();m++) {
			StereoMolecule bb = buildingBlocks.get(m);
			for(int a=0;a<bb.getAtoms();a++) {
				if(bb.getAtomLabel(a).equals("U")) {
					rgrps.putIfAbsent(1, new ArrayList<int[]>());
					rgrps.get(1).add(new int[]{m,bb.getConnAtom(a, 0)});
					bb.markAtomForDeletion(a);
				}
				else if(bb.getAtomLabel(a).equals("Np")) {
					rgrps.putIfAbsent(2, new ArrayList<int[]>());
					rgrps.get(2).add(new int[]{m,bb.getConnAtom(a, 0)});
					bb.markAtomForDeletion(a);
				}
				else if(bb.getAtomLabel(a).equals("Pu")) {
					rgrps.putIfAbsent(3, new ArrayList<int[]>());
					rgrps.get(3).add(new int[]{m,bb.getConnAtom(a, 0)});
					bb.markAtomForDeletion(a);
				}
				else if(bb.getAtomLabel(a).equals("Am")) {
					rgrps.putIfAbsent(4, new ArrayList<int[]>());
					rgrps.get(4).add(new int[]{m,bb.getConnAtom(a, 0)});
					bb.markAtomForDeletion(a);
				}
				
			}
		}
		// remove R groups
		for(int m=0;m<buildingBlocks.size();m++) {
			StereoMolecule bb = buildingBlocks.get(m);
			int[] map = bb.deleteMarkedAtomsAndBonds();
			deletionMap.add(map);
			final int m_ = m;
			rgrps.forEach((k,v) -> {
				v.stream().forEach(e -> {
					if(e[0]==m_)
						e[1] = map[e[1]];
				});
				
			});
			
		}

		StereoMolecule product = new StereoMolecule();
		atomMap.add(product.addMolecule(buildingBlocks.get(0)));
		for(int m=1;m<buildingBlocks.size();m++) {
			int[] map = product.addMolecule(buildingBlocks.get(m));
			atomMap.add(map);
			final int m_ = m;
			rgrps.forEach((k,v) -> {
				v.stream().forEach(e -> {
					if(e[0]==m_)
						e[1] = map[e[1]];
				});
				
			});
		
		}
		rgrps.forEach((k,v) -> {
			int atom1 = v.get(0)[1];
			int atom2 = v.get(1)[1];
			int bond = product.addBond(atom1, atom2);
			product.setBondOrder(bond, bonds[k-1]);
			bondLabels.put(bond, k);
		});
		int[] map = product.getHandleHydrogenMap();

		atomMap.stream().forEach(e -> {
			Arrays.stream(e).map(i -> map[i]);});

		product.ensureHelperArrays(Molecule.cHelperCIP);
		return product;
		
	}
	



	public List<StereoMolecule> getSynthons() {
		return synthons;
	}


	public Map<Integer, Set<Integer>> getRingWithLinks() {
		return ringWithLinks;
	}


	public Map<Integer, Set<Integer>> getRingReactantIndeces() {
		return ringReactantIndeces;
	}


	public StereoMolecule getGenericProduct() {
		return genericProduct;
	}


	public List<int[]> getSynthonsToGenericProductMap() {
		return synthonsToGenericProductMap;
	}
	
	public int mapSynthonToGenericProductIndex(int synthonID, int index) {
		int newIndex = deletionMap.get(synthonID)[index];
		return synthonsToGenericProductMap.get(synthonID)[newIndex];
	}


	public List<Set<Integer>> getFormedRings() {
		return productRings;
	}


	public Map<Integer, List<Integer>> getReactantsWithLinkers() {
		return reactantsWithLinkers;
	}
	
	public int[] getBondOrders() {
		return bonds;
	}
	

}
