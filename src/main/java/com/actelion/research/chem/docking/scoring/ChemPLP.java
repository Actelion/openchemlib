package com.actelion.research.chem.docking.scoring;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.actelion.research.calc.combinatorics.CombinationGenerator;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.alignment3d.KabschAlignment;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.docking.LigandPose;
import com.actelion.research.chem.docking.scoring.chemscore.HBTerm;
import com.actelion.research.chem.docking.scoring.chemscore.MetalTerm;
import com.actelion.research.chem.docking.scoring.chemscore.SimpleMetalTerm;
import com.actelion.research.chem.docking.scoring.plp.PLPTerm;
import com.actelion.research.chem.docking.scoring.plp.REPTerm;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import com.actelion.research.chem.io.pdb.converter.MoleculeGrid;
import com.actelion.research.chem.phesa.pharmacophore.ChargedGroupDetector;
import com.actelion.research.chem.phesa.pharmacophore.IonizableGroupDetector;
import com.actelion.research.chem.phesa.pharmacophore.PharmacophoreCalculator;
import com.actelion.research.chem.phesa.pharmacophore.pp.ChargePoint;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;

/**
 * Implementation of ChemPLP scoring function as described in: doi: 10.1021/ci800298z
 * THIS SCORING FUNCTION REQUIRES EXPLICIT HYDROGENS TO BE PRESENT!
 * @author joel
 *
 */

public class ChemPLP extends AbstractScoringEngine {
	
	private static final double METAL_INTERACTION_CUTOFF = 2.6;
	private static final double STRAIN_CUTOFF = 10;
	private static final double METAL_OPTIMAL_DIST = 2.2;
			
	
	private Set<Integer> receptorAcceptors;
	private Set<Integer> receptorDonorHs;
	private Set<Integer> receptorDonors;
	private Map<Integer, Double>receptorDonorHPos;
	private Map<Integer, Double> receptorAcceptorNeg;
	private Set<Integer> receptorMetals;
	
	private Set<Integer> ligandAcceptors;
	private Set<Integer> ligandDonorHs;
	private Set<Integer> ligandDonors;
	private Map<Integer,Double> ligandDonorHPos;
	private Map<Integer, Double> ligandAcceptorNeg;
	
	private List<PotentialEnergyTerm> plp;
	private List<PotentialEnergyTerm> chemscoreHbond;
	private List<PotentialEnergyTerm> chemscoreMetal;
	private ForceFieldMMFF94 ff;
	private double e0;
	
	private Map<Integer,List<Coordinates>> metalInteractionSites;
	
	private static final Set<Integer> SIMPLE_METAL_ATOMS = new HashSet(Arrays.asList(12,20)); //Mg and Ca
	

	public ChemPLP(Molecule3D receptor, Set<Integer> bindingSiteAtoms, MoleculeGrid grid) {
		super(receptor, bindingSiteAtoms, grid);
		receptorAcceptors = new HashSet<>();
		
		receptorDonorHs = new HashSet<>();
		receptorDonorHPos = new HashMap<Integer, Double>();
		receptorAcceptorNeg = new HashMap<Integer, Double>();
		receptorMetals = new HashSet<>();
		receptorDonors = new HashSet<>();
		identifyHBondFunctionality(receptor,receptorAcceptors,receptorDonorHs, receptorDonors, receptorMetals,receptorAcceptorNeg,
				receptorDonorHPos);
		
		metalInteractionSites = new HashMap<Integer,List<Coordinates>>();
		
		for(int met : receptorMetals) 
			metalInteractionSites.put(met, processMetalCoordination(receptorConf,met,receptorAcceptors));
			
		
		
	}
	

	@Override
	public double getFGValue(double[] grad) {
		double energy = getBumpTerm();
		for(PotentialEnergyTerm term : chemscoreHbond)
			energy+=term.getFGValue(grad);
		for(PotentialEnergyTerm term : getMetalTerm(chemscoreMetal)) 	
			energy+=term.getFGValue(grad);
		for(PotentialEnergyTerm term : plp) 
			energy+=term.getFGValue(grad);
		ff.setState(candidatePose.getCartState());
		double ffEnergy = ff.getTotalEnergy();
		if((ffEnergy-e0)>STRAIN_CUTOFF) {
			energy+=ffEnergy-e0;
			ff.addGradient(grad);
		}
		for(PotentialEnergyTerm term : constraints)
			energy+=term.getFGValue(grad);

		return energy;
	}
	
	@Override 
	public void updateState() {
		ff.setState(candidatePose.getCartState());
	}
	
	@Override
	public double getScore() {
		double[] grad = new double[3*candidatePose.getLigConf().getMolecule().getAllAtoms()];
		double energy = getBumpTerm();
		for(PotentialEnergyTerm term : chemscoreHbond)
			energy+=term.getFGValue(grad);
		for(PotentialEnergyTerm term : getMetalTerm(chemscoreMetal)) 
			energy+=term.getFGValue(grad);
		for(PotentialEnergyTerm term : plp) 
			energy+=term.getFGValue(grad);
		Map<String, Object> ffOptions = new HashMap<String, Object>();
		ffOptions.put("dielectric constant", 80.0);

		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		StereoMolecule toOptimize = new StereoMolecule(candidatePose.getLigConf().getMolecule());
		toOptimize.ensureHelperArrays(Molecule.cHelperCIP);
		ForceFieldMMFF94 ff = new ForceFieldMMFF94(toOptimize, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
		double e = ff.getTotalEnergy();
		ff.minimise();
		double e0 = ff.getTotalEnergy();
		double deltaE = e-e0;
		if(deltaE>STRAIN_CUTOFF) {
			energy += deltaE-STRAIN_CUTOFF;
		}
	

		return energy;
	}
	

	@Override
	public void init(LigandPose candidatePose, double e0) {
		this.e0 = e0;
		this.candidatePose = candidatePose;
		
		plp = new ArrayList<>();
		chemscoreHbond = new ArrayList<>();
		chemscoreMetal = new ArrayList<>();
		
		ligandAcceptors = new HashSet<>();
		ligandDonorHs = new HashSet<>();
		ligandDonors = new HashSet<>();
		ligandDonorHPos = new HashMap<Integer, Double>();
		ligandAcceptorNeg = new HashMap<Integer, Double>();
		constraints = new ArrayList<>();
		
		
		StereoMolecule ligand = candidatePose.getLigConf().getMolecule();
		
		
		Map<String, Object> ffOptions = new HashMap<String, Object>();
		ffOptions.put("dielectric constant", 80.0);

		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		ff = new ForceFieldMMFF94(ligand, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
		StereoMolecule receptor = receptorConf.getMolecule();
		identifyHBondFunctionality(ligand,ligandAcceptors,ligandDonorHs, ligandDonors, new HashSet<Integer>(),ligandAcceptorNeg,
				ligandDonorHPos);

		for(int p : bindingSiteAtoms) {
			if(receptor.getAtomicNo(p)==1) { // receptor hydrogen atom
				if(receptorDonorHs.contains(p)) {
					// receptor atom is donor hydrogen
					int d = receptor.getConnAtom(p, 0);
					boolean chargedP = receptorDonorHPos.keySet().contains(p);
					for(int l=0;l<ligand.getAtoms();l++) {
						if(ligandAcceptors.contains(l)) {
							final int li = l;
							int[] acceptorNeighbours = IntStream.range(0, ligand.getConnAtoms(li)).map(i -> ligand.getConnAtom(li, i)).toArray();
							boolean chargedL = ligandAcceptorNeg.keySet().contains(l);
							double scale = 1.0;
							if(chargedP && chargedL) {
								scale += (ligandAcceptorNeg.get(l) * receptorDonorHPos.get(p));
							}
							HBTerm hbTerm = HBTerm.create(receptorConf, candidatePose.getLigConf(), l, d, p, true, false, acceptorNeighbours, scale);
							chemscoreHbond.add(hbTerm);
						}	
				}
			}
			}
			else { //receptor heavy atom
				if(receptorDonors.contains(p)) { // receptor donor heavy atom -> only plp terms 
					for(int l=0;l<ligand.getAtoms();l++) { //only consider ligand heavy atoms
						if(ligandAcceptors.contains(l)) {  //plp hbond donor-acceptor
							PLPTerm plpTerm = PLPTerm.create(receptorConf, candidatePose.getLigConf(), p, l, PLPTerm.HBOND_TERM);
							plp.add(plpTerm);
						}
						else if(ligandDonors.contains(l)) {  //repulsive donor-donor
							REPTerm repTerm = REPTerm.create(receptorConf, candidatePose.getLigConf(), p, l);
							plp.add(repTerm);
						}
						else { //buried
							PLPTerm plpTerm = PLPTerm.create(receptorConf, candidatePose.getLigConf(), p, l, PLPTerm.BURIED_TERM);
							plp.add(plpTerm);
						}					
					}
				}
				else if(receptorAcceptors.contains(p)) { // receptor acceptor heavy atom
					int[] acceptorNeighbours = IntStream.range(0, receptor.getConnAtoms(p)).map(i -> receptor.getConnAtom(p, i)).toArray();
					boolean chargedP = receptorAcceptorNeg.keySet().contains(p);
					for(int l=0;l<ligand.getAllAtoms();l++) { 
						if(ligand.getAtomicNo(l)==1) { //ligand hydrogen atom
							if(ligandDonorHs.contains(l)) {
								boolean chargedL = ligandDonorHPos.keySet().contains(l);
								int d = ligand.getConnAtom(l, 0);
								double scale = 1.0;
								if(chargedP && chargedL)
									scale+= (receptorAcceptorNeg.get(p) * ligandDonorHPos.get(l));
								HBTerm hbTerm = HBTerm.create(receptorConf, candidatePose.getLigConf(), p, d,l, false, true, acceptorNeighbours, scale);
								chemscoreHbond.add(hbTerm);
							}
						}
						else { //ligand heavy atom
							if(ligandDonors.contains(l)) {  //plp hbond donor-acceptor
								PLPTerm plpTerm = PLPTerm.create(receptorConf, candidatePose.getLigConf(), p, l, PLPTerm.HBOND_TERM);
								plp.add(plpTerm);

							}
							else if(ligandAcceptors.contains(l)) {  //repulsive donor-donor
								REPTerm repTerm = REPTerm.create(receptorConf, candidatePose.getLigConf(), p, l);
								plp.add(repTerm);
							}
							else { //buried
								PLPTerm plpTerm = PLPTerm.create(receptorConf, candidatePose.getLigConf(), p, l, PLPTerm.BURIED_TERM);
								plp.add(plpTerm);
							}			
							
						}
					}	
					
				}
				else if(receptorMetals.contains(p)) {
						for(int l=0;l<ligand.getAtoms();l++) { 
							if(ligandDonors.contains(l)) {  //met-donor -> repulsive
								REPTerm repTerm = REPTerm.create(receptorConf, candidatePose.getLigConf(), p, l);
								plp.add(repTerm);
							}
							else if(ligandAcceptors.contains(l)) { //attractive met-acc interaction
								PLPTerm plpTerm = PLPTerm.create(receptorConf, candidatePose.getLigConf(), p, l, PLPTerm.METAL_TERM);
								plp.add(plpTerm);
							}
							else { //buried met-nonp interaction
								PLPTerm plpTerm = PLPTerm.create(receptorConf, candidatePose.getLigConf(), p, l, PLPTerm.BURIED_TERM);
								plp.add(plpTerm);
							}	
						}
					
					if(SIMPLE_METAL_ATOMS.contains(receptor.getAtomicNo(p))) {
						for(int l : ligandAcceptors) {
							double scale = 1.0;
							if(ligandAcceptorNeg.keySet().contains(l))
								scale += ligandAcceptorNeg.get(l);
							int[] acceptorNeighbours = IntStream.range(0, ligand.getConnAtoms(l)).map(i -> ligand.getConnAtom(l, i)).toArray();
							SimpleMetalTerm metTerm = SimpleMetalTerm.create(receptorConf, candidatePose.getLigConf(), 
									l, p, acceptorNeighbours, scale);
							chemscoreMetal.add(metTerm);
						}
					}
	
					
					
					else { //standard metal term;
						List<Coordinates> interactionSites = metalInteractionSites.get(p);
						for(int l : ligandAcceptors) {
							double scale = 1.0;
							if(ligandAcceptorNeg.keySet().contains(l))
								scale += ligandAcceptorNeg.get(l);
							int[] acceptorNeighbours = IntStream.range(0, ligand.getConnAtoms(l)).map(i -> ligand.getConnAtom(l, i)).toArray();
							for(Coordinates site : interactionSites) {
								MetalTerm metTerm = MetalTerm.create(candidatePose.getLigConf(), l, receptorConf,p, acceptorNeighbours, site,scale);
								chemscoreMetal.add(metTerm);
							}
					}
					
					
					}
					

					
				}
				else { // non-polar heavy atom
					for(int l=0;l<ligand.getAtoms();l++) { 
						if(ligandDonors.contains(l)) {  //nonpolar-donor -> buried
							PLPTerm plpTerm = PLPTerm.create(receptorConf, candidatePose.getLigConf(), p, l, PLPTerm.BURIED_TERM);
							plp.add(plpTerm);
						}
						else if(ligandAcceptors.contains(l)) { //nonpolar-acceptor -> buried
							PLPTerm plpTerm = PLPTerm.create(receptorConf, candidatePose.getLigConf(), p, l, PLPTerm.BURIED_TERM);
							plp.add(plpTerm);
						}
						else { // nonpolar-nonpolar -> nonpolar
							PLPTerm plpTerm = PLPTerm.create(receptorConf, candidatePose.getLigConf(), p, l, PLPTerm.NONPOLAR_TERM);
							plp.add(plpTerm);
						}	
				}
				}
			}
			}
		
		
	}
	
	//tries tetrahedral or octahedral coordination at metal and choses the one that gives the better fit with the alignment
	private static List<Coordinates> processMetalCoordination(Conformer receptor, int metalAtom, Set<Integer> receptorAcceptors) {
		List<Coordinates> interactionPoints = new ArrayList<Coordinates>();
		Coordinates metalCoordinates = receptor.getCoordinates(metalAtom);
		List<Coordinates> interactionSites = new ArrayList<Coordinates>();
		interactionSites.add(new Coordinates(metalCoordinates));
		for(int acceptor : receptorAcceptors) {
			Coordinates acceptorCoords = receptor.getCoordinates(acceptor);
			if(acceptorCoords.distance(metalCoordinates)<METAL_INTERACTION_CUTOFF)
				interactionSites.add(acceptorCoords);
		}
		Coordinates[] occupiedSites = interactionSites.toArray(new Coordinates[interactionSites.size()]);
		//try tetrahedron 
		double overallMinimalRMSD = Double.MAX_VALUE;
		Coordinates[] fittedGeometry = null;;
		int[][] finalMapping = null;
		if(occupiedSites.length<5) {
			Coordinates[] idealGeom = getTetrahedron();
			int[][] bestMapping = new int[occupiedSites.length][2];
			overallMinimalRMSD = getBestMetalFit(idealGeom, occupiedSites,bestMapping,4);	
			fittedGeometry = idealGeom;
			finalMapping = bestMapping;
		}
		//try trigonal bipyramidal
		if(occupiedSites.length<6) {
			Coordinates[] idealGeom3 = getTrigonalBipyramidal();
			int[][] bestMapping3 = new int[occupiedSites.length][2];
			if(getBestMetalFit(idealGeom3, occupiedSites,bestMapping3,5)<overallMinimalRMSD) {
				fittedGeometry = idealGeom3;
				finalMapping = bestMapping3;
			}
		}
		//try octahedron 
		if(occupiedSites.length<7) {
			Coordinates[] idealGeom2 = getOctahedron();
			int[][] bestMapping2 = new int[occupiedSites.length][2];
			if(getBestMetalFit(idealGeom2, occupiedSites,bestMapping2,6)<overallMinimalRMSD) {
				fittedGeometry = idealGeom2;
				finalMapping = bestMapping2;
			}
		}
		// identify unoccpied sites
		if(fittedGeometry!=null) {

			KabschAlignment alignment = new KabschAlignment(occupiedSites,fittedGeometry,finalMapping);
			alignment.align();
			List<Integer> mappedPoints = Arrays.stream(finalMapping).map(e -> e[1]).collect(Collectors.toList());
			for(int i=0;i<fittedGeometry.length;i++) {
				if(!mappedPoints.contains(i))
					interactionPoints.add(fittedGeometry[i]);
			}
		}
		//place interaction points 2.2A from Metal
		for(Coordinates point : interactionPoints) {
			Coordinates v = point.subC(metalCoordinates);
			v.unit();
			v.scale(METAL_OPTIMAL_DIST);
			Coordinates newPoint = metalCoordinates.addC(v);
			point.x = newPoint.x;
			point.y = newPoint.y;
			point.z = newPoint.z;
		}
		
		return interactionPoints;
			
			
				
		}
	//for metal terms with coordination: only the most contributing term per interaction site is kept!
	private List<PotentialEnergyTerm> getMetalTerm(List<PotentialEnergyTerm> metalTerms) {
		double[] grad = new double[3*candidatePose.getLigConf().getMolecule().getAllAtoms()];
		Map<Coordinates,List<PotentialEnergyTerm>> interactionPt2Term = new HashMap<Coordinates,List<PotentialEnergyTerm>>();
		List<PotentialEnergyTerm> curatedTerms = new ArrayList<>();
		for(PotentialEnergyTerm term : metalTerms) {
			if(!(term instanceof MetalTerm)) {
				curatedTerms.add(term);
			}
			else {
				MetalTerm metalTerm = (MetalTerm) term;
				interactionPt2Term.putIfAbsent(metalTerm.getFitPoint(),new ArrayList<>());
				interactionPt2Term.get(metalTerm.getFitPoint()).add(metalTerm);
			}
		}
		for(Coordinates key : interactionPt2Term.keySet()) {
			PotentialEnergyTerm bestTerm = interactionPt2Term.get(key).stream().sorted((e1,e2) -> {
				return Double.compare(e1.getFGValue(grad), e2.getFGValue(grad));
			}).collect(Collectors.toList()).get(0);
			curatedTerms.add(bestTerm);
			
		}
		return curatedTerms;
		
	}


	
	private static double getBestMetalFit(Coordinates[] idealGeom, Coordinates[] occupiedSites,
			int[][] bestMapping, int coordinationNumber) {
		double overallMinimalRMSD = Double.MAX_VALUE;
		if((occupiedSites.length-1)<coordinationNumber) {
			List<int[]> assignments = enumerateMetalPosAssignments(coordinationNumber,occupiedSites.length-1);
			for(int[] ass : assignments) {
				int[][] mapping = new int[occupiedSites.length][2];
				int counter = 0;
				mapping[0][0] =0;
				mapping[0][1] =0;
				for(int i=1;i<occupiedSites.length;i++) {
					int[] m = new int[]{counter,ass[counter]};
					for(int j=0;j<m.length;j++) {
						m[j] = m[j] + 1;
					}
					mapping[i] = m;
					counter++;
				}
				Coordinates[] idealGeomCopy = Arrays.stream(idealGeom).map(e -> new Coordinates(e)).toArray(Coordinates[]::new);
				//shift mapping by one because metal is always on first position
				
				KabschAlignment alignment = new KabschAlignment(occupiedSites,idealGeomCopy,mapping);
				alignment.align();
				double rmsd = getRMSD(occupiedSites,idealGeomCopy,mapping);
				if(rmsd<overallMinimalRMSD) {
					overallMinimalRMSD = rmsd;
					for (int i = 0; i<mapping.length; i++){
					     for (int j = 0; j<mapping[i].length; j++){
					    	 bestMapping[i][j] = mapping[i][j];
					     }
					}
					
				}
			}
			}
		return overallMinimalRMSD;
			

	}
	
	
	private static double getRMSD(Coordinates[] c1, Coordinates[] c2, int[][] mapping) {
		double rmsd = 0.0;
		for(int[] m : mapping) {
			Coordinates coor1 = c1[m[0]];
			Coordinates coor2 = c2[m[1]];
			double dx = coor1.x-coor2.x;
			double dy = coor1.y-coor2.y;
			double dz = coor1.z-coor2.z;
			rmsd+=dx*dx + dy*dy + dz*dz;
		}
		rmsd/=mapping.length;
		return Math.sqrt(rmsd);
	}
	
	private static List<int[]> enumerateMetalPosAssignments(int coordinationNumber, int occupiedSites) {
		List<int[]> allAssignments = new ArrayList<int[]>();
		if(occupiedSites!=0) {
			if(occupiedSites==1) {
				int[] assignment = new int[] {0,0};
				allAssignments.add(assignment);
			}
			else {
				List<int[]> liComb = CombinationGenerator.getAllOutOf(coordinationNumber-1,occupiedSites-1);
				for(int[] r : liComb) {
					List<int[]> permutations = CombinationGenerator.getPermutations(r,r.length);
					for(int[] per : permutations) {
						int[] arr = new int[per.length+1];
						arr[0] = 0;
						IntStream.range(0, per.length).forEach(e -> {
						arr[e+1] = per[e]+1;});
						allAssignments.add(arr);
					}
				}
				}
		}
		return allAssignments;
	}
	
	
	public static void identifyHBondFunctionality(StereoMolecule mol, Set<Integer> acceptors, Set<Integer> donorHs, Set<Integer> donors,
			Set<Integer> metals, Map<Integer,Double> chargedAcceptors, Map<Integer,Double> chargedDonorHs) {

		for(int a=0;a<mol.getAllAtoms();a++) {
			if (mol.getAtomicNo(a)==7 || mol.getAtomicNo(a)==8) {
				if(PharmacophoreCalculator.isAcceptor(mol, a))
					acceptors.add(a);
				}
			else if(PharmacophoreCalculator.isDonorHydrogen(mol, a)) {
				donorHs.add(a);
				donors.add(mol.getConnAtom(a, 0));
			}
			else if (mol.isMetalAtom(a))
				metals.add(a);
		}
		ChargedGroupDetector detector = new ChargedGroupDetector(mol);
		detector.detect();
		List<ArrayList<Integer>> chargedGroups = detector.getChargedGroups();
		for(int a : acceptors) {
			for(List<Integer> chargedGroup : chargedGroups) {
				int invScale = chargedGroup.stream().mapToInt(e -> {
					if(mol.getAtomicNo(e)!=6 && mol.getAtomicNo(e)!=15 && mol.getAtomicNo(e)!=16)
						return 1;
					else 
						return 0;
				}).sum();
				double scale = 1.0/invScale;
				if(chargedGroup.contains(a)) {
					chargedAcceptors.put(a, scale);
				}
			}
		
		}
		for(int h : donorHs) {
			for(List<Integer> chargedGroup : chargedGroups) {
				int invScale = chargedGroup.stream().mapToInt(e -> {
					if(mol.getAtomicNo(e)!=6 && mol.getAtomicNo(e)!=15 && mol.getAtomicNo(e)!=16)
						return 1;
					else 
						return 0;
				}).sum();
				double scale = 1.0/invScale;;
				int d = mol.getConnAtom(h, 0);
				if(chargedGroup.contains(d)) {
					chargedDonorHs.put(h, scale);
				
				}
			}
		
		}
	}
	

	
	
	/*
	 * first coordinate always is the metal
	 */
	private static Coordinates[] getTetrahedron() {
		 return new Coordinates[] {
				    new Coordinates(0.0,0.0,0.0),
					new Coordinates(2.074,0,-0.733), new Coordinates(-1.037,1.796,-0.733),
					new Coordinates(-1.037,-1.796,-0.733), new Coordinates(0.0,0.0,2.2)
					};
	}
	
	private static Coordinates[] getOctahedron() {
		return new Coordinates[] {
				new Coordinates(0.0,0.0,0.0),
				new Coordinates(0,0,2.2), new Coordinates(0,0,-2.2), 
				new Coordinates(1.555,1.555,0.0), new Coordinates(1.555,-1.555,0.0),
				new Coordinates(-1.555,1.555,0.0), new Coordinates(-1.555,-1.555,0.0)
		};
	}
	
	private static Coordinates[] getTrigonalBipyramidal() {
		return new Coordinates[] {
				new Coordinates(0.0,0.0,0.0),
				new Coordinates(0,0,2.2), new Coordinates(0,0,-2.2), 
				new Coordinates(-1.1,1.905,0.0), new Coordinates(-1.1,-1.905,0.0),
				new Coordinates(2.2,0.0,0.0)
		};
	}


	@Override
	public Map<String, Double> getContributions() {
		Map<String,Double> contributions = new HashMap<String,Double>();
		double[] grad = new double[3*candidatePose.getLigConf().getMolecule().getAllAtoms()];
		double hbond = 0.0;
		for(PotentialEnergyTerm term : chemscoreHbond)
			hbond+=term.getFGValue(grad);
		contributions.put("HBOND", hbond);
		double metal = 0.0;
		for(PotentialEnergyTerm term : getMetalTerm(chemscoreMetal)) {
			metal+=term.getFGValue(grad);
		}
		contributions.put("METAL", metal);
		double plpContr = 0.0;
		for(PotentialEnergyTerm term : plp) 
			plpContr+=term.getFGValue(grad);
		contributions.put("PLP", plpContr);
		double strain = 0.0;
		Map<String, Object> ffOptions = new HashMap<String, Object>();
		ffOptions.put("dielectric constant", 80.0);

		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		StereoMolecule toOptimize = new StereoMolecule(candidatePose.getLigConf().getMolecule());
		toOptimize.ensureHelperArrays(Molecule.cHelperCIP);
		ForceFieldMMFF94 ff = new ForceFieldMMFF94(toOptimize, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
		double e = ff.getTotalEnergy();
		ff.minimise();
		double e0 = ff.getTotalEnergy();
		double deltaE = e-e0;
		if(deltaE>STRAIN_CUTOFF) {
			strain += deltaE-STRAIN_CUTOFF;
		}
		contributions.put("STRAIN", strain);
		return contributions;
	}






	
	
		

}
