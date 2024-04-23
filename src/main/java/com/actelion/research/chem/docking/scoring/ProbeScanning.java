package com.actelion.research.chem.docking.scoring;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.IntStream;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.docking.LigandPose;
import com.actelion.research.chem.docking.scoring.chemscore.HBTerm;
import com.actelion.research.chem.docking.scoring.chemscore.SimpleMetalTerm;
import com.actelion.research.chem.docking.scoring.plp.PLPTerm;
import com.actelion.research.chem.docking.scoring.plp.REPTerm;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import com.actelion.research.chem.io.pdb.converter.MoleculeGrid;
import com.actelion.research.chem.phesa.pharmacophore.pp.IPharmacophorePoint;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;

/**
 * scans a receptor binding site for interaction points using a probe atom
 * @author wahljo1
 *
 */

public class ProbeScanning {

	private IPharmacophorePoint.Functionality probeType;
	private Set<Integer> receptorAcceptors;	
	private Set<Integer> receptorDonorHs;
	private Map<Integer,Double> receptorDonorHPos;
	private Map<Integer,Double> receptorAcceptorNeg;
	private Set<Integer> receptorMetals;
	private Set<Integer> receptorDonors;
	
	private Set<Integer> ligandAcceptors;
	private Set<Integer> ligandDonors;
	private Map<Integer,Double> ligandDonorPos;
	private Map<Integer,Double> ligandAcceptorNeg;
	
	private StereoMolecule receptor;
	private Conformer receptorConf;
	private Probe probe;
	private List<PotentialEnergyTerm> plp;
	private List<PotentialEnergyTerm> chemscoreMetal;
	private List<HBProbeTerm> chemscoreHbond;
	private Set<Integer> bindingSiteAtoms;

	
	public ProbeScanning(StereoMolecule receptor, Set<Integer> bindingSiteAtoms, MoleculeGrid grid) {
		this.bindingSiteAtoms = bindingSiteAtoms;
		this.receptorConf = new Conformer(receptor);
		receptorAcceptors = new HashSet<>();	
		receptorDonorHs = new HashSet<>();
		receptorDonorHPos = new HashMap<Integer,Double>();
		receptorAcceptorNeg = new HashMap<Integer,Double>();
		receptorMetals = new HashSet<>();
		receptorDonors = new HashSet<>();
		this.receptor = receptor;
		ChemPLP.identifyHBondFunctionality(receptor,receptorAcceptors,receptorDonorHs, receptorDonors, receptorMetals,receptorAcceptorNeg,
				receptorDonorHPos);
	
	}
	
	public void init(Probe probe) {
		this.probe = probe;
		plp = new ArrayList<>();
		chemscoreHbond = new ArrayList<>();
		chemscoreMetal = new ArrayList<>();
		ligandAcceptors = new HashSet<>();
		ligandDonors = new HashSet<>();
		ligandDonorPos = new HashMap<Integer,Double>();
		ligandAcceptorNeg = new HashMap<Integer,Double>();
		IPharmacophorePoint.Functionality type = probe.type;
		switch(type) {
			case ACCEPTOR:
				ligandAcceptors.add(0);
				break;
			case DONOR:
				ligandDonors.add(0);
				break;
			case NEG_CHARGE:
				ligandAcceptorNeg.put(0,1.0);
				ligandAcceptors.add(0);
				break;
			case POS_CHARGE:
				ligandDonorPos.put(0,1.0);
				ligandDonors.add(0);
				break;
			default:
				break;
		}
		
		for(int p : bindingSiteAtoms) {
			if(receptor.getAtomicNo(p)==1) { // receptor hydrogen atom
				if(receptorDonorHs.contains(p)) {
					// receptor atom is donor hydrogen
					int d = receptor.getConnAtom(p, 0);
					boolean chargedP = receptorDonorHPos.keySet().contains(p);
					for(int l : ligandAcceptors) {
							int[] acceptorNeighbours = new int[0];
							boolean chargedL = ligandAcceptorNeg.keySet().contains(l);
							double scale = 1.0;
							if(chargedP && chargedL)
								scale = 2.0;
							HBProbeTerm hbTerm = HBProbeTerm.create(receptorConf, probe.probeConf, l, d, p, true, false, acceptorNeighbours, scale);
							chemscoreHbond.add(hbTerm);
						}	
				}
			}
			else { //receptor heavy atom
				if(receptorDonors.contains(p)) { // receptor donor heavy atom -> only plp terms 
					for(int l=0;l<probe.probeConf.getMolecule().getAtoms();l++) { //only consider ligand heavy atoms
						if(ligandAcceptors.contains(l)) {  //plp hbond donor-acceptor
							PLPTerm plpTerm = PLPTerm.create(receptorConf, probe.probeConf, p, l, PLPTerm.HBOND_TERM);
							plp.add(plpTerm);
						}
						else if(ligandDonors.contains(l)) {  //repulsive donor-donor
							REPTerm repTerm = REPTerm.create(receptorConf, probe.probeConf, p, l);
							plp.add(repTerm);
						}
						else { //buried
							PLPTerm plpTerm = PLPTerm.create(receptorConf,probe.probeConf, p, l, PLPTerm.BURIED_TERM);
							plp.add(plpTerm);
						}					
					}
				}
				else if(receptorAcceptors.contains(p)) { // receptor acceptor heavy atom
					int[] acceptorNeighbours = IntStream.range(0, receptor.getConnAtoms(p)).map(i -> receptor.getConnAtom(p, i)).toArray();
					boolean chargedP = receptorAcceptorNeg.keySet().contains(p);
					for(int l=0;l<probe.probeConf.getMolecule().getAtoms();l++) { 
							if(ligandDonors.contains(l)) {
								boolean chargedL = ligandDonorPos.keySet().contains(l);
								double scale = 1.0;
								if(chargedP && chargedL)
									scale = 2.0;
								HBProbeTerm hbTerm = HBProbeTerm.create(receptorConf, probe.probeConf, p, l,-1, false, true, acceptorNeighbours, scale);
								chemscoreHbond.add(hbTerm);
								PLPTerm plpTerm = PLPTerm.create(receptorConf, probe.probeConf, p, l, PLPTerm.HBOND_TERM);
								plp.add(plpTerm);
							}
													
							else if(ligandAcceptors.contains(l)) {  //repulsive acceptor-acceptor
								REPTerm repTerm = REPTerm.create(receptorConf, probe.probeConf, p, l);
								plp.add(repTerm);
							}
							else { //buried
								PLPTerm plpTerm = PLPTerm.create(receptorConf, probe.probeConf, p, l, PLPTerm.BURIED_TERM);
								plp.add(plpTerm);
							}			
							
						
					}	
					
				}
				else if(receptorMetals.contains(p)) {
						for(int l=0;l<probe.probeConf.getMolecule().getAtoms();l++) { 
							if(ligandDonors.contains(l)) {  //met-donor -> repulsive
								REPTerm repTerm = REPTerm.create(receptorConf, probe.probeConf, p, l);
								plp.add(repTerm);
							}
							else if(ligandAcceptors.contains(l)) { //attractive met-acc interaction
								PLPTerm plpTerm = PLPTerm.create(receptorConf, probe.probeConf, p, l, PLPTerm.METAL_TERM);
								plp.add(plpTerm);
							}
							else { //buried met-nonp interaction
								PLPTerm plpTerm = PLPTerm.create(receptorConf, probe.probeConf, p, l, PLPTerm.BURIED_TERM);
								plp.add(plpTerm);
							}	
						}
					
					//if(SIMPLE_METAL_ATOMS.contains(receptor.getAtomicNo(p))) {
						for(int l : ligandAcceptors) {
							double scale = 1.0;
							if(ligandAcceptorNeg.keySet().contains(l))
								scale = 2.0;
							int[] acceptorNeighbours = new int[0];
							SimpleMetalTerm metTerm = SimpleMetalTerm.create(receptorConf, probe.probeConf, 
									l, p, acceptorNeighbours, scale);
							chemscoreMetal.add(metTerm);
						}
	
						//}
						/*
					else { //standard metal term;
						List<Coordinates> interactionSites = metalInteractionSites.get(p);
						for(int l : ligandAcceptors) {
							double scale = 1.0;
							if(ligandAcceptorNeg.contains(l))
								scale = 2.0;
							int[] acceptorNeighbours = IntStream.range(0, ligand.getConnAtoms(l)).map(i -> ligand.getConnAtom(l, i)).toArray();
							for(Coordinates site : interactionSites) {
								MetalTerm metTerm = MetalTerm.create(candidatePose, l, acceptorNeighbours, site,scale);
								chemscoreMetal.add(metTerm);
							}
					}
					
					
					}
					*/
				}
				else { // non-polar heavy atom
					for(int l=0;l<probe.probeConf.getMolecule().getAtoms();l++) { 
						if(ligandDonors.contains(l)) {  //nonpolar-donor -> buried
							PLPTerm plpTerm = PLPTerm.create(receptorConf, probe.probeConf, p, l, PLPTerm.BURIED_TERM);
							plp.add(plpTerm);
						}
						else if(ligandAcceptors.contains(l)) { //nonpolar-acceptor -> buried
							PLPTerm plpTerm = PLPTerm.create(receptorConf, probe.probeConf, p, l, PLPTerm.BURIED_TERM);
							plp.add(plpTerm);
						}
						else { // nonpolar-nonpolar -> nonpolar
							PLPTerm plpTerm = PLPTerm.create(receptorConf, probe.probeConf, p, l, PLPTerm.NONPOLAR_TERM);
							plp.add(plpTerm);
						}	
				}
				}
			}
			}
		
		
		
		
	}
	
	public double getScore() {
	
		double energy = 0.0;
		double[] grad = new double[3];
		for(HBProbeTerm term : chemscoreHbond)
			energy+=term.getEnergy();
		for(PotentialEnergyTerm term : chemscoreMetal) 
			energy+=term.getFGValue(grad);
		for(PotentialEnergyTerm term : plp) 
			energy+=term.getFGValue(grad);
		return energy;
	}
	

	
	public static class Probe {
		Conformer probeConf;
		IPharmacophorePoint.Functionality type;
		
		public Probe(Coordinates c, IPharmacophorePoint.Functionality type) {
			StereoMolecule probeMol = new StereoMolecule();
			probeMol.addAtom(6);
			probeMol.ensureHelperArrays(Molecule.cHelperNeighbours);
			probeConf = new Conformer(probeMol);
			probeConf.setCoordinates(0, c);
			this.type = type;
		}
		
		public void updateCoordinates(Coordinates c) {
			probeConf.setCoordinates(0, c);
		}
	}
	
	/*
	 * based on ChemScore HBTerm, but adapted for probe scanning
	 */
	private static class HBProbeTerm {
		
		private static final double D0 = 1.85;
		private static final double D1 = 0.25;
		private static final double D2 = 0.65;

		
		private static final double PHI0 = Math.PI;
		private static final double PHI1 = 30*Math.PI/180.0;
		private static final double PHI2 = 80*Math.PI/180.0;

		
		private static final double PSI0 = Math.PI;
		private static final double PSI1 = 80*Math.PI/180.0;
		private static final double PSI2 = 100*Math.PI/180.0;
		
		private static final double ENERGY = -3.0;
		
		private Conformer receptor;
		private Conformer probe;
		private int acceptor;
		private int donor;
		private int hydrogen;
		private double scale;
		private boolean isProbeAcceptor;
		private boolean isProbeDonor;
		private int[] acceptorNeighbours;

		
		
		private HBProbeTerm(Conformer receptor, Conformer probe, int acceptor, int donor, int hydrogen,
				boolean isProbeAcceptor, boolean isProbeDonor, int[] acceptorNeighbours, double scale) {
			this.receptor = receptor;
			this.probe = probe;
			this.acceptor = acceptor;
			this.donor = donor;
			this.hydrogen = hydrogen;
			this.scale = scale;
			this.isProbeAcceptor = isProbeAcceptor;
			this.isProbeDonor = isProbeDonor;
			this.acceptorNeighbours = acceptorNeighbours;
			assert(isProbeAcceptor!=isProbeDonor);
		}
		
		
		public static HBProbeTerm create(Conformer receptor, Conformer ligand, int acceptor, int donor, int hydrogen, 
				boolean isProbeAcceptor, boolean isProbeDonor, int[] acceptorNeighbours, double scale) {
			return new HBProbeTerm(receptor, ligand, acceptor, donor, hydrogen, isProbeAcceptor, isProbeDonor,
					acceptorNeighbours, scale);
		}
		
		
		/** d0 = 1.85;
		 *  d1 = 0.25;
		 *  d2 = 0.65
		 *  x = r-d0     r is the hydrogen-acceptor distance
		 *  if(|x|<d1 f = 1;
		 *  if(|x|>d2 f = 0;
		 *  else f(x) = (x2-x)/(x2-x1)
		 * @param gradient
		 * @return
		 */
		private double getDistTerm() {
			Coordinates a;
			Coordinates h;
			double energy = 0.0;
			if(isProbeAcceptor)
				a = probe.getCoordinates(acceptor);
			else 
				a = receptor.getCoordinates(acceptor);
			if(isProbeDonor) {
				//if the probe is a donor functionality, the dummy hydrogen is optimally placed along the A--D vector, so that the 
				// angle A--H--D equals 180 degrees
				Coordinates d = probe.getCoordinates(donor);
				Coordinates v = a.subC(d);
				double scale = 1.0/v.dist();
				h = d.addC(v.scaleC(scale));
				assert(h.subC(d).dist()-1.0<0.001);
			}
			else {
				h = receptor.getCoordinates(hydrogen);
			}
			Coordinates r = a.subC(h);
			double d = r.dist();
			double distTerm = d - D0;
			if(distTerm<0)
				distTerm = -distTerm;
			if(distTerm<D1) {
				energy = 1.0;
			}
			else if(distTerm>D2) {
				energy = 0.0;
			}
			else {
				energy = (D2-distTerm)/(D2-D1);
			}

					
			return energy;
			
		}
		
		private double getAngleTerm(int a1, int a2, int a3, boolean isLigAtomA1, boolean
				isLigAtomA2, boolean isLigAtomA3, double x0, double x1, double x2) {
			
			//if two atoms come from the receptor, the angles are calculated, else, perfect geometry is assumed 
			// (equal to optimally placing the probe)
			
			boolean[] b = new boolean[] {isLigAtomA1, isLigAtomA2, isLigAtomA3};
			int sum = 0;
			for(boolean bool : b) {
				if(bool)
					sum++;
			}
			
			double energy = 0.0;
			if(sum==2)
				energy = 1.0;
			else {
				Coordinates c1, c2, c3;
				if(isLigAtomA1)
					c1 = probe.getCoordinates(a1);
				else
					c1 = receptor.getCoordinates(a1);
				if(isLigAtomA2)
					c2 = probe.getCoordinates(a2);
				else
					c2 = receptor.getCoordinates(a2);
				if(isLigAtomA3)
					c3 = probe.getCoordinates(a3);
				else
					c3 = receptor.getCoordinates(a3);
				
			    Coordinates r0 = c1.subC(c2).unit();
			    Coordinates r1 = c3.subC(c2).unit();

			    double cosTheta = r0.cosAngle(r1);
	
			    double angleTerm = Math.acos(cosTheta) - x0;
			    if(angleTerm<0) {
			    	angleTerm=-angleTerm;
			    }
			    	
			    
			    if(angleTerm<x1)
			    	energy = 1.0;
			    else if(angleTerm>x2)
			    	energy = 0.0;
			    else {
				    energy = (x2-angleTerm)/(x2-x1);
	
			    }
			}
			return energy;
			
		}
		


		public double getEnergy() {
			List<Double> energies = new ArrayList<>();
			double energy = 0.0;
			energy = getDistTerm();
			if(energy!=0.0) {
				energies.add(energy);
				//D--H--A angle term
				energy = getAngleTerm(donor,hydrogen,acceptor,
						isProbeDonor,isProbeDonor,isProbeAcceptor,
						PHI0,PHI1,PHI2);
				energies.add(energy);

				if(isProbeAcceptor) {
					energy=getAngleTerm(-1,acceptor,hydrogen,isProbeAcceptor,isProbeAcceptor,
							isProbeDonor,PSI0,PSI1,PSI2);
					energies.add(energy);
				}
				else {
					for(int aa : acceptorNeighbours) {
						energy=getAngleTerm(aa,acceptor,donor,isProbeAcceptor,isProbeAcceptor,
								isProbeDonor,PSI0,PSI1,PSI2);
						energies.add(energy);
					}
					
				}
			}
			else 
				energies.add(0.0);
			double totEnergy = scale*ENERGY;
			for(double eng : energies)
				totEnergy*=eng;
			return totEnergy;
		}

	}
	

	

	
	
}