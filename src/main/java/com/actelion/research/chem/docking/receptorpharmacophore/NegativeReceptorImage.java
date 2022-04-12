package com.actelion.research.chem.docking.receptorpharmacophore;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.docking.DockingUtils;
import com.actelion.research.chem.docking.scoring.ProbeScanning;
import com.actelion.research.chem.io.pdb.converter.MoleculeGrid;
import com.actelion.research.chem.phesa.AtomicGaussian;
import com.actelion.research.chem.phesa.Gaussian3D;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.ShapeVolume;
import com.actelion.research.chem.phesa.pharmacophore.PharmacophoreCalculator;
import com.actelion.research.chem.phesa.pharmacophore.pp.IPharmacophorePoint;
import com.actelion.research.chem.phesa.pharmacophore.pp.PPGaussian;
import com.actelion.research.chem.phesa.pharmacophore.pp.SimplePharmacophorePoint;
import smile.clustering.KMeans;

import java.util.*;
import java.util.stream.IntStream;

/**
 * creates a negative receptor image in a binding site, using atomic gaussians for shape and pharmacophore gaussians to
 * indicate optimal sites of interaction for a ligand
 * - probe atoms are used around receptor pharmacophore features to localize optimal sites of interaction
 * - interaction sites are clustered  
 * @author wahljo1
 *
 */


public class NegativeReceptorImage extends MoleculeGrid {
	
	public enum InteractionProbe {NEG_CHARGE, POS_CHARGE, HB_DONOR, HB_ACCEPTOR };
	
	private static final long SEED = 12345L;
	private static final int STARTING_POINTS_CAVITY_DETECTION = 10;
	private static final int RAYS = 120;
	private static final double RAY_LENGTH = 8.0;
	private static final double BURIEDNESS_RATIO_CUTOFF = 0.4;
	private static final int BUMP_RADIUS = 3;
	private static final double BUMP_RADIUS2 = 2.0; //to check buriedness
	private static final double INTERACTION_CUTOFF = 0.0;
	private static final double INTERACTION_CUTOFF_CHARGE = -1.0;
	private static final double NONPOLAR_RADIUS = 2.5;
	private static final double GAUSSIAN_DISTANCE = 1.3;
	private static final int MAX_NR_INTERACTION_POINTS = 8; //maximum number of interaction points per site
	private static double MIN_INTERACTION_POINT_DIST = 1.5; //minimum distance between polar interaction points

	private boolean[][][] bumpGrid;	
	private Set<Integer> receptorAtoms;
	private StereoMolecule receptor;
	private Map<Integer,Map<Integer,List<Coordinates>>> receptorInteractionSites;
	private ProbeScanning probeScanning;
	private MoleculeGrid receptorGrid;

	public NegativeReceptorImage(StereoMolecule ligand, StereoMolecule receptor) {
		this(ligand, receptor,0.4,new Coordinates(5.0,5.0,5.0));
	}
	

	public NegativeReceptorImage(StereoMolecule ligand, StereoMolecule receptor,double gridWidth, Coordinates extension) {
		super(ligand,gridWidth,extension);
		receptorGrid = new MoleculeGrid(receptor,gridWidth,new Coordinates(0.0,0.0,0.0));
		this.receptor = receptor;
		receptorAtoms = new HashSet<Integer>();
		receptorInteractionSites = new HashMap<Integer,Map<Integer,List<Coordinates>>>();
		int[] recGridSize = receptorGrid.getGridSize();
		bumpGrid = new boolean[recGridSize[0]][recGridSize[1]][recGridSize[2]];
		for(int i=0;i<receptor.getAllAtoms();i++) {
			int[] gridC = getGridCoordinates(receptor.getCoordinates(i));
			int x = gridC[0];
			int y = gridC[1];
			int z = gridC[2];	
			if(x>0 && x<grid.length) {
				if(y>0 && y<grid[0].length) {
					if(z>0 && z<grid[0][0].length) {
						receptorAtoms.add(i);
					}
			}
		}
		}
		probeScanning = new ProbeScanning(receptor, receptorAtoms, this);
	}
	
	public ShapeVolume calculate() {
		List<PPGaussian> ppGaussians = new ArrayList<>();
		List<AtomicGaussian> shapeGaussians = new ArrayList<>();
		analyzeBindingSiteAtoms();
		ShapeVolume recVol = new ShapeVolume();
		analyzeBumps();
		createPolarInteractionSites(ppGaussians);
		createShapeAtoms(shapeGaussians);
		List<Coordinates> startingPoints = new ArrayList<>();
		for(int a=0;a<mol.getAtoms() && a<STARTING_POINTS_CAVITY_DETECTION;a++) {
			startingPoints.add(mol.getCoordinates(a));
		}
		prunePoints(startingPoints, ppGaussians, shapeGaussians, 2.0);
		//assign ids to SimplePPs
		for(int i=0;i<ppGaussians.size();i++) {
			ppGaussians.get(i).setAtomId(i);
		}
		recVol.setAtomicGaussians(shapeGaussians);
		recVol.setPPGaussians(ppGaussians);
		return recVol;
	}
	
	private void createPolarInteractionSites(List<PPGaussian> ppGaussians) {
		for(int i : receptorInteractionSites.keySet()) {
			for(int type : receptorInteractionSites.get(i).keySet()) {
			for(Coordinates c : receptorInteractionSites.get(i).get(type)) {
				IPharmacophorePoint.Functionality functionality = null;
				double energyCutoff = INTERACTION_CUTOFF;
				int coarseFactor = 1;
				switch(type) {
					case PharmacophoreCalculator.CHARGE_NEG_ID:
						functionality = IPharmacophorePoint.Functionality.POS_CHARGE;
						energyCutoff = INTERACTION_CUTOFF_CHARGE;
						coarseFactor = 2;
						break;
					case PharmacophoreCalculator.CHARGE_POS_ID:
						functionality = IPharmacophorePoint.Functionality.NEG_CHARGE;
						energyCutoff = INTERACTION_CUTOFF_CHARGE;
						coarseFactor = 2;
						break;
					case PharmacophoreCalculator.ACCEPTOR_ID:
						functionality = IPharmacophorePoint.Functionality.DONOR;
						break;
					case PharmacophoreCalculator.DONOR_ID:
						functionality = IPharmacophorePoint.Functionality.ACCEPTOR;
						break;
					default:
						functionality = null;
				}
				if(functionality==null)
					continue;
				double[][] centroids = scanProbe(c,functionality,i, energyCutoff, coarseFactor);
				
				if(centroids!=null) {
					for(double[] centroid : centroids) {
						SimplePharmacophorePoint spp = new SimplePharmacophorePoint(-1,new Coordinates(centroid[0],
								centroid[1], centroid[2]), functionality);
						PPGaussian ppg = new PPGaussian(6, spp);
						ppGaussians.add(ppg);
					}
				}
				}
			}
		}
		

	}
	/**
	 * fill accessible parts of the binding site with shape gaussian, without incorporating surface points (check buriedness)
	 * @param shapeGaussians
	 */
	private void createShapeAtoms(List<AtomicGaussian> shapeGaussians) {
		List<AtomicGaussian> gaussians = new ArrayList<>();
		double radiusSq = NONPOLAR_RADIUS*NONPOLAR_RADIUS;
		double gaussDistSq = GAUSSIAN_DISTANCE*GAUSSIAN_DISTANCE;
		int cutoff = (int) (NONPOLAR_RADIUS / this.gridWidth);
		for(int x=cutoff;x<this.gridSize[0];x++) {
			for(int y=cutoff;y<this.gridSize[1];y++) {
				for(int z=cutoff;z<this.gridSize[2];z++) {
					boolean clash = false;
					Coordinates probeCoords = this.getCartCoordinates(new int[] {x,y,z});
					for(AtomicGaussian ag : gaussians) {
						double dx = (ag.getCenter().x-probeCoords.x);
						double dy = (ag.getCenter().y-probeCoords.y);
						double dz = (ag.getCenter().z-probeCoords.z);
						double rSq = dx*dx + dy*dy + dz*dz;
						if(rSq>gaussDistSq ) {
							continue;
						}
						else {
							double r = Math.sqrt(rSq);
							if(r<GAUSSIAN_DISTANCE) { // too close to another automic gaussian
								clash = true;
								break;
							}
						}
					}
					if(!clash) {
						for(int atom : receptorAtoms) {
							Coordinates c = receptor.getCoordinates(atom);
							double dx = (c.x-probeCoords.x);
							double dy = (c.y-probeCoords.y);
							double dz = (c.z-probeCoords.z);
							double rSq = dx*dx + dy*dy + dz*dz;
							if(rSq>radiusSq) {
								continue;
							}
							else {
								double r = Math.sqrt(rSq);
								if(r<NONPOLAR_RADIUS) {
									clash = true;
									break;
								}
							}
						}
					}					
					boolean isBuried = getBuriedness(new int[] {x,y,z});
					if(!clash && isBuried) {
						AtomicGaussian ag = new AtomicGaussian(-1,6,probeCoords);
						gaussians.add(ag);
					}
				}		 
			}
		}
		gaussians.stream().forEach(ag -> {
			shapeGaussians.add(ag);
		});
	}
	
	
	
	/**
	 * using the given probe, the vicinity of the given receptor atom is scanned for interaction
	 * hot spots
	 * @param probeType
	 * @param receptorAtom
	 */
	private double[][] scanProbe(Coordinates c, IPharmacophorePoint.Functionality functionality, int receptorAtom, double energyCutoff, 
			int coarseFactor) {
		
		ProbeScanning.Probe probe = new ProbeScanning.Probe(new Coordinates(),functionality);
		probeScanning.init(probe);
		List<double[]> interactionPoints = new ArrayList<>(); // input for clustering procedure
		List<Double> interactionEnergies = new ArrayList<>();
		int bound = coarseFactor;
		int[] probeCenterCoords = this.getGridCoordinates(c);
		int cutoff = (int) (BUMP_RADIUS/this.gridWidth);
		for(int dx=-bound;dx<bound+1;dx++) {
			for(int dy=-bound;dy<bound+1;dy++) {
				for(int dz=-bound;dz<bound+1;dz++) {
					int x = probeCenterCoords[0]+dx*coarseFactor;
					int y = probeCenterCoords[1]+dy*coarseFactor;
					int z = probeCenterCoords[2]+dz*coarseFactor;
					if(x<cutoff || x>(this.gridSize[0]-cutoff))
						continue;
					if(y<cutoff || y>(this.gridSize[1]-cutoff))
						continue;
					if(z<cutoff || z>(this.gridSize[2]-cutoff))
						continue;
					boolean isBuried = getBuriedness(new int[] {x,y,z});
					if(!isBuried)
						continue;
					Coordinates probeCoords = this.getCartCoordinates(new int[] {x,y,z});
					probe.updateCoordinates(probeCoords);
					double score = probeScanning.getScore();
					if(score<energyCutoff) {
						interactionPoints.add(new double[] {probeCoords.x,probeCoords.y,probeCoords.z});
						interactionEnergies.add(score);
					}
					
					
					
				}
			}
		}
		/*
		double[][] centroids = new double[interactionPoints.size()][3];
		for(int i=0;i<interactionPoints.size();i++) {
			centroids[i] = interactionPoints.get(i);
		}
		*/
		
		double[][] centroids;
		if(interactionPoints.size()==0)
			centroids = null;
		else {
			double[][] data = new double[interactionPoints.size()][3];
			IntStream.range(0, interactionPoints.size()).forEach(i -> data[i] = interactionPoints.get(i));
			final double minCentroidDistSq = MIN_INTERACTION_POINT_DIST*MIN_INTERACTION_POINT_DIST;
			int maxNrPoints = Math.min(data.length, MAX_NR_INTERACTION_POINTS);
			int k = maxNrPoints;
			if(maxNrPoints==1) {
				centroids = new double[][] {{data[0][0],data[0][1],data[0][2]}};
			}
			else {	
				boolean breakCondition = false;
				for(int i=2;i<maxNrPoints+1 && !breakCondition;i++) {
					KMeans kmeans = new KMeans(data,i);
					centroids = kmeans.centroids();
					for(int j=0;j<centroids.length && !breakCondition;j++) {
						double[] centroid1 = centroids[j];
						for(int l=j+1;l<centroids.length && !breakCondition;l++) {
							double[] centroid2 = centroids[l];
							double dx = centroid1[0] - centroid2[0];
							double dy = centroid1[1] - centroid2[1];
							double dz = centroid1[2] - centroid2[2];
							double distSq = dx*dx + dy*dy + dz*dz;
							if(distSq>0.001 && distSq < minCentroidDistSq) { //clusters too close
								k = i-1;
								breakCondition = true;
							}
						}
					}
				}
				if(k==1) { //only one interaction point, take geometric mean
					double[] centroid = new double[] {0.0,0.0,0.0};
					int n = 0;
					for(double[] p : data) {
						centroid[0]+=p[0];
						centroid[1]+=p[1];
						centroid[2]+=p[2];
						n+=1;
					}
					centroid[0]/=n;
					centroid[1]/=n;
					centroid[2]/=n;
					centroids = new double[][] {centroid};
				}
				else {
					KMeans kmeans = new KMeans(data,k);
					centroids = kmeans.centroids();
				}
			}
		}
		
		return centroids;
		
	}
	
	private void analyzeBumps() {
		double radiusSq = BUMP_RADIUS2*BUMP_RADIUS2;

// Replaced because of out-of-bounds exceptions; TLS 17Oct2021
//		for(int x=0;x<this.gridSize[0];x++) {
//			for(int y=0;y<this.gridSize[1];y++) {
//				for(int z=0;z<this.gridSize[2];z++) {

		int xmax = Math.min(this.gridSize[0], bumpGrid.length);
		int ymax = Math.min(this.gridSize[1], bumpGrid[0].length);
		int zmax = Math.min(this.gridSize[2], bumpGrid[0][0].length);
		for(int x=0;x<xmax;x++) {
			for(int y=0;y<ymax;y++) {
				for(int z=0;z<zmax;z++) {

					Coordinates probeCoords = this.getCartCoordinates(new int[] {x,y,z});
					for(int atom=0;atom<receptor.getAtoms();atom++) {
						Coordinates c = receptor.getCoordinates(atom);
						double dx = (c.x-probeCoords.x);
						double dy = (c.y-probeCoords.y);
						double dz = (c.z-probeCoords.z);
						double rSq = dx*dx + dy*dy + dz*dz;
						if(rSq>radiusSq) {
							continue;
						}
						else {
							double r = Math.sqrt(rSq);
							if(r<BUMP_RADIUS2) {
								bumpGrid[x][y][z] = true;
								break;
							}
						}
					}
				}
			}
		}

	}
	
	public boolean getBuriedness(int[] gridCoords) {
		Random rnd = new Random(SEED);
		int steps = (int) (RAY_LENGTH/this.gridWidth);
		int intersections = 0;
		Coordinates center = this.getCartCoordinates(gridCoords);
		for(int i=0;i<RAYS;i++) {
			Coordinates ray = DockingUtils.randomVectorInSphere(rnd).scale(this.gridWidth);
			for(int j=1;j<steps+1;j++) {
				Coordinates c = center.addC(ray.scaleC(j));
				int[] gridC = this.getGridCoordinates(c);
				try {
					boolean intersect = bumpGrid[gridC[0]][gridC[1]][gridC[2]];
					if(intersect) {
						intersections++;
						break;
					}
				}
				catch(Exception e) {
					continue;
				}

				}
			}
		double buriedness = (double)intersections/RAYS;
		boolean isBuried = buriedness>BURIEDNESS_RATIO_CUTOFF ? true : false;
		return isBuried;
	}
	
	
	/**
	 * go through list of receptor pharmacophore features and calculate preferable sites of interaction
	 */
	private void analyzeBindingSiteAtoms() {
		MolecularVolume molVol = new MolecularVolume(receptor);
		final double hbondDistance = 2.8;
		for(PPGaussian ppg : molVol.getPPGaussians()) {
			int a = ppg.getAtomId();
			IPharmacophorePoint pp = ppg.getPharmacophorePoint();
			int funcIndex = pp.getFunctionalityIndex();
			receptorInteractionSites.putIfAbsent(a, new HashMap<Integer,List<Coordinates>>());
			Map<Integer,List<Coordinates>> sites = receptorInteractionSites.get(a);
			sites.putIfAbsent(funcIndex, new ArrayList<>());
			List<Coordinates> siteCoords = sites.get(funcIndex);
			if(ppg.getPharmacophorePoint().getFunctionalityIndex()==IPharmacophorePoint.Functionality.ACCEPTOR.getIndex()) {
				if(receptor.getConnAtoms(a)==1) {
					if(receptor.getBondOrder(receptor.getBond(a, receptor.getConnAtom(a,0)))==2) {//sp2 oxygen
						if(siteCoords.size()==0) { //add interaction point along C=0 axis
							Coordinates v = receptor.getCoordinates(a).subC(receptor.getCoordinates( receptor.getConnAtom(a,0))).unitC();
							Coordinates c =  receptor.getCoordinates(a).addC(v.scaleC(hbondDistance));
							siteCoords.add(c);
						}
					}
				}
				Coordinates c = receptor.getCoordinates(a).addC(pp.getDirectionality().scaleC(hbondDistance )); // interaction points at lone pairs
				siteCoords.add(c);
			}
			else if(ppg.getPharmacophorePoint().getFunctionalityIndex()==IPharmacophorePoint.Functionality.NEG_CHARGE.getIndex()) {
				//either -O(-) or -C(=O)-O(-)
				if(receptor.getConnAtoms(a)==1) {
					Coordinates v = receptor.getCoordinates(a).subC(receptor.getCoordinates( receptor.getConnAtom(a,0))).unitC();
					Coordinates c =  receptor.getCoordinates(a).addC(v.scaleC(hbondDistance ));
					siteCoords.add(c);
				}
				else if(receptor.getConnAtoms(a)==3) {
					int aa = -1;
					for(int i=0;i<receptor.getConnAtoms(a);i++) {
						if(receptor.getAtomicNo(receptor.getConnAtom(a, i))==6) {
							aa = receptor.getConnAtom(a, i);
							break;
						}
					}
					if(aa!=-1) {
						Coordinates v = receptor.getCoordinates(a).subC(receptor.getCoordinates(aa)).unitC();
						Coordinates c =  receptor.getCoordinates(a).addC(v.scaleC(4.0));
						siteCoords.add(c);
	
					}
				}
			}
			
			else if(ppg.getPharmacophorePoint().getFunctionalityIndex()==IPharmacophorePoint.Functionality.POS_CHARGE.getIndex()) {
				//either -O(-) or -C(=O)-O(-)
				if(receptor.getConnAtoms(a)==1) {
					Coordinates v = receptor.getCoordinates(a).subC(receptor.getCoordinates( receptor.getConnAtom(a,0))).unitC();
					Coordinates c =  receptor.getCoordinates(a).addC(v.scaleC(4.0));
					siteCoords.add(c);
				}
				else if(receptor.getConnAtoms(a)==3) { //guanidinium
					int aa = -1;
					for(int i=0;i<receptor.getConnAtoms(a);i++) {
						int neighbour = receptor.getConnAtom(a, i);
						if(receptor.getConnAtoms(neighbour)==2) {
							aa=neighbour;
						}
					}
					if(aa!=-1) {
						Coordinates v = receptor.getCoordinates(a).subC(receptor.getCoordinates(aa)).unitC();
						Coordinates c =  receptor.getCoordinates(a).addC(v.scaleC(4.0));
						siteCoords.add(c);
					}
				}
			}
			else if(ppg.getPharmacophorePoint().getFunctionalityIndex()==IPharmacophorePoint.Functionality.DONOR.getIndex()) { 
				Coordinates c = receptor.getCoordinates(a).addC(pp.getDirectionality().scaleC(hbondDistance-1.0));
				siteCoords.add(c);
			}
				
		}
		
	}
	
	/**
	 * algorithm to remove disconnected clusters of Shape Gaussians and PP Gaussians and only keep one single cavity
	 * @param center
	 */
	private void prunePoints(List<Coordinates> centers, List<PPGaussian> ppGaussians,  List<AtomicGaussian> atomicGaussians, double stepSize) {
		Set<Gaussian3D> cavityPoints = new HashSet<Gaussian3D>();
		Queue<Gaussian3D> pq = new LinkedList<Gaussian3D>();
		Set<Gaussian3D> visited = new HashSet<Gaussian3D>();
		Set<Gaussian3D> neighbours = new HashSet<>();
		for(Coordinates center:centers) {
			neighbours.addAll(getNeighbourPoints(center, ppGaussians, atomicGaussians, stepSize));
		}
	
		pq.addAll(neighbours);
		cavityPoints.addAll(neighbours);
		while(!pq.isEmpty()) {
			Gaussian3D point = pq.poll();
			if(visited.contains(point))
				continue;
			visited.add(point);
			neighbours = getNeighbourPoints(point.getCenter(), ppGaussians, atomicGaussians, stepSize);
			pq.addAll(neighbours);
			cavityPoints.addAll(neighbours);
		}
		List<PPGaussian> ppgToDelete = new ArrayList<PPGaussian>();
		for(PPGaussian ppg:ppGaussians) {
			if(!cavityPoints.contains(ppg))
				ppgToDelete.add(ppg);
		}
		ppGaussians.removeAll(ppgToDelete);
		
		List<AtomicGaussian> agToDelete = new ArrayList<AtomicGaussian>();
		for(AtomicGaussian ag:atomicGaussians) {
			if(!cavityPoints.contains(ag))
				agToDelete.add(ag);
		}
		atomicGaussians.removeAll(agToDelete);
		
	}
	
	private Set<Gaussian3D> getNeighbourPoints(Coordinates center, List<PPGaussian> ppGaussians,  List<AtomicGaussian> atomicGaussians, double distCutoff) {
		double distSqCutoff = distCutoff*distCutoff;
		Set<Gaussian3D> neighbours = new HashSet<Gaussian3D>();
		Set<Gaussian3D> allGauss = new HashSet<Gaussian3D>();
		allGauss.addAll(atomicGaussians);
		allGauss.addAll(ppGaussians);
		for(Gaussian3D gauss2 : allGauss) {
			double distSq = gauss2.getCenter().distSquareTo(center);
			if(distSq>distSqCutoff)
				continue;
			else 
				neighbours.add(gauss2);
		}
		return neighbours;
		
	}
	
	

	

}
