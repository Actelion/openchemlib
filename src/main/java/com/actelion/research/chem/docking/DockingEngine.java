package com.actelion.research.chem.docking;


import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.*;
import com.actelion.research.chem.alignment3d.KabschAlignment;
import com.actelion.research.chem.alignment3d.transformation.Rotation;
import com.actelion.research.chem.alignment3d.transformation.TransformationSequence;
import com.actelion.research.chem.alignment3d.transformation.Translation;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.conf.ConformerSetGenerator;
import com.actelion.research.chem.docking.receptorpharmacophore.NegativeReceptorImageCreator;
import com.actelion.research.chem.docking.scoring.AbstractScoringEngine;
import com.actelion.research.chem.docking.scoring.ChemPLP;
import com.actelion.research.chem.docking.scoring.IdoScore;
import com.actelion.research.chem.docking.shape.ShapeDocking;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import com.actelion.research.chem.forcefield.mmff.MMFFPositionConstraint;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import com.actelion.research.chem.io.pdb.converter.MoleculeGrid;
import com.actelion.research.chem.mcs.MCS;
import com.actelion.research.chem.optimization.OptimizerLBFGS;
import com.actelion.research.chem.phesa.EncodeFunctions;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.PheSAAlignment;
import com.actelion.research.chem.phesa.ShapeVolume;
import com.actelion.research.chem.potentialenergy.PositionConstraint;
import org.openmolecules.chem.conf.gen.ConformerGenerator;

import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.Base64.Decoder;
import java.util.Base64.Encoder;
import java.util.Map.Entry;
import java.util.stream.Collectors;

public class DockingEngine {
/**
 * applies molecular docking to find the binding pose of a ligand molecule into the binding site of the protein
 * nativeLigand: defines the location of the binding site
 * this class is not thread safe! every thread requires it's own instance
 * @author wahljo1
 *
 */
	
	public enum ScoringFunction {CHEMPLP,IDOSCORE;}
	private static final int DEFAULT_NR_MC_STEPS = 50;
	private static final int DEFAULT_START_POSITIONS = 10;
	private static final double BOLTZMANN_FACTOR = 1.2; //as for AutoDock Vina
	public static final double GRID_DIMENSION = 6.0;
	public static final double GRID_RESOLUTION = 0.5;
	public static final double MINI_CUTOFF = 100; //if energy higher after MC step, don't minimize
	public static final int MCS_EXHAUSTIVENESS = 3; //MCS docking requires more starting positions, but less MC sampling (reduced DOF)
	

	private Rotation rotation; //for initial prealignment to native ligand
	private Coordinates origCOM;
	private int mcSteps;
	private Random random;
	private AbstractScoringEngine engine;
	private int startPositions;
	private StereoMolecule nativeLigand;
	private ShapeDocking shapeDocking;
	private ThreadMaster threadMaster;
	private StereoMolecule mcsRef;
	private List<Integer> mcsConstrainedBonds;
	private List<Integer> mcsConstrainedAtoms;
	private double gridDimension;

	
	public DockingEngine(StereoMolecule rec, StereoMolecule nativeLig, int mcSteps, int startPositions, double gridDimension,
			ScoringFunction scoringFunction) throws DockingFailedException {
		for(int ra=0;ra<rec.getAtoms();ra++) {
			if(rec.getImplicitHydrogens(ra)>0)
				throw new DockingFailedException("please add hydrogen atoms to receptor structure!");
		}
		nativeLigand = new Molecule3D(nativeLig);
		nativeLigand.ensureHelperArrays(Molecule.cHelperCIP);
		Molecule3D receptor = new Molecule3D(rec);
		receptor.ensureHelperArrays(Molecule.cHelperCIP);
		MolecularVolume molVol = new MolecularVolume(nativeLigand);
		origCOM  = new Coordinates(molVol.getCOM());
		Conformer conf = new Conformer(nativeLigand);
		rotation = molVol.preProcess(conf);
		this.startPositions = startPositions;
		preprocess(receptor,nativeLigand);
		
		TransformationSequence transform = new TransformationSequence();
		ShapeVolume bsVolume = NegativeReceptorImageCreator.create(nativeLigand, receptor, transform);
		shapeDocking = new ShapeDocking(bsVolume,transform);
		
		MoleculeGrid grid = new MoleculeGrid(nativeLigand,GRID_RESOLUTION,
				new Coordinates(gridDimension,gridDimension,
						gridDimension));
		
		Set<Integer> bindingSiteAtoms = new HashSet<Integer>();
		
		if(scoringFunction==ScoringFunction.CHEMPLP) {
			DockingEngine.getBindingSiteAtoms(receptor, bindingSiteAtoms, grid, true);
			engine = new ChemPLP(receptor,bindingSiteAtoms,grid);
		}
		else if(scoringFunction==ScoringFunction.IDOSCORE) {
			DockingEngine.getBindingSiteAtoms(receptor, bindingSiteAtoms, grid, false);
			engine = new IdoScore(receptor,bindingSiteAtoms, getReceptorAtomTypes(receptor),grid);
		}
		
		this.mcSteps = mcSteps;
		this.random = new Random(LigandPose.SEED);


	}
	
	public DockingEngine(StereoMolecule receptor, StereoMolecule nativeLigand, double gridDimension) throws DockingFailedException {
		this(receptor,nativeLigand,DEFAULT_NR_MC_STEPS,DEFAULT_START_POSITIONS, gridDimension,ScoringFunction.CHEMPLP);
	}
	
	public DockingEngine(StereoMolecule receptor, StereoMolecule nativeLigand) throws DockingFailedException {
		this(receptor,nativeLigand,DEFAULT_NR_MC_STEPS,DEFAULT_START_POSITIONS, GRID_DIMENSION,ScoringFunction.CHEMPLP);
	}
	
	public void setThreadMaster(ThreadMaster tm) {
		threadMaster = tm;
		shapeDocking.setThreadMaster(tm);
	}
	

	
	/**
	 * generate initial poses: 
	 * 1) shape docking into the negative receptor image
	 * 2) constrained optimization of initial poses to reduce strain energy
	 * @param mol
	 * @param initialPos
	 * @return
	 * @throws DockingFailedException
	 */
	private double getStartingPositions(StereoMolecule mol, List<Conformer> initialPos) throws DockingFailedException {

		double eMin = Double.MAX_VALUE;
		Map<String, Object> ffOptions = new HashMap<String, Object>();
		ffOptions.put("dielectric constant", 80.0);

		ConformerSet confSet = new ConformerSet();
		List<StereoMolecule> alignedMol = shapeDocking.dock(mol);
		alignedMol.stream().forEach(e -> confSet.add(new Conformer(e)));
		for(Conformer c : confSet) {
			if(c!=null) {
				StereoMolecule conf = c.toMolecule();
				conf.ensureHelperArrays(Molecule.cHelperParities);
				
				ForceFieldMMFF94 mmff;
				try {
					mmff = new ForceFieldMMFF94(conf, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
				}
				catch(Exception e) {
					throw new DockingFailedException("could not assess atom types");
				}
				MMFFPositionConstraint constraint = new MMFFPositionConstraint(conf,50,0.5);
				mmff.addEnergyTerm(constraint);
				mmff.minimise();
				Conformer ligConf = new Conformer(conf);
				initialPos.add(ligConf);
				if(initialPos.size()>=startPositions)
					break;
				double e = mmff.getTotalEnergy();
				if(e<eMin)
					eMin = e;
				if(threadMaster!=null && threadMaster.threadMustDie())
					break;
			}
		}
		
		
		return eMin;
		
	}
	/**
	 * 1) find MCS between reference molecule and candidate molecule
	 * 2) 
	 * @param mol
	 */
	private double getStartingPositionsMCS(StereoMolecule mol, List<Conformer> initialPos) {
		double eMin = Double.MAX_VALUE;
		int rotBondsMol = mol.getRotatableBondCount();
		int rotBondsRef = nativeLigand.getRotatableBondCount();
		MCS mcs = new MCS();
		//in MCS, the first molecule should have more rotatable bonds than the second 
		boolean[] bondMCSMol = null;
		boolean[] bondMCSFrag = null;
		boolean isNativeLigFrag;//if native ligand is smaller than the candidate molecule
		if(rotBondsMol > rotBondsRef) {
			mcs.set(mol,nativeLigand);
			bondMCSMol = new boolean[mol.getAllBonds()];
			bondMCSFrag = new boolean[nativeLigand.getAllBonds()];
			isNativeLigFrag = true;
		}
		else {
			mcs.set(mol,nativeLigand);
			bondMCSMol = new boolean[mol.getAllBonds()];
			bondMCSFrag = new boolean[nativeLigand.getAllBonds()];
			isNativeLigFrag = false;
		}

		StereoMolecule mcsMol = mcs.getMCS();
		SSSearcher sss = new SSSearcher();
		sss.setFragment(mcsMol);
		sss.setMolecule(nativeLigand);
		sss.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cMatchDBondToDelocalized);
		int[] map1 = sss.getMatchList().get(0);
		sss.setMolecule(mol);
		sss.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cMatchDBondToDelocalized);
		int[] map2 = sss.getMatchList().get(0);
		mcsConstrainedAtoms = new ArrayList<>();
		for(int a: map2) {
			mcsConstrainedAtoms.add(a);
		}
		//map from reference to mol
		Map<Integer,Integer> map = new HashMap<Integer,Integer>();
		for(int i=0;i<map1.length;i++) {
			map.put(map1[i], map2[i]);
		}
		//find bonds not part of MCS
		mcs.getMCSBondArray(bondMCSMol, bondMCSFrag);
		mcsConstrainedBonds = new ArrayList<>();
		boolean[] bondArray = null;
		if(isNativeLigFrag) 
			bondArray = bondMCSMol;
		else 
			bondArray = bondMCSFrag;
		for(int b=0;b<bondArray.length;b++) {
				if(bondArray[b])
					mcsConstrainedBonds.add(b);
		}
		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		Map<String, Object> ffOptions = new HashMap<String, Object>();
		ffOptions.put("dielectric constant", 80.0);
		ConformerSetGenerator csGen = new ConformerSetGenerator();
		ConformerSet confSet = csGen.generateConformerSet(mol);
		Map<Conformer,Double> confsWithRMSDs = new HashMap<Conformer,Double>();
		for(Conformer conf : confSet) {
			//align MCS using Kabsch algorithm
			Coordinates[] coords1 = new Coordinates[map1.length];
			Coordinates[] coords2 = new Coordinates[map1.length];
			int counter = 0;
			for(int key : map.keySet()) {
				coords1[counter] = new Coordinates(nativeLigand.getCoordinates(key));
				coords2[counter] = new Coordinates(conf.getCoordinates(map.get(key)));
				counter++;
			}
			int[][] mapping = new int[coords1.length][2];
			counter = 0;
			for(int[] m : mapping) {
				m[0] = counter;
				m[1] = counter;
				counter++;
			}
			Coordinates trans1 = new Coordinates();
			Matrix rot = new Matrix(3,3); 
			Coordinates trans2 = new Coordinates();
			KabschAlignment alignment = new KabschAlignment(coords1,coords2,mapping);
			alignment.align(trans1,rot,trans2);
			double rmsd = getCoreRMSD(coords1,coords2);
			for(Coordinates coord : conf.getCoordinates()) {
				coord.add(trans1);
				coord.rotate(rot.getArray());
				coord.add(trans2);
			}
			confsWithRMSDs.put(conf, rmsd);
		}
		
		LinkedHashMap<Conformer,Double> sortedConfs = confsWithRMSDs.entrySet().stream().sorted(Map.Entry.<Conformer,Double>comparingByValue()).collect(Collectors.toMap(Map.Entry::getKey,Map.Entry::getValue,
				(e1, e2) -> e1, LinkedHashMap::new));
		Iterator<Entry<Conformer,Double>> iterator = sortedConfs.entrySet().iterator();
		boolean done = false;
		Entry<Conformer,Double> entry;
		int c=0;
		while(!done && c<startPositions*MCS_EXHAUSTIVENESS) {
			c++;
			try	{
				entry = iterator.next();
			}
			catch(Exception e) {
				done = true;
				continue;
			}
			Conformer conf = entry.getKey();
			StereoMolecule aligned = new StereoMolecule();
			aligned = conf.toMolecule(null);
			aligned.ensureHelperArrays(Molecule.cHelperParities);
			ForceFieldMMFF94 mmff;
			MMFFPositionConstraint constraint = new MMFFPositionConstraint(aligned,50,0.2);
			mmff = new ForceFieldMMFF94(aligned, ForceFieldMMFF94.MMFF94SPLUS,ffOptions);
			mmff.addEnergyTerm(constraint);
			mmff.minimise();
			double e = mmff.getTotalEnergy();
			for(int a=0;a<aligned.getAllAtoms();a++) {
				conf.setCoordinates(a, new Coordinates(aligned.getCoordinates(a)));
			}
			initialPos.add(conf);
			if(e<eMin)
				eMin = e;
		}
		return eMin;
	
	}
	
	
	public DockingResult dockMolecule(StereoMolecule mol) throws DockingFailedException {
		Conformer bestPose = null;
		double bestEnergy = Double.MAX_VALUE;
		if(ForceFieldMMFF94.table(ForceFieldMMFF94.MMFF94SPLUS)==null)
			ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		List<Conformer> startPoints = new ArrayList<>();
		double eMin = 0.0;
		int steps = mcSteps;
		if(mcsRef!=null) {
			eMin = getStartingPositionsMCS(mol, startPoints);
			steps=steps/MCS_EXHAUSTIVENESS;
		}
		else {
			eMin = getStartingPositions(mol, startPoints);
		}

		Map<String,Double> contributions = null;
		for(Conformer ligConf : startPoints) {
			Conformer newLigConf = new Conformer(ligConf);		
			LigandPose pose = initiate(newLigConf,eMin);
			if(mcsRef!=null) {
				pose.setMCSBondConstraints(mcsConstrainedBonds);
				for(int a : mcsConstrainedAtoms) {
					PositionConstraint constr = new PositionConstraint(newLigConf,a,50,1.0);
					pose.addConstraint(constr);
				}
				
			}
			double energy = mcSearch(pose,steps);
			if(energy<bestEnergy) {
				bestEnergy = pose.getScore();
				bestPose = pose.getLigConf();
				contributions = pose.getContributions();
			}
			if(threadMaster!=null && threadMaster.threadMustDie())
				break;
		}
		
		if(bestPose!=null) {
			StereoMolecule best = bestPose.toMolecule();
			Rotation rot = rotation.getInvert();
			Translation translate = new Translation(new double[] {origCOM.x, origCOM.y, origCOM.z});
			rot.apply(best);
			translate.apply(best);
			return new DockingResult(mol, best,bestEnergy,contributions);
		}
		else {
			throw new DockingFailedException("docking failed");
		}
		
	}
	
	

	/**
	 * use monte carlo steps to permute molecular rotation, translation, torsion angles
	 * promising poses (below a certain cutoff) are optimized
	 * @param pose
	 * @return
	 */
	private double mcSearch(LigandPose pose, int steps) {
		double[] bestState = new double[pose.getState().length];
		double[] oldState = new double[pose.getState().length];
		double[] state = new double[pose.getState().length];
		double bestEnergy = -Float.MAX_VALUE;
		double oldEnergy = -Float.MAX_VALUE;
		double energy = -Float.MAX_VALUE;
		OptimizerLBFGS optimizer = new OptimizerLBFGS(200,0.001);
		oldState = optimizer.optimize(pose);
		bestState = oldState;
		oldEnergy = pose.getFGValue(new double[bestState.length]);
		bestEnergy = oldEnergy;
	
		for(int i=0;i<steps;i++) {
			pose.randomPerturbation();
			double energyMC = pose.getFGValue(new double[bestState.length]);
			if(energyMC<MINI_CUTOFF) {
				state = optimizer.optimize(pose);	
				energy = pose.getFGValue(new double[bestState.length]);
			}
			else {
				state = pose.getState();
				energy=energyMC;
			}


			if(energy<oldEnergy) { //accept new pose
				oldEnergy = energy;
				oldState = state;
				if(energy<bestEnergy) {
					bestEnergy = energy;
					bestState = state;
				}
			}
			else {
				double dE = energy-oldEnergy;
				double randNr = random.nextDouble();
				double probability = Math.exp(-dE/BOLTZMANN_FACTOR); // Metropolis-Hastings
				if(randNr<probability) { //accept
					oldEnergy = energy;
					oldState = state;
				}
				else { //reject
					pose.setState(oldState);
				}
				
			}
			
		}

		pose.setState(bestState);
		pose.removeConstraints();
		bestEnergy = pose.getScore();
		return bestEnergy;
		
	}
	
	private LigandPose initiate(Conformer ligConf, double e0) {
		LigandPose pose = new LigandPose(ligConf, engine, e0);
		
		return pose;
		
		
	}
	
	public void setMCSReference(StereoMolecule referencePose) {
		mcsRef = referencePose;
	}
	
	public static void getBindingSiteAtoms(StereoMolecule receptor, Set<Integer> bindingSiteAtoms, MoleculeGrid grid,
			boolean includeHydrogens) {
		int[] gridSize = grid.getGridSize();
		int atoms = receptor.getAtoms();
		if(includeHydrogens)
			atoms = receptor.getAllAtoms();
		
		for(int i=0;i<atoms;i++) {
			int[] gridC = grid.getGridCoordinates(receptor.getCoordinates(i));
			int x = gridC[0];
			int y = gridC[1];
			int z = gridC[2];	
			if(x>0 && x<gridSize[0]) {
				if(y>0 && y<gridSize[1]) {
					if(z>0 && z<gridSize[2]) {
						bindingSiteAtoms.add(i);
					}
			}
		}
		}
		
	}
	
	public static int[] getReceptorAtomTypes(StereoMolecule receptor) {
		int[] receptorAtomTypes = new int[receptor.getAtoms()];
		for(int i=0;i<receptor.getAtoms();i++) {
			receptorAtomTypes[i] = InteractionAtomTypeCalculator.getAtomType(receptor, i);
		}
		return receptorAtomTypes;
		
		
	}
	/**
	 * rotates receptor and ligand to principal moments of inertia of ligand, for efficient grid creation
	 * @param receptor
	 * @param ligand
	 * @return
	 */
	private void preprocess(StereoMolecule receptor, StereoMolecule ligand) {
		Translation translate = new Translation(new double[] {-origCOM.x,-origCOM.y,-origCOM.z});
		translate.apply(ligand);
		rotation.apply(ligand);
		translate.apply(receptor);
		rotation.apply(receptor);
		
	}
	/**
	 * the parameter d defines how much the atoms are allowed to move from their original position
	 * @param d
	 * @return
	 */
	
	public double refineNativePose(double d, double[] coords) {
		Map<String, Object> ffOptions = new HashMap<String, Object>();
		ffOptions.put("dielectric constant", 80.0);
		if(ForceFieldMMFF94.table(ForceFieldMMFF94.MMFF94SPLUS)==null)
			ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		Molecule3D nativePose = new Molecule3D(nativeLigand);
		new Canonizer(nativePose);
		ForceFieldMMFF94 mmff = new ForceFieldMMFF94(nativePose, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
		MMFFPositionConstraint constraint = new MMFFPositionConstraint(nativePose,50,0.2);
		mmff.addEnergyTerm(constraint);
		mmff.minimise();
		ConformerSetGenerator confSetGen = new ConformerSetGenerator(100,ConformerGenerator.STRATEGY_LIKELY_RANDOM, false,
				LigandPose.SEED);
		ConformerSet confSet = confSetGen.generateConformerSet(nativePose);
		double eMin = Double.MAX_VALUE;

		ConformerSet initialPos = new ConformerSet();
		for(Conformer conformer : confSet) {
			if(conformer!=null) {
				StereoMolecule conf = conformer.toMolecule();
				conf.ensureHelperArrays(Molecule.cHelperParities);
				mmff = new ForceFieldMMFF94(conf, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
				constraint = new MMFFPositionConstraint(conf,50,0.5);
				mmff.addEnergyTerm(constraint);
				mmff.minimise();
				Conformer ligConf = new Conformer(conf);
				initialPos.add(ligConf);
				if(initialPos.size()>=startPositions)
					break;
				double e = mmff.getTotalEnergy();
				if(e<eMin)
					eMin = e;
			}
		}
		LigandPose pose = new LigandPose(new Conformer(nativePose), engine, eMin);
		pose.addPositionalConstraints(d);
		OptimizerLBFGS optimizer = new OptimizerLBFGS(200,0.001);
		optimizer.optimize(pose);
		double energy = pose.getFGValue(new double[pose.getState().length]);
		StereoMolecule best = pose.getLigConf().toMolecule();
		double[][] rot = rotation.getInvert().getRotation();
		PheSAAlignment.rotateMol(best, rot);
		PheSAAlignment.translateMol(best, new double[] {origCOM.x, origCOM.y, origCOM.z} );
		for(int a=0;a<best.getAllAtoms();a++) {
			coords[3*a] = best.getAtomX(a);
			coords[3*a+1] = best.getAtomY(a);
			coords[3*a+2] = best.getAtomZ(a);
		}
	
		
		return energy;
		
		
	}
	
	private static double getCoreRMSD(Coordinates[] coords1, Coordinates[] coords2) {
		double rmsd = 0.0;
		for(int i=0;i<coords1.length;i++) { 
			Coordinates c1 = coords1[i];
			Coordinates c2 = coords2[i];
			rmsd+=c1.distanceSquared(c2);
		}
		rmsd/=coords1.length;
		rmsd = Math.sqrt(rmsd);
		return rmsd;
	}
	
	public static class DockingResult implements Comparable<DockingResult>  {
		private double score;
		private StereoMolecule input; //might be a different enantiomer/protomer than the pose
		private StereoMolecule pose;
		private Map<String,Double> contributions;
		private static final String DELIMITER = ";";
		private static final String DELIMITER2 = ":";
		private static final String DELIMITER3 = "%";
		private static final String NULL_CONTRIBUTION = "#";
		
		public DockingResult(StereoMolecule input,
				StereoMolecule pose, double score, Map<String,Double> contributions) {
			this.score = score;
			this.pose = pose;
			this.contributions = contributions;
			this.input = input;
		}
		
		
		
		public void setInput(StereoMolecule input) {
			this.input = input;
		}

		public StereoMolecule getInput() {
			return input;
		}


		public double getScore() {
			return score;
		}
		
		public StereoMolecule getPose() {
			return pose;
		}
		
		public Map<String,Double> getContributions() {
			return contributions;
		}
		
		
		public String encode() {
			Encoder encoder = Base64.getEncoder();
			StringBuilder sb = new StringBuilder();
			Canonizer can = new Canonizer(pose, Canonizer.COORDS_ARE_3D);
			String idcoords = can.getEncodedCoordinates(true);
			String idcode = can.getIDCode();
			sb.append(idcode);
			sb.append(DELIMITER);
			sb.append(idcoords);
			sb.append(DELIMITER);
			sb.append(input.getIDCode());
			sb.append(DELIMITER);
			sb.append(encoder.encodeToString(EncodeFunctions.doubleToByteArray(score)));
			sb.append(DELIMITER);
			if(contributions==null || contributions.keySet().size()==0)
				sb.append(NULL_CONTRIBUTION);
			else {
				for(String name : contributions.keySet()) {
					sb.append(name);
					sb.append(DELIMITER3);
					sb.append(encoder.encodeToString(EncodeFunctions.doubleToByteArray(contributions.get(name))));
					sb.append(DELIMITER2);
				}
				sb.setLength(sb.length() - 1);
			}
			
			return sb.toString();
		}
		
		public static DockingResult decode(String resultString) {
			Decoder decoder = Base64.getDecoder();
			String[] s = resultString.split(DELIMITER);
			String idcode = s[0];
			String idcoords = s[1];
			StereoMolecule pose = new StereoMolecule();
			IDCodeParserWithoutCoordinateInvention parser = new IDCodeParserWithoutCoordinateInvention();
			parser.parse(pose, idcode, idcoords);
			pose.ensureHelperArrays(Molecule.cHelperCIP);
			String idcodeInput = s[2];
			StereoMolecule input = new StereoMolecule();
			new IDCodeParser().parse(input, idcodeInput);
			double score = EncodeFunctions.byteArrayToDouble(decoder.decode(s[3].getBytes(StandardCharsets.UTF_8)));
			Map<String,Double> contributions = null;
			if(!s[4].equals(NULL_CONTRIBUTION)) {
				contributions = new HashMap<String,Double>();
				String[] splitted = s[4].split(DELIMITER2);
				for(String contr : splitted) {
					String[] splitted2 = contr.split(DELIMITER3);
					String name = splitted2[0];
					double value = EncodeFunctions.byteArrayToDouble(decoder.decode(splitted2[1].getBytes(StandardCharsets.UTF_8)));
					contributions.put(name, value);
				}
			}
				
			DockingResult dockingResult = new DockingResult(input,pose,score,contributions);
			return dockingResult;
		}

		@Override
		public int compareTo(DockingResult o) {
            if(Double.isNaN(this.score)&& Double.isNaN(o.score)) { return 0; }
            if(Double.isNaN(this.score)) { return -1; }
            if(Double.isNaN(o.score)) { return  1; }

            return Double.compare( this.score, o.score);
           
		}
		
	}
	


	
	
	
	

}
