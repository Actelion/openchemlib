package com.actelion.research.chem.docking;


import java.util.ArrayList;
import java.util.Base64;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.Base64.Decoder;
import java.util.Base64.Encoder;

import org.openmolecules.chem.conf.gen.ConformerGenerator;

import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.IDCodeParserWithoutCoordinateInvention;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.alignment3d.transformation.ExponentialMap;
import com.actelion.research.chem.alignment3d.transformation.Quaternion;
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
import com.actelion.research.chem.forcefield.mmff.PositionConstraint;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import com.actelion.research.chem.io.pdb.converter.MoleculeGrid;
import com.actelion.research.chem.optimization.OptimizerLBFGS;
import com.actelion.research.chem.phesa.EncodeFunctions;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.PheSAAlignment;
import com.actelion.research.chem.phesa.ShapeVolume;
import com.actelion.research.chem.phesa.PheSAAlignment.PheSAResult;

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
	

	private Rotation rotation; //for initial prealignment to native ligand
	private Coordinates origCOM;
	private int mcSteps;
	private Random random;
	private AbstractScoringEngine engine;
	private int startPositions;
	private StereoMolecule nativeLigand;
	private ShapeDocking shapeDocking;
	private ThreadMaster threadMaster;

	
	public DockingEngine(StereoMolecule rec, StereoMolecule nativeLig, int mcSteps, int startPositions,
			ScoringFunction scoringFunction) {

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
		
		MoleculeGrid grid = new MoleculeGrid(nativeLigand,DockingEngine.GRID_RESOLUTION,
				new Coordinates(DockingEngine.GRID_DIMENSION,DockingEngine.GRID_DIMENSION,
						DockingEngine.GRID_DIMENSION));
		
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
	
	public DockingEngine(StereoMolecule receptor, StereoMolecule nativeLigand) {
		this(receptor,nativeLigand,DEFAULT_NR_MC_STEPS,DEFAULT_START_POSITIONS,ScoringFunction.CHEMPLP);
	}
	
	public void setThreadMaster(ThreadMaster tm) {
		threadMaster = tm;
	}
	
	private double getStartingPositions(StereoMolecule mol, List<Conformer> initialPos) throws DockingFailedException {

		double eMin = Double.MAX_VALUE;
		Map<String, Object> ffOptions = new HashMap<String, Object>();
		ffOptions.put("dielectric constant", 80.0);
		//ffOptions.put("angle bend", false);
		//ffOptions.put("stretch bend", false);
		//ffOptions.put("bond stretch", false);
		//ffOptions.put("out of plane", false);

		ConformerSet confSet = new ConformerSet();
		long t0 = System.currentTimeMillis();
		List<StereoMolecule> alignedMol = shapeDocking.dock(mol);
		long t1 = System.currentTimeMillis();
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
				PositionConstraint constraint = new PositionConstraint(conf,50,0.2);
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
	
	
	
	public DockingResult dockMolecule(StereoMolecule mol) throws DockingFailedException {
		Conformer bestPose = null;
		double bestEnergy = Double.MAX_VALUE;
		if(ForceFieldMMFF94.table(ForceFieldMMFF94.MMFF94SPLUS)==null)
			ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		List<Conformer> startPoints = new ArrayList<>();
		double eMin = getStartingPositions(mol, startPoints);
		for(Conformer ligConf : startPoints) {
			for(double[] transform : PheSAAlignment.initialTransform(0)) {
				Conformer newLigConf = new Conformer(ligConf);
				ExponentialMap eMap = new ExponentialMap(transform[0],transform[1],transform[2]);
				Quaternion q = eMap.toQuaternion();
				Coordinates com = DockingUtils.getCOM(newLigConf);
				Rotation rot = new Rotation(q.getRotMatrix().getArray());
				Translation trans = new Translation(new double[] {transform[3],transform[4],transform[5]});
				TransformationSequence t = new TransformationSequence();
				Translation t1 = new Translation(-com.x,-com.y,-com.z);
				Translation t2 = new Translation(com.x,com.y,com.z);
				t.addTransformation(t1);
				t.addTransformation(rot);
				t.addTransformation(t2);
				t.addTransformation(trans);
				t.apply(newLigConf);
				LigandPose pose = initiate(newLigConf,eMin);
				double energy = mcSearch(pose);
				if(energy<bestEnergy) {
					bestEnergy = energy;
					bestPose = pose.getLigConf();
				}
				if(threadMaster!=null && threadMaster.threadMustDie())
					break;
			}
		}
		if(bestPose!=null) {
			StereoMolecule best = bestPose.toMolecule();
			Rotation rot = rotation.getInvert();
			Translation translate = new Translation(new double[] {origCOM.x, origCOM.y, origCOM.z});
			rot.apply(best);
			translate.apply(best);

			return new DockingResult(best,bestEnergy);
		}
		else {
			throw new DockingFailedException("docking failed");
		}
		
	}
	
	

	
	private double mcSearch(LigandPose pose) {
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
	
		for(int i=0;i<mcSteps;i++) {
			pose.randomPerturbation(random);
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
		engine.getScore();
		return bestEnergy;
		
	}
	
	private LigandPose initiate(Conformer ligConf, double e0) {
		LigandPose pose = new LigandPose(ligConf, engine, e0);
		
		return pose;
		
		
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
		//ffOptions.put("angle bend", false);
		//ffOptions.put("stretch bend", false);
		//ffOptions.put("bond stretch", false);
		if(ForceFieldMMFF94.table(ForceFieldMMFF94.MMFF94SPLUS)==null)
			ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		Molecule3D nativePose = new Molecule3D(nativeLigand);
		new Canonizer(nativePose);
		ForceFieldMMFF94 mmff = new ForceFieldMMFF94(nativePose, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
		PositionConstraint constraint = new PositionConstraint(nativePose,50,0.2);
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
				constraint = new PositionConstraint(conf,50,0.2);
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
	
	public static class DockingResult implements Comparable<DockingResult>  {
		private double score;
		private StereoMolecule pose;
		private static final String DELIMITER = ";";
		
		public DockingResult(StereoMolecule pose, double score) {
			this.score = score;
			this.pose = pose;
		}
		
		public double getScore() {
			return score;
		}
		
		public StereoMolecule getPose() {
			return pose;
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
			sb.append(encoder.encodeToString(EncodeFunctions.doubleToByteArray(score)));
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
			double score = EncodeFunctions.byteArrayToDouble(decoder.decode(s[2].getBytes()));
			DockingResult dockingResult = new DockingResult(pose,score);
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
