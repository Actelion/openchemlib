/*
 * Copyright 2013-2020 Thomas Sander, openmolecules.org
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.so;

import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.*;
import com.actelion.research.util.DoubleFormat;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

public class ConformationSelfOrganizer {
	private static final int    INITIAL_POOL_SIZE = 4;
	private static final int    MAX_CONFORMER_TRIES = 6;
	private static final int    MAX_BREAKOUT_ROUNDS = 3;
	private static final int    PREPARATION_CYCLES = 40;		// without ambiguous rules (i.e. torsion)
	private static final int    PRE_OPTIMIZATION_CYCLES = 20;	// with ambiguous rules before breakout
	private static final int    BREAKOUT_CYCLES = 20;
	private static final int    OPTIMIZATION_CYCLES = 100;
	private static final int    MINIMIZATION_CYCLES = 20;
	private static final double	STANDARD_CYCLE_FACTOR = 1.0;
	private static final double	MINIMIZATION_REDUCTION = 20.0;
	private static final double ATOM_FLAT_RING_BREAKOUT_STRAIN = 0.25;
	private static final double ATOM_CAGE_BREAKOUT_STRAIN = 2.0;
	private static final double	BREAKOUT_DISTANCE = 8.0;
	private static final double MAX_AVERAGE_ATOM_STRAIN = 0.025;
	private static final double MAX_HIGHEST_ATOM_STRAIN = 0.05;
	private static final double MAX_STRAIN_TOLERANCE = 1.5;

public static boolean KEEP_INITIAL_COORDINATES = false;
public static boolean WRITE_DW_FILE = false;
private static final String DATAWARRIOR_DEBUG_FILE = "/home/thomas/data/debug/conformationSampler.dwar";
private BufferedWriter mDWWriter;
private Conformer mLastDWConformer;
private int mDWCycle;
private double[] mDWStrain; 	// TODO get rid of this

	private StereoMolecule		mMol;
    private Random				mRandom;
    private int					mMaxConformers;
    private boolean				mPoolIsClosed;
	private ArrayList<ConformationRule> mRuleList;
	private ArrayList<SelfOrganizedConformer> mConformerList;
	private double              mMinAverageAtomStrainInPool,mMinHighestAtomStrainInPool;
	private int[]				mRuleCount;
	private boolean[]			mSkipRule;
	private int[]				mRotatableBondForDescriptor;
	private ThreadMaster        mThreadMaster;

	/**
	 * Generates a new ConformationSelfOrganizer from the given molecule.
	 * Explicit hydrogens are removed from the molecule, unless the
	 * keepHydrogen flag is set.<br>
	 * One conformer can be generated with the getOneConformer()
	 * or getOneConformerInPlace(). <br>
	 * Multiple different conformers can be generated with initializeConformers()
	 * and getNextConformer(). In this case conformers are considered different
	 * if at least one dihedral angle of a rotatable bond is substantially different.
	 * If atoms of the molecule are marked, these are not considered part of the molecule,
	 * when the rotatable bonds for the difference check are located.
	 * @param mol
	 */
	public ConformationSelfOrganizer(final StereoMolecule mol, boolean keepHydrogen) {

/*// winkel zwischen zwei vektoren:
final double[] v1 = { Math.sqrt(3.0)/2.0, 0.5, 0 };
final double[] v2 = { -1, 0, 1 };
double cosa = (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
			/ (Math.sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])
			 * Math.sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]));
double a = Math.acos(cosa);
System.out.println("angle:"+a+"  in degrees:"+(a*180/Math.PI));
*/
		mMol = mol;
		if (!keepHydrogen)
		    mMol.removeExplicitHydrogens();
		mMol.ensureHelperArrays(Molecule.cHelperParities);

		mRuleList = new ArrayList<ConformationRule>();
		mSkipRule = new boolean[ConformationRule.RULE_NAME.length];
		mRuleCount = new int[ConformationRule.RULE_NAME.length];

		DistanceRule.calculateRules(mRuleList, mol);
		mRuleCount[ConformationRule.RULE_TYPE_DISTANCE] = mRuleList.size();

		mRuleCount[ConformationRule.RULE_TYPE_PLANE] = -mRuleList.size();
		PlaneRule.calculateRules(mRuleList, mol);
		mRuleCount[ConformationRule.RULE_TYPE_PLANE] += mRuleList.size();

		mRuleCount[ConformationRule.RULE_TYPE_LINE] = -mRuleList.size();
		StraightLineRule.calculateRules(mRuleList, mol);
		mRuleCount[ConformationRule.RULE_TYPE_LINE] += mRuleList.size();

		mRuleCount[ConformationRule.RULE_TYPE_STEREO] = -mRuleList.size();
		TetrahedralStereoRule.calculateRules(mRuleList, mol);
		mRuleCount[ConformationRule.RULE_TYPE_STEREO] += mRuleList.size();

		mRuleCount[ConformationRule.RULE_TYPE_BINAP] = -mRuleList.size();
		AxialStereoRule.calculateRules(mRuleList, mol);
		mRuleCount[ConformationRule.RULE_TYPE_BINAP] += mRuleList.size();

		mRuleCount[ConformationRule.RULE_TYPE_TORSION] = -mRuleList.size();
		TorsionRule.calculateRules(mRuleList, mol);
		mRuleCount[ConformationRule.RULE_TYPE_TORSION] += mRuleList.size();

//		listRules();
		}

/*	private void listRules() {
		System.out.println("---------------------------------------------------------------------");
		for (int i=0; i<mMol.getAllAtoms(); i++) {
			System.out.print(""+i+" "+Molecule.cAtomLabel[mMol.getAtomicNo(i)]);
			for (int j=0; j<mMol.getAllConnAtoms(i); j++) {
				int connBond = mMol.getConnBond(i, j);
				if (mMol.isAromaticBond(connBond))
					System.out.print(" .");
				else if (mMol.getBondOrder(connBond) == 0)
					System.out.print(" ~");
				else if (mMol.getBondOrder(connBond) == 1)
					System.out.print(" -");
				else if (mMol.getBondOrder(connBond) == 2)
					System.out.print(" =");
				else if (mMol.getBondOrder(connBond) == 3)
					System.out.print(" #");
				System.out.print(""+mMol.getConnAtom(i, j));
				}
			if (mMol.getAtomParity(i) != 0)
				System.out.print(" parity:"+mMol.getAtomParity(i));
			System.out.println();
			}
		System.out.print("CanRanks:");
		mMol.ensureHelperArrays(Molecule.cHelperBitSymmetrySimple);
		for (int i=0; i<mMol.getAtoms(); i++) {
			System.out.print(" "+i+":"+mMol.getSymmetryRank(i));
			if (mMol.getAtomParity(i) != 0)
				System.out.print(" absTHP:"+mMol.getAbsoluteAtomParity(i));
			}
		System.out.println();
    	BondLengthSet bondLengthSet = new BondLengthSet(mMol);
		BondAngleSet bondAngleSet = new BondAngleSet(mMol, bondLengthSet);
		System.out.println("Angles:");
		for (int i=0; i<mMol.getAtoms(); i++) {
		for (int j=1; j<mMol.getAllConnAtoms(i); j++) {
		for (int k=0; k<j; k++) {
		System.out.print(bondAngleSet.getConnAngle(i, j, k)+" ");
		}}
		System.out.println();
		}

		for (ConformationRule rule:mRuleList)
			System.out.println(rule.toString());
		}
*/

	public ArrayList<ConformationRule> getRuleList() {
		return mRuleList;
	}

	/**
	 * @return returns the molecule that was passed to the constructor.
	 */
	public StereoMolecule getMolecule() {
		return mMol;
	}

	public void setThreadMaster(ThreadMaster tm) {
		mThreadMaster = tm;
		}

	/**
	 * This convenience method returns the StereoMolecule that has been passed
	 * to the constructor after modifying its atom coordinates
	 * to reflect the conformer internally created by generateOneConformer().
	 * @return 
	 */
	public StereoMolecule generateOneConformerInPlace(long randomSeed) {
		SelfOrganizedConformer conformer = generateOneConformer(randomSeed);
		return (conformer == null) ? null : conformer.toMolecule(mMol);
		}

	/**
	 * Generates the coordinates for one conformer in the calling thread.
	 * This is done by trying MAX_CONFORMER_TRIES times to create a random
	 * conformer that meets MAX_ATOM_STRAIN and MAX_TOTAL_STRAIN criteria.
	 * If one is found it is returned. Otherwise the conformer with the lowest
	 * total strain is returned.<br>
	 * <b>Note:</b> If randomSeed is different from 0, then only one conformer
	 * is done produced and returned, no matter whether it meets any strain criteria.
     * @param randomSeed 0 or specific seed
	 */
	public SelfOrganizedConformer generateOneConformer(long randomSeed) {
        mRandom = (randomSeed == 0) ? new Random() : new Random(randomSeed);

        SelfOrganizedConformer conformer = new SelfOrganizedConformer(mMol);

		if (WRITE_DW_FILE) {
			try {
				writeDWFileStart();
				mDWCycle = 0;
				mLastDWConformer = null;
				tryGenerateConformer(conformer);
				writeDWFileEnd();
				mDWWriter.close();
				return conformer;
				}
			catch (IOException e) {
				e.printStackTrace();
				return null;
				}
			}

		SelfOrganizedConformer bestConformer = null;
		for (int i=0; i<MAX_CONFORMER_TRIES; i++) {
			if (mThreadMaster != null && mThreadMaster.threadMustDie())
				break;

			if (tryGenerateConformer(conformer) || randomSeed != 0L)
				return conformer;	// sufficiently low strain, we take this and don't try generating better ones

			if (bestConformer == null) {
				bestConformer = conformer;
				conformer = new SelfOrganizedConformer(mMol);
				}
			else if (bestConformer.isWorseThan(conformer)) {
				SelfOrganizedConformer tempConformer = bestConformer;
				bestConformer = conformer;
				conformer = tempConformer;
				}
			}

		return bestConformer;
		}

	/**
	 * Needs to be called, before getting individual conformers of the same molecule by
	 * getNextConformer(). Depending on the flexibility of the molecule, this method creates
	 * a small pool of random conformers, of which getNextConformer() always picks that
	 * conformer that is an optimum concerning atom strain and diversity to the already picked ones.
	 * If the pool is getting low on conformers, then new ones are generated on-the-fly.
	 * @param randomSeed use a value != 0 for a reproducible random number sequence
	 * @param maxConformers -1 to automatically generate maximum on degrees of freedom
	 */
	public void initializeConformers(long randomSeed, int maxConformers) {
        mRandom = (randomSeed == 0) ? new Random() : new Random(randomSeed);

        mConformerList = new ArrayList<SelfOrganizedConformer>();
		mMinHighestAtomStrainInPool = MAX_STRAIN_TOLERANCE * MAX_HIGHEST_ATOM_STRAIN;
		mMinAverageAtomStrainInPool = MAX_STRAIN_TOLERANCE * MAX_AVERAGE_ATOM_STRAIN;
        mPoolIsClosed = false;

		int freeBondCount = 0;
		int ringBondCount = 0;
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (!mMol.isAromaticBond(bond)
			 && mMol.getBondOrder(bond) == 1
			 && mMol.getAllConnAtoms(mMol.getBondAtom(0, bond)) > 1
			 && mMol.getAllConnAtoms(mMol.getBondAtom(1, bond)) > 1) {
				if (!mMol.isRingBond(bond))
					freeBondCount++;
				else if (mMol.getBondRingSize(bond) > 4)
					ringBondCount++;
				}
			}

		mMaxConformers = (maxConformers == -1) ? (1 << (1+freeBondCount+ringBondCount/2)) : maxConformers;
		increaseConformerPool(Math.min(INITIAL_POOL_SIZE, mMaxConformers));
		}

	private void increaseConformerPool(int newConformerCount) {
		SelfOrganizedConformer bestRefusedConformer = null;
		SelfOrganizedConformer conformer = null;
		int finalPoolSize = mConformerList.size() + newConformerCount;
		int tryCount = newConformerCount * MAX_CONFORMER_TRIES;
		for (int i=0; i<tryCount && mConformerList.size()<finalPoolSize; i++) {
			if (mThreadMaster != null && mThreadMaster.threadMustDie())
				break;

			if (conformer == null)
				conformer = new SelfOrganizedConformer(mMol);
			if (tryGenerateConformer(conformer)) {
				if (addConformerIfNew(conformer))
					conformer = null;
				}
			else if (conformer.getTotalStrain() / conformer.getSize() < MAX_STRAIN_TOLERANCE * mMinAverageAtomStrainInPool
				  && conformer.getHighestAtomStrain() < MAX_STRAIN_TOLERANCE * mMinHighestAtomStrainInPool) {
				if (addConformerIfNew(conformer))
					conformer = null;
				}
			else {
				if (bestRefusedConformer == null) {
					bestRefusedConformer = conformer;
					conformer = null;
					}
				else if (bestRefusedConformer.isWorseThan(conformer)) {
					SelfOrganizedConformer tempConformer = bestRefusedConformer;
					bestRefusedConformer = conformer;
					conformer = tempConformer;
					}
				}
			}
		if (mConformerList.isEmpty() && bestRefusedConformer != null)
			mConformerList.add(bestRefusedConformer);
		if (mConformerList.size() < finalPoolSize
		 || mConformerList.size() == mMaxConformers)
			mPoolIsClosed = true;
		}

	private boolean addConformerIfNew(SelfOrganizedConformer conformer) {
		if (mRotatableBondForDescriptor == null)
			mRotatableBondForDescriptor = TorsionDescriptorHelper.findRotatableBonds(getMolecule());

		conformer.calculateDescriptor(mRotatableBondForDescriptor);
		boolean isNew = true;
		for (SelfOrganizedConformer c:mConformerList) {
			if (conformer.equals(c)) {
				isNew = false;
				break;
			}
		}
		if (!isNew)
			return false;

		mConformerList.add(conformer);

		double averageStrain = conformer.getTotalStrain() / conformer.getSize();
		double highestStrain = conformer.getHighestAtomStrain();
		if ((mMinAverageAtomStrainInPool > averageStrain)
		 || (mMinHighestAtomStrainInPool > highestStrain)) {
			if (mMinAverageAtomStrainInPool > averageStrain)
				mMinAverageAtomStrainInPool = averageStrain;
			if (mMinHighestAtomStrainInPool > highestStrain)
				mMinHighestAtomStrainInPool = highestStrain;
			for (int j=mConformerList.size()-1; j>=0; j--) {
				SelfOrganizedConformer soc = mConformerList.get(j);
				if (!soc.isAcceptable(mRuleList)
						&& (soc.getTotalStrain() / soc.getSize() > MAX_STRAIN_TOLERANCE * mMinAverageAtomStrainInPool
						|| soc.getHighestAtomStrain() > MAX_STRAIN_TOLERANCE * mMinHighestAtomStrainInPool))
					mConformerList.remove(j);
				}
			}

		return true;
		}

	/**
	 * Picks a new conformer from the conformer pool created by initializeConformers().
	 * Low strain conformers and conformers being most different from already selected ones
	 * are selected first. If the pool is getting short on conformers, new conformers are
	 * created as long as molecule flexibility allows. If a representative set of low strain
	 * molecules have been picked, this method returns null, provided that at least one conformer
	 * was returned.
	 * @return
	 */
	public SelfOrganizedConformer getNextConformer() {
		if (mConformerList == null)
			return null;

		if (!mPoolIsClosed)
			increaseConformerPool(1);

		SelfOrganizedConformer bestConformer = null;
		for (SelfOrganizedConformer conformer:mConformerList) {
			if (!conformer.isUsed()
			 && (bestConformer == null || bestConformer.isWorseThan(conformer))) {
				bestConformer = conformer;
				}
			}

		if (bestConformer != null)
			bestConformer.setUsed(true);
		else
			mConformerList = null;

		return bestConformer;
		}

	private void writeDWFileStart() throws IOException {
        mDWWriter = new BufferedWriter(new FileWriter(DATAWARRIOR_DEBUG_FILE));
        mDWWriter.write("<column properties>");
        mDWWriter.newLine();
        mDWWriter.write("<columnName=\"Structure\">");
        mDWWriter.newLine();
        mDWWriter.write("<columnProperty=\"specialType\tidcode\">");
        mDWWriter.newLine();
        mDWWriter.write("<columnName=\"before\">");
        mDWWriter.newLine();
        mDWWriter.write("<columnProperty=\"specialType\tidcoordinates3D\">");
        mDWWriter.newLine();
        mDWWriter.write("<columnProperty=\"parent\tStructure\">");
        mDWWriter.newLine();
        mDWWriter.write("<columnName=\"after\">");
        mDWWriter.newLine();
        mDWWriter.write("<columnProperty=\"specialType\tidcoordinates3D\">");
        mDWWriter.newLine();
        mDWWriter.write("<columnProperty=\"parent\tStructure\">");
        mDWWriter.newLine();
        mDWWriter.write("</column properties>");
        mDWWriter.newLine();
        mDWWriter.write("Structure\tbefore\tafter\tcycle\truleName\truleAtoms\truleDetail");
        for (int i=0; i<ConformationRule.RULE_NAME.length; i++)
            mDWWriter.write("\t"+ConformationRule.RULE_NAME[i]);
        mDWWriter.write("\ttotalStrain\tstrainGain\truleStrainBefore\truleStrainAfter\truleStrainGain");
        mDWWriter.newLine();
        }


	private void writeDWFileEnd() throws IOException {
		mDWWriter.write("<datawarrior properties>");
        mDWWriter.newLine();
		mDWWriter.write("<axisColumn_2D View_0=\"cycle\">");
        mDWWriter.newLine();
		mDWWriter.write("<axisColumn_2D View_1=\"totalStrain\">");
        mDWWriter.newLine();
		mDWWriter.write("<chartType_2D View=\"scatter\">");
        mDWWriter.newLine();
		mDWWriter.write("<colorColumn_2D View=\"ruleName\">");
        mDWWriter.newLine();
		mDWWriter.write("<colorCount_2D View=\"3\">");
        mDWWriter.newLine();
		mDWWriter.write("<colorListMode_2D View=\"Categories\">");
        mDWWriter.newLine();
		mDWWriter.write("<color_2D View_0=\"-11992833\">");
        mDWWriter.newLine();
		mDWWriter.write("<color_2D View_1=\"-65494\">");
        mDWWriter.newLine();
		mDWWriter.write("<color_2D View_2=\"-16732826\">");
        mDWWriter.newLine();
		mDWWriter.write("<detailView=\"height[Data]=0.4;height[before]=0.3;height[after]=0.3\">");
        mDWWriter.newLine();
		mDWWriter.write("<mainSplitting=\"0.71712\">");
        mDWWriter.newLine();
		mDWWriter.write("<mainView=\"2D View\">");
        mDWWriter.newLine();
		mDWWriter.write("<mainViewCount=\"2\">");
        mDWWriter.newLine();
		mDWWriter.write("<mainViewDockInfo0=\"root\">");
        mDWWriter.newLine();
		mDWWriter.write("<mainViewDockInfo1=\"Table	center\">");
        mDWWriter.newLine();
		mDWWriter.write("<mainViewName0=\"Table\">");
        mDWWriter.newLine();
		mDWWriter.write("<mainViewName1=\"2D View\">");
        mDWWriter.newLine();
		mDWWriter.write("<mainViewType0=\"tableView\">");
        mDWWriter.newLine();
		mDWWriter.write("<mainViewType1=\"2Dview\">");
        mDWWriter.newLine();
		mDWWriter.write("<rightSplitting=\"0\">");
        mDWWriter.newLine();
		mDWWriter.write("<rowHeight_Table=\"80\">");
        mDWWriter.newLine();
        mDWWriter.write("<filter0=\"#category#\truleName\">");
        mDWWriter.newLine();
        mDWWriter.write("<connectionColumn_2D View=\"<connectAll>\">");
        mDWWriter.newLine();
        mDWWriter.write("<connectionLineWidth_2D View=\"0.17640000581741333\">");
        mDWWriter.newLine();
   		mDWWriter.write("<logarithmicView=\"totalStrain\">");
        mDWWriter.newLine();
		mDWWriter.write("<markersize_2D View=\"0.1936\">");
        mDWWriter.newLine();
		mDWWriter.write("<sizeAdaption_2D View=\"false\">");
        mDWWriter.newLine();
		mDWWriter.write("</datawarrior properties>");
        mDWWriter.newLine();
		}

	/**
	 * Tries to create one random conformer in current thread by once going through
	 * sequence of preparation, break-out, optimization and fine-tuning phases.
	 * Stops successfully early, if none of the atoms' strain exceeds MAX_ATOM_STRAIN
	 * and the conformer's total strain is below MAX_TOTAL_STRAIN.
	 * Returns false, if after going though all phases still one of these two conditions
	 * is not met.
	 * @param conformer receives coordinates of the conformer
	 * @return true if conformer with satisfactory low strain could be generated
	 */
	private boolean tryGenerateConformer(SelfOrganizedConformer conformer) {
		if (mMol.getAllAtoms() < 2)
			return true;

		if (!KEEP_INITIAL_COORDINATES)
			jumbleAtoms(conformer);

		mSkipRule[ConformationRule.RULE_TYPE_TORSION] = true;

		optimize(conformer, PREPARATION_CYCLES, STANDARD_CYCLE_FACTOR, 1.0);

		boolean done = false;

		if (mRuleCount[ConformationRule.RULE_TYPE_TORSION] != 0) {
			mSkipRule[ConformationRule.RULE_TYPE_TORSION] = false;
			done = optimize(conformer, PRE_OPTIMIZATION_CYCLES, STANDARD_CYCLE_FACTOR, 1.0);
			}

		for (int i=0; !done && i<MAX_BREAKOUT_ROUNDS; i++) {
			if (mThreadMaster != null && mThreadMaster.threadMustDie())
				break;

			if (jumbleStrainedAtoms(conformer) == 0)
				break;

			done = optimize(conformer, BREAKOUT_CYCLES, STANDARD_CYCLE_FACTOR, 1.0);
			}

//		if (!done && disableCollidingTorsionRules(conformer))
//			done = optimize(conformer, OPTIMIZATION_CYCLES, STANDARD_CYCLE_FACTOR, 1.0);

		if (!done)
			done = optimize(conformer, OPTIMIZATION_CYCLES, STANDARD_CYCLE_FACTOR, 1.0);

		if (!done)
			done = optimize(conformer, MINIMIZATION_CYCLES, STANDARD_CYCLE_FACTOR, MINIMIZATION_REDUCTION);

		return done;
		}

	public boolean optimize(SelfOrganizedConformer conformer, int cycles, double startFactor, double factorReduction) {
		int atomsSquare = mMol.getAllAtoms() * mMol.getAllAtoms();

		double k = Math.log(factorReduction)/cycles;
//double[] dummy_ = new double[mMol.getAllAtoms()];

		for (int outerCycle=0; outerCycle<cycles; outerCycle++) {
			if (mThreadMaster != null && mThreadMaster.threadMustDie())
				break;

			double cycleFactor = startFactor * Math.exp(-k*outerCycle);

			for (int innerCycle=0; innerCycle<atomsSquare; innerCycle++) {
				if (mThreadMaster != null && mThreadMaster.threadMustDie())
					break;

				ConformationRule rule = mRuleList.get((int)(mRandom.nextDouble() * mRuleList.size()));

/*				// Always use maximum strain constraint.
				ConformationRule rule = null;
				double maxStrain = -1;
				for (ConformationRule r:mRuleList) {
					if (r.isEnabled() && !mSkipRule[r.getRuleType()]) {
						double strain = r.addStrain(conformer, dummy_) / r.getAtomList().length;
						if (maxStrain < strain) {
							maxStrain = strain;
							rule = r;
							}
						}
					}
*/
				if (rule.isEnabled() && !mSkipRule[rule.getRuleType()]) {

//System.out.println("#1 rule:"+rule.toString());
					boolean conformerChanged = rule.apply(conformer, cycleFactor);
//System.out.println("atom 2  x:"+conformer.x[2]+" y:"+conformer.y[2]+" z:"+conformer.z[2]);
//System.out.println("atom 5  x:"+conformer.x[5]+" y:"+conformer.y[5]+" z:"+conformer.z[5]);
//System.out.println("atom 10 x:"+conformer.x[10]+" y:"+conformer.y[10]+" z:"+conformer.z[10]);
//System.out.println("#2");

					if (conformerChanged)
						conformer.invalidateStrain();

if (mDWWriter != null && conformerChanged) {
 try {
  double[] dummy = new double[mMol.getAllAtoms()];
  double s1 = (mLastDWConformer == null) ? 0 : rule.addStrain(mLastDWConformer, dummy);
  double s2 = rule.addStrain(conformer, dummy);
  writeStrains(conformer, rule, null, s1, s2);
 } catch (Exception e) { e.printStackTrace(); } }
					}
				}

			if (conformer.isAcceptable(mRuleList))
				return true;
			}

		return false;
		}

	private void writeStrains(SelfOrganizedConformer newConformer, ConformationRule rule, String stepName, double strain1, double strain2) throws Exception {
		newConformer.calculateStrain(mRuleList);
		double[] strain = new double[ConformationRule.RULE_NAME.length];
		double strainSum = 0f;
		for (int i=0; i<ConformationRule.RULE_NAME.length; i++) {
			strain[i] = newConformer.getRuleStrain(i);
			strainSum += strain[i];
			}

		double oldStrainSum = 0f;
		if (mDWStrain != null)
			for (int i=0; i<mDWStrain.length; i++)
				oldStrainSum += mDWStrain[i];

        String ruleName = (rule != null) ? ConformationRule.RULE_NAME[rule.getRuleType()] : stepName;

		StereoMolecule mol = mMol.getCompactCopy();
		for (int atom=0; atom<mol.getAllAtoms(); atom++) {
			if (mol.getAtomicNo(atom) == 1)
				mol.setAtomicNo(atom, 9);
			}

		String atoms = "";
	    if (rule != null) {
	    	int[] atomList = rule.getAtomList();
	        for (int i=0; i<atomList.length; i++) {
	            if (i!=0) atoms = atoms + ",";
	            atoms = atoms+atomList[i];
	            if (atomList[i] != -1)
	            	mol.setAtomicNo(atomList[i], 5);
	        	}
	    	}

		newConformer.toMolecule(mol);
		Canonizer newCanonizer = new Canonizer(mol);
		String newCoords = newCanonizer.getEncodedCoordinates();

		if (mLastDWConformer != null) {
			mLastDWConformer.toMolecule(mol);
			Canonizer oldCanonizer = new Canonizer(mol);
			String idcode = oldCanonizer.getIDCode();
			String oldCoords = oldCanonizer.getEncodedCoordinates();

			mDWWriter.write(idcode + "\t" + oldCoords + "\t" + newCoords + "\t" + mDWCycle + "\t" + ruleName + "\t" + atoms + "\t" + (rule != null ? rule.toString() : stepName));
			for (double s : strain)
				mDWWriter.write("\t" + s);
			mDWWriter.write("\t" + strainSum + "\t" + (oldStrainSum - strainSum) + "\t" + DoubleFormat.toString(strain1) + "\t" + DoubleFormat.toString(strain2) + "\t" + DoubleFormat.toString(strain1 - strain2));
			mDWWriter.newLine();
			mDWStrain = strain;
			mDWCycle++;
			}

		mLastDWConformer = new Conformer(newConformer);
	    }

	private void jumbleAtoms(SelfOrganizedConformer conformer) {
		double boxSize = 1.0 + 3.0 * Math.sqrt(mMol.getAllAtoms());
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			conformer.setX(atom, boxSize * mRandom.nextDouble() - boxSize / 2);
			conformer.setY(atom, boxSize * mRandom.nextDouble() - boxSize / 2);
			conformer.setZ(atom, boxSize * mRandom.nextDouble() - boxSize / 2);
			}

		conformer.invalidateStrain();
		}

	private int jumbleStrainedAtoms(SelfOrganizedConformer conformer) {
		conformer.calculateStrain(mRuleList);

		int atomCount = 0;
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			double atomStrain = conformer.getAtomStrain(atom);
			if (atomStrain > ATOM_FLAT_RING_BREAKOUT_STRAIN) {
				if (tryEscapeFromFlatRingTrap(conformer, atom)) {
					atomCount++;
					continue;
					}
				}
			if (atomStrain > ATOM_CAGE_BREAKOUT_STRAIN) {
				if (mDWWriter != null) {
					try {
						writeStrains(conformer, null, "escapeCage", atomStrain, Double.NaN);
						}
					catch (Exception e) { e.printStackTrace(); }
					}

				Coordinates c = conformer.getCoordinates(atom);
				c.add(BREAKOUT_DISTANCE * mRandom.nextDouble() - BREAKOUT_DISTANCE / 2,
					  BREAKOUT_DISTANCE * mRandom.nextDouble() - BREAKOUT_DISTANCE / 2,
					  BREAKOUT_DISTANCE * mRandom.nextDouble() - BREAKOUT_DISTANCE / 2);
				atomCount++;
				}
			}

		if (atomCount != 0)
			conformer.invalidateStrain();

		return atomCount;
		}

	/**
	 * Sometimes individual exocyclic atoms end up trapped inside a flat ring,
	 * because a plane rule combined with a distance rule stabilize the situation.
	 * This method checks, whether atom<br>
	 *     - has only one neighbour<br>
	 *     - is connected to a small ring atom<br>
	 *     - all angles atom-neighbour-otherRingAtom are below 90 degrees<br>
	 * If all conditions are true then atom is moved to the opposite side of the neighbour atom.
	 * @param conformer
	 * @param atom
	 * @return whether atom was moved from trapped state
	 */
	private boolean tryEscapeFromFlatRingTrap(SelfOrganizedConformer conformer, int atom) {
		if (mMol.getAllConnAtoms(atom) == 1) {
			int connAtom = mMol.getConnAtom(atom, 0);
			if (mMol.isSmallRingAtom(connAtom)) {
				Coordinates ca = conformer.getCoordinates(atom);
				Coordinates cc = conformer.getCoordinates(connAtom);
				Coordinates va = cc.subC(ca);
				for (int i=0; i<mMol.getConnAtoms(connAtom); i++) {
					int ringAtom = mMol.getConnAtom(connAtom, i);
					if (mMol.isRingAtom(ringAtom)) {
						Coordinates vr = cc.subC(conformer.getCoordinates(ringAtom));
						if (va.getAngle(vr) > Math.PI / 2)
							return false;
						}
					}

				ca.add(va).add(va);

				if (mDWWriter != null) {
					try {
						writeStrains(conformer, null, "escapeFlatRing", conformer.getAtomStrain(atom), Double.NaN);
						}
					catch (Exception e) { e.printStackTrace(); }
					}
				return true;
				}
			}
		return false;
		}

	public boolean disableCollidingTorsionRules(SelfOrganizedConformer conformer) {
		boolean found = false;
		conformer.calculateStrain(mRuleList);
		StereoMolecule mol = mMol;
		boolean[] isInvolvedAtom = new boolean[mol.getAllAtoms()];
		for (ConformationRule rule:mRuleList) {
			if (rule instanceof TorsionRule) {
				if (((TorsionRule)rule).disableIfColliding(conformer)) {
					int[] atom = rule.getAtomList();
					for (int i=1; i<=2; i++) {
					    for (int j=0; j<mol.getAllConnAtoms(atom[i]); j++) {
					    	int connAtom = mol.getConnAtom(atom[i], j);
					        if (connAtom != atom[3-i])
					        	isInvolvedAtom[connAtom] = true;
					    	}
						}
					found = true;
					}
				}
			}
		if (found)
			for (int atom=0; atom<mMol.getAllAtoms(); atom++)
				if (isInvolvedAtom[atom])
					conformer.getCoordinates(atom).add(0.6 * mRandom.nextDouble() - 0.3,
													  0.6 * mRandom.nextDouble() - 0.3,
													  0.6 * mRandom.nextDouble() - 0.3);
		return found;
		}

	public void disableTorsionRules() {
		for (ConformationRule rule:mRuleList)
			if (rule instanceof TorsionRule)
				rule.setEnabled(false);
		}

	public boolean disablePlaneRules() {
		boolean found = false;
		for (ConformationRule rule:mRuleList) {
			if (rule instanceof PlaneRule) {
				rule.setEnabled(false);
				found = true;
				}
			}
		return found;
		}

	public void enableTorsionRules() {
		for (ConformationRule rule:mRuleList)
			if (rule instanceof TorsionRule)
				rule.setEnabled(true);
		}
	}
