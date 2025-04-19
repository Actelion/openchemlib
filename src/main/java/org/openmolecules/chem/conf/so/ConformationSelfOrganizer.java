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
import com.actelion.research.chem.*;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.TorsionDescriptorHelper;
import com.actelion.research.util.DoubleFormat;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Random;

public class ConformationSelfOrganizer {
	private static final int PHASE_PREPARATION = 0;         // without ambiguous rules (i.e. torsion)
	private static final int PHASE_PRE_OPTIMIZATION = 1;    // with ambiguous rules before breakout
	private static final int PHASE_BREAKOUT = 2;
	private static final int PHASE_OPTIMIZATION = 3;
	private static final int PHASE_MINIMIZATION = 4;
	private static final int[] PHASE_CYCLES = { 40, 20, 20, 100, 20 };
	private static final double[] PHASE_FACTOR = { 40, 20, 20, 100, 40 };
	private static final String[] PHASE_NAME = { "preparation", "pre-optimization", "breakout", "optimization", "minimization" };

	private static final int    INITIAL_POOL_SIZE = 8;
	private static final int    MAX_CONFORMER_TRIES = 12;
	private static final int    MAX_BREAKOUT_ROUNDS = 3;
	private static final boolean PREFER_HIGH_STRAIN_RULES = false;
	private static final boolean INITIALLY_SKIP_TORSION_RULES = false;
	private static final int    CONSIDERED_PREVIOUS_CYCLE_COUNT = 10;    // must be much smaller than all the ...CYCLES above
	private static final double	STANDARD_CYCLE_FACTOR = 1.0;
	private static final double	MINIMIZATION_END_FACTOR = 0.01;
	private static final double	ATOM_ACCEPTABLE_STRAIN = 2.72; // A conformer is considered acceptable if all atom strains are below this value
	private static final double ATOM_FLAT_RING_BREAKOUT_STRAIN = 1000;    // TODO check
	private static final double ATOM_CAGE_BREAKOUT_STRAIN = 2000;  // TODO check Max allowed strain for 4-neighbour atoms
	private static final double TWIST_BOAT_ESCAPE_ANGLE = 0.6;  // angle to rotate ring member out of plane (from mid point between ring center and atom), after rotating it into the plane from the other side
	private static final int    TWIST_BOAT_ESCAPE_FREQUENCY = 10;   // in every tenth cycle we try escaping trapped twist boats
	private static final double	BREAKOUT_DISTANCE = 8.0;
	private static final double	MAX_POOL_STRAIN_DIF = 2.72;     // 2 * 1.36 kcal/mol, which is factor 100

public static boolean KEEP_INITIAL_COORDINATES = false;
public static boolean WRITE_DW_FILE = false;
private static final String DATAWARRIOR_DEBUG_FILE = "/home/thomas/data/debug/conformationSampler.dwar";
private static BufferedWriter mDWWriter;
private static Conformer mLastDWConformer;
private static int mDWCycle;
private String mDWMode;
private double[] mDWStrain; 	// TODO get rid of this section

	private final StereoMolecule mMol;
    private Random				mRandom;
    private int					mMaxConformers;
    private boolean				mPoolIsClosed;
	private final ArrayList<ConformationRule> mRuleList;
	private ArrayList<SelfOrganizedConformer> mConformerList;
	private double              mMinStrainInPool;
	private final int[]			mRuleCount;
	private final boolean[]		mSkipRule;
	private int[]				mRotatableBondForDescriptor;
	private ThreadMaster        mThreadMaster;
	private long                mStopMillis;

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
		    mMol.removeExplicitHydrogens(true);
		mMol.ensureHelperArrays(Molecule.cHelperParities);

		mRuleList = new ArrayList<>();
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
	 * @param millis time point as system millis after which to gracefully stop self organization even if not successful
	 */
	public void setStopTime(long millis) {
		mStopMillis = millis;
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
		SelfOrganizedConformer bestConformer = null;
		for (int i=0; i<MAX_CONFORMER_TRIES; i++) {
			if (mustStop())
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

        mConformerList = new ArrayList<>();
		mMinStrainInPool = Double.MAX_VALUE;
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
		for (int i=0; i<tryCount && mConformerList.size()<finalPoolSize && !mustStop(); i++) {
			if (conformer == null)
				conformer = new SelfOrganizedConformer(mMol);
			if (tryGenerateConformer(conformer)) {
				// if optimization is successful, we have a good conformer
				if (addIfNewOrReplaceIfBetter(conformer))
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

	/**
	 * @param conformer
	 * @return true if the conformer was added or replaced another one
	 */
	private boolean addIfNewOrReplaceIfBetter(SelfOrganizedConformer conformer) {
		if (mRotatableBondForDescriptor == null)
			mRotatableBondForDescriptor = TorsionDescriptorHelper.findRotatableBonds(getMolecule());

		if (conformer.getTotalStrain() > mMinStrainInPool + MAX_POOL_STRAIN_DIF)
			return false;

		conformer.calculateDescriptor(mRotatableBondForDescriptor);

		for (int i=mConformerList.size()-1; i>=0; i--) {
			SelfOrganizedConformer soc = mConformerList.get(i);
			if (conformer.equals(soc)) {
				if (soc.isWorseThan(conformer)) {
					mConformerList.remove(i);
					mConformerList.add(conformer);
					if (mMinStrainInPool > conformer.getTotalStrain())
						mMinStrainInPool = conformer.getTotalStrain();
					return true;
					}

				return false;
				}
			}

		mConformerList.add(conformer);
		if (mMinStrainInPool > conformer.getTotalStrain()) {
			mMinStrainInPool = conformer.getTotalStrain();

			for (int i=mConformerList.size()-1; i>=0; i--) {
				SelfOrganizedConformer soc = mConformerList.get(i);
				if (soc.getTotalStrain() > mMinStrainInPool + MAX_POOL_STRAIN_DIF)
					mConformerList.remove(i);
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
	 * Internally, TorsionDescriptors are compared to decide whether a conformer is new.
	 * These TorsionDescriptors are created on the fly using the default method to determine all
	 * considered rotatable bonds. Bonds that contain marked atoms are excluded from being considered
	 * rotatable.
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

	public void printStrainDetails(SelfOrganizedConformer conformer) {
		System.out.println("############ new conformer ###############");
		for (ConformationRule rule:mRuleList)
			if (rule instanceof DistanceRule)
				((DistanceRule)rule).printStrain(conformer);
		}

	public static void writeDWFileStart() {
		mDWCycle = 0;
		mLastDWConformer = null;

		try {
			mDWWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(DATAWARRIOR_DEBUG_FILE), StandardCharsets.UTF_8));
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
			mDWWriter.write("Structure\tbefore\tafter\tcycle\tmode\truleName\truleAtoms\truleDetail");
			for (int i = 0; i<ConformationRule.RULE_NAME.length; i++)
				mDWWriter.write("\t" + ConformationRule.RULE_NAME[i]);
			mDWWriter.write("\ttotalStrain\tstrainGain\truleStrainBefore\truleStrainAfter\truleStrainGain\tatomStrain");
			mDWWriter.newLine();
			}
		catch (IOException ioe) {}
        }

	public static void writeDWFileEnd() {
		try {
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
			mDWWriter.write("<mainViewCount=\"3\">");
	        mDWWriter.newLine();
			mDWWriter.write("<mainViewDockInfo0=\"root\">");
	        mDWWriter.newLine();
			mDWWriter.write("<mainViewDockInfo1=\"Table	center\">");
			mDWWriter.newLine();
			mDWWriter.write("<mainViewDockInfo2=\"2D View\tbottom\t0.501\">");
	        mDWWriter.newLine();
			mDWWriter.write("<mainViewName0=\"Table\">");
	        mDWWriter.newLine();
			mDWWriter.write("<mainViewName1=\"2D View\">");
	        mDWWriter.newLine();
			mDWWriter.write("<mainViewName2=\"Form View\">");
			mDWWriter.newLine();
			mDWWriter.write("<mainViewType0=\"tableView\">");
	        mDWWriter.newLine();
			mDWWriter.write("<mainViewType1=\"2Dview\">");
	        mDWWriter.newLine();
			mDWWriter.write("<mainViewType2=\"formView\">");
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
			mDWWriter.write("<shapeColumn_2D View=\"mode\">");
			mDWWriter.newLine();
			mDWWriter.write("<sizeAdaption_2D View=\"false\">");
	        mDWWriter.newLine();
			mDWWriter.write("<autoZoomFactor_2D View=\""+mDWCycle/120+".0;0.0\">");
			mDWWriter.newLine();
			mDWWriter.write("<formLayout_Form View=\"TableLayout,7,11.0,-1.0,11.0,-1.0,11.0,-1.0,11.0,27,7.0,-1.0,7.0,-1.0,7.0,-1.0,7.0,-1.0,7.0,-1.0,7.0,-1.0,7.0,-1.0,7.0,-1.0,7.0,-1.0,7.0,-1.0,7.0,-1.0,7.0,-1.0,7.0,-1.0,7.0\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectCount_Form View=\"19\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_0=\"before\tstructure3D\t1, 5, 1, 21, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_1=\"after\tstructure3D\t3, 5, 3, 21, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_10=\"torsion\ttextLine\t5, 21, 5, 21, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_11=\"stereo\ttextLine\t5, 23, 5, 23, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_12=\"binap\ttextLine\t5, 25, 5, 25, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_13=\"totalStrain\ttextLine\t5, 13, 5, 13, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_14=\"strainGain\ttextLine\t5, 1, 5, 1, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_15=\"ruleStrainBefore\ttextLine\t5, 3, 5, 3, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_16=\"ruleStrainAfter\ttextLine\t5, 5, 5, 5, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_17=\"ruleStrainGain\ttextLine\t5, 7, 5, 7, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_18=\"atomStrain\ttextLine\t1, 23, 3, 25, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_2=\"cycle\ttextLine\t1, 1, 1, 1, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_3=\"mode\ttextLine\t3, 1, 3, 1, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_4=\"ruleName\ttextLine\t1, 3, 1, 3, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_5=\"ruleAtoms\ttextLine\t3, 3, 3, 3, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_6=\"ruleDetail\ttextLine\t5, 9, 5, 11, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_7=\"distance\ttextLine\t5, 15, 5, 15, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_8=\"plane\ttextLine\t5, 17, 5, 17, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("<formObjectInfo_Form View_9=\"line\ttextLine\t5, 19, 5, 19, full, full\">");
			mDWWriter.newLine();
			mDWWriter.write("</datawarrior properties>");
	        mDWWriter.newLine();
			mDWWriter.close();
			mDWWriter = null;
			}
		catch (IOException ioe) {}
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

		// Make sure no rule were disabled in previous run.
		for (ConformationRule rule:mRuleList)
			rule.setEnabled(true);
		for (int i=0; i<mSkipRule.length; i++)
			mSkipRule[i] = false;

		if (!KEEP_INITIAL_COORDINATES)
			jumbleAtoms(conformer);

		if (INITIALLY_SKIP_TORSION_RULES)
			mSkipRule[ConformationRule.RULE_TYPE_TORSION] = true;

		SelfOrganizedConformer bestConformer = new SelfOrganizedConformer(conformer);

		optimize(conformer, bestConformer, PHASE_PREPARATION);

		if (INITIALLY_SKIP_TORSION_RULES
		 && mRuleCount[ConformationRule.RULE_TYPE_TORSION] != 0) {
			mSkipRule[ConformationRule.RULE_TYPE_TORSION] = false;

			optimize(conformer, bestConformer, PHASE_PRE_OPTIMIZATION);
			}

		for (int i=0; containsTrappedAtom(conformer) && i<MAX_BREAKOUT_ROUNDS && !mustStop(); i++) {
			if (escapeFromTrappedStates(conformer) == 0)
				break;
			optimize(conformer, bestConformer, PHASE_BREAKOUT);
			}

		optimize(conformer, bestConformer, PHASE_OPTIMIZATION);
		optimize(conformer, bestConformer, PHASE_MINIMIZATION);

		return isAcceptable(conformer);
		}

	private void optimize(SelfOrganizedConformer conformer, SelfOrganizedConformer bestConformer, int phase) {
		int cycles = PHASE_CYCLES[phase];
		double cycleFactor = STANDARD_CYCLE_FACTOR;
		double reductionFactor = (phase == PHASE_MINIMIZATION) ? Math.exp(Math.log(MINIMIZATION_END_FACTOR)/cycles) : 1.0;
mDWMode = PHASE_NAME[phase];

		int atomsSquare = mMol.getAllAtoms() * mMol.getAllAtoms();

		double[] previousTotalStrains = new double[CONSIDERED_PREVIOUS_CYCLE_COUNT];
		int previousTotalStrainIndex = 0;

		for (int outerCycle=0; outerCycle<cycles && !mustStop(); outerCycle++) {
			if (phase != PHASE_PREPARATION && (outerCycle % TWIST_BOAT_ESCAPE_FREQUENCY) == 0)
				tryEscapeTwistBoats(conformer);

			for (int innerCycle=0; innerCycle<atomsSquare && !mustStop(); innerCycle++) {
				ConformationRule rule = null;
				if (PREFER_HIGH_STRAIN_RULES) {
					// Select the rule based on a weighted random algorithm preferring rules with larger strain
					double[] ruleStrain = new double[1+mRuleList.size()];
					int index = 0;
					for (ConformationRule r:mRuleList) {
						if (r.isEnabled() && !mSkipRule[r.getRuleType()]) {
							index++;
							ruleStrain[index] = ruleStrain[index - 1] + r.addStrain(conformer, null);
							}
						}
					double random = ruleStrain[index] * mRandom.nextDouble();
					index = 0;
					for (ConformationRule r:mRuleList) {
						if (r.isEnabled() && !mSkipRule[r.getRuleType()]) {
							index++;
							if (random < ruleStrain[index]) {
								rule = r;
								break;
								}
							}
						}
					}
				else {
					// Purely random rule selection.
					rule = mRuleList.get((int)(mRandom.nextDouble() * mRuleList.size()));
					}

				// Always use maximum strain constraint.
/*				ConformationRule rule = null;
				double maxStrain = -1;
				for (ConformationRule r:mRuleList) {
					if (r.isEnabled() && !mSkipRule[r.getRuleType()]) {
						double strain = r.addStrain(conformer, null);
						if (maxStrain < strain) {
							maxStrain = strain;
							rule = r;
							}
						}
					}*/

				if (rule.isEnabled() && !mSkipRule[rule.getRuleType()]) {
					boolean conformerChanged = rule.apply(conformer, cycleFactor);

					if (conformerChanged)
						conformer.invalidateStrain();

if (mDWWriter != null && conformerChanged) {
 try {
  double s1 = (mLastDWConformer == null) ? 0 : rule.addStrain(mLastDWConformer, null);
  double s2 = rule.addStrain(conformer, null);
  writeStrains(conformer, rule, null, s1, s2);
 } catch (Exception e) { e.printStackTrace(); } }
					}
				}

			// we break if the current total strain is not lower than the average from the previous X cycles.
			conformer.calculateStrain(mRuleList);

			previousTotalStrains[previousTotalStrainIndex++] = conformer.getTotalStrain();
			if (previousTotalStrainIndex == CONSIDERED_PREVIOUS_CYCLE_COUNT)
				previousTotalStrainIndex = 0;

			if (outerCycle > CONSIDERED_PREVIOUS_CYCLE_COUNT) {
				double averagePreviousStrain = 0.0;
				for (double strain:previousTotalStrains)
					averagePreviousStrain += strain;
				averagePreviousStrain /= CONSIDERED_PREVIOUS_CYCLE_COUNT;

				if (conformer.getTotalStrain() > averagePreviousStrain)
					break;
				}

			if (bestConformer.getTotalStrain() > conformer.getTotalStrain())
				bestConformer.copyFrom(conformer);

			cycleFactor *= reductionFactor;
			}

		if (conformer.isWorseThan(bestConformer))
			conformer.copyFrom(bestConformer);
		}

	private boolean isAcceptable(SelfOrganizedConformer conformer) {
		for (int atom=0; atom<conformer.getMolecule().getAllAtoms(); atom++)
			if (conformer.getAtomStrain(atom) > ATOM_ACCEPTABLE_STRAIN)
				return false;

		return true;
		}

	private boolean containsTrappedAtom(SelfOrganizedConformer conformer) {
		for (int atom=0; atom<conformer.getMolecule().getAllAtoms(); atom++)
			if (conformer.getAtomStrain(atom) > ATOM_FLAT_RING_BREAKOUT_STRAIN
			 || conformer.getAtomStrain(atom) > ATOM_CAGE_BREAKOUT_STRAIN)
				return true;

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

			mDWWriter.write(idcode + "\t" + oldCoords + "\t" + newCoords + "\t" + mDWCycle + "\t" + mDWMode + "\t" + ruleName + "\t" + atoms + "\t" + (rule != null ? rule.toString() : stepName));
			for (double s : strain)
				mDWWriter.write("\t" + DoubleFormat.toString(s));
			mDWWriter.write("\t" + DoubleFormat.toString(strainSum) + "\t" + DoubleFormat.toString(oldStrainSum - strainSum) + "\t" + DoubleFormat.toString(strain1) + "\t" + DoubleFormat.toString(strain2) + "\t" + DoubleFormat.toString(strain1 - strain2));
			for (int atom=0; atom<mol.getAllAtoms(); atom++) {
				mDWWriter.write(atom == 0 ? "\t" : " ");
				mDWWriter.write(atom+":"+DoubleFormat.toString(newConformer.getAtomStrain(atom), 2));
				}
			mDWWriter.newLine();
			mDWStrain = strain;
			mDWCycle++;
			}

		mLastDWConformer = new Conformer(newConformer);
	    }

	private void jumbleAtoms(SelfOrganizedConformer conformer) {
		double boxSize = 2.0 + 2.0 * Math.sqrt(mMol.getAllAtoms());
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			if (mMol.getAllConnAtoms(atom) != 1) {
				conformer.setX(atom, boxSize * mRandom.nextDouble() - boxSize / 2);
				conformer.setY(atom, boxSize * mRandom.nextDouble() - boxSize / 2);
				conformer.setZ(atom, boxSize * mRandom.nextDouble() - boxSize / 2);
				}
			}
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			if (mMol.getAllConnAtoms(atom) == 1) {
				int connAtom = mMol.getConnAtom(atom, 0);
				conformer.setX(atom, conformer.getX(connAtom) + 4 * mRandom.nextDouble() - 2);
				conformer.setY(atom, conformer.getY(connAtom) + 4 * mRandom.nextDouble() - 2);
				conformer.setZ(atom, conformer.getZ(connAtom) + 4 * mRandom.nextDouble() - 2);
				}
			}

		conformer.invalidateStrain();
		}

	private int escapeFromTrappedStates(SelfOrganizedConformer conformer) {
		conformer.calculateStrain(mRuleList);

		int atomCount = 0;
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			double atomStrain = conformer.getAtomStrain(atom);
			if (atomStrain > ATOM_FLAT_RING_BREAKOUT_STRAIN
			 && tryEscapeFromFlatRingTrap(conformer, atom))
					atomCount++;
			}

		if (atomCount == 0) {
			int neighbourCount = 16;
			for (int atom=0; atom<mMol.getAllAtoms(); atom++)
				if (neighbourCount > mMol.getAllConnAtoms(atom)
				 && conformer.getAtomStrain(atom) > maxCageBreakoutStrain(atom))
					neighbourCount = mMol.getAllConnAtoms(atom);

			if (neighbourCount != 16) {
				for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
					if (neighbourCount == mMol.getAllConnAtoms(atom)
					 && conformer.getAtomStrain(atom) > maxCageBreakoutStrain(atom)) {
//System.out.println("escape "+neighbourCount+" neighbours");
//System.out.print("strains: "); for (int i=0; i<mMol.getAllAtoms(); i++) System.out.print(" "+i+":"+conformer.getAtomStrain(i)); System.out.println();
						if (mDWWriter != null) {
							try {
								writeStrains(conformer, null, "escapeCage", conformer.getAtomStrain(atom), Double.NaN);
								}
							catch (Exception e) { e.printStackTrace(); }
							}

						Coordinates c = conformer.getCoordinates(atom);
						if (mMol.getAllConnAtoms(atom) == 1) {
							Coordinates cn = conformer.getCoordinates(mMol.getConnAtom(atom, 0));
							c.add(2.0 * (cn.x - c.x), 2.0 * (cn.y - c.y), 2.0 * (cn.z - c.z));
							}
						else {
							double distance = (mMol.getAllConnAtoms(atom) == 0) ? 2.0 * BREAKOUT_DISTANCE : BREAKOUT_DISTANCE;
							c.add(distance * mRandom.nextDouble() - distance / 2,
								  distance * mRandom.nextDouble() - distance / 2,
								  distance * mRandom.nextDouble() - distance / 2);
							}
						atomCount++;
						}
					}
				}
			}

		if (atomCount != 0)
			conformer.invalidateStrain();

		return atomCount;
		}

	private double maxCageBreakoutStrain(int atom) {
		int connAtoms = Math.min(4, mMol.getAllConnAtoms(atom));
		return ATOM_CAGE_BREAKOUT_STRAIN / 5 + connAtoms * ATOM_CAGE_BREAKOUT_STRAIN / 5;
		}

	/**
	 * Sometimes individual exo-cyclic atoms end up trapped inside a flat ring,
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

	private boolean tryEscapeTwistBoats(SelfOrganizedConformer conformer) {
		boolean changeDone = false;
		RingCollection ringSet = mMol.getRingSet();
		for (int r=0; r<ringSet.getSize(); r++) {
			int[] ringAtom = ringSet.getRingAtoms(r);
			if (ringAtom.length == 6 && !ringSet.isDelocalized(r)) {
				boolean isCandidate = true;
				for (int rb:ringSet.getRingBonds(r)) {
					if (mMol.getBondOrder(rb) != 1 || mMol.isDelocalizedBond(rb)) {
						isCandidate = false;
						break;
						}
					}
				if (isCandidate) {
					Coordinates cog = new Coordinates();
					Coordinates n = new Coordinates();
					double[][] coords = new double[6][3];
					ConformationRule.calculateNearestPlane(conformer, ringAtom, cog, n, coords);

					double[] distance = new double[ringAtom.length];
					int sideMatchCount = 0;
//					int pattern = 0;
					for (int i=0; i<ringAtom.length; i++) {
						distance[i] = -(n.x * coords[i][0] + n.y * coords[i][1] + n.z * coords[i][2]);
//						pattern += (distance[i] < 0) ? 1 << i : 0;
						if ((distance[i] < 0) ^ ((i & 1) == 1))
							sideMatchCount++;
						}

					// We don't want to escape from boat conformations!
//					if (pattern == 9 || pattern == 18 || pattern == 36
//					 || pattern == 54 || pattern == 45 || pattern == 27)
//						continue;

					if (sideMatchCount != 0 && sideMatchCount != 6) {
						for (int i=0; i<ringAtom.length; i++) {
							if ((distance[i] < 0) ^ ((i & 1) == 1) ^ (sideMatchCount >= 3)) {
								int[] notAtom = new int[2];
								notAtom[0] = ringAtom[i == 0 ? 5 : i-1];
								notAtom[1] = ringAtom[i == 5 ? 0 : i+1];
								Coordinates pAtom = conformer.getCoordinates(ringAtom[i]);
								Coordinates axis = conformer.getCoordinates(ringAtom[i]).subC(cog).cross(n).unit();
								double angle = TWIST_BOAT_ESCAPE_ANGLE + Math.asin(Math.abs(distance[i])/pAtom.distance(cog));
								double theta = (distance[i] < 0) ? angle : -angle;
								Coordinates center = new Coordinates(pAtom).center(cog);
								boolean isBridgedRingAtom = false;
								if (mMol.getConnAtoms(ringAtom[i]) > 2) {
									for (int j=0; j<mMol.getConnAtoms(ringAtom[i]); j++) {
										int connAtom = mMol.getConnAtom(ringAtom[i], j);
										if (connAtom != notAtom[0]
										 && connAtom != notAtom[1]
										 && mMol.isRingBond(mMol.getConnBond(ringAtom[i], j))) {
											isBridgedRingAtom = true;
											break;
										}
									}
								}
								if (!isBridgedRingAtom)
									ConformationRule.rotateGroup(conformer, ringAtom[i], notAtom, center, axis, theta);
								}
							}
						changeDone = true;
						}
					}
				}
			}

		if (changeDone) {
			double strain1 = conformer.getTotalStrain();

			conformer.invalidateStrain();
			conformer.calculateStrain(mRuleList);

			if (mDWWriter != null) {
				try {
					writeStrains(conformer, null, "escapeTwistBoat", strain1, conformer.getTotalStrain());
				}
				catch (Exception e) { e.printStackTrace(); }
			}
		}

			return changeDone;
		}

	private boolean mustStop() {
		return (mThreadMaster != null && mThreadMaster.threadMustDie())
			|| (mStopMillis != 0 && System.currentTimeMillis() > mStopMillis);
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
