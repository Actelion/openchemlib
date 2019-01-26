/*
 * @(#)TorsionRule.java
 *
 * Copyright 2013 openmolecules.org, Inc. All Rights Reserved.
 *
 * NOTICE: All information contained herein is, and remains the property
 * of openmolecules.org.  The intellectual and technical concepts contained
 * herein are proprietary to openmolecules.org.
 * Actelion Pharmaceuticals Ltd. is granted a non-exclusive, non-transferable
 * and timely unlimited usage license.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.so;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.*;

import java.util.ArrayList;

public class TorsionRule extends ConformationRule {
	final static double COLLIDING_ATOM_STRAIN = 0.1;

	private int		    mSmallerSubstituentIndex;
	private int[]	    mAtomToRotate;
	private short[] 	mTorsion,mFrequency;
	private short[][]   mRange;

	public TorsionRule(short[] torsion, short[] frequency, short[][] range, int[] torsionAtom, int[] atomToRotate, int smallerSubstituentIndex) {
		super(torsionAtom);
		mTorsion = torsion;
		mFrequency = frequency;
		mRange = range;
		mAtomToRotate = atomToRotate;
		mSmallerSubstituentIndex = smallerSubstituentIndex;
		}

	@Override
	public int getRuleType() {
		return RULE_TYPE_TORSION;
		}

    public static void calculateRules(ArrayList<ConformationRule> ruleList, StereoMolecule mol) {
		TorsionDB.initialize(TorsionDB.MODE_ANGLES);

		boolean[] isRotatableBond = new boolean[mol.getAllBonds()];
		TorsionDB.findRotatableBonds(mol, false, isRotatableBond);
		for (int bond=0; bond<mol.getBonds(); bond++) {
			if (isRotatableBond[bond]) {
				int[] torsionAtom = new int[4];
				TorsionDetail torsionDetail = new TorsionDetail();
			    String torsionID = TorsionDB.getTorsionID(mol, bond, torsionAtom, torsionDetail);
			    if (torsionID != null && !conflictWithPlaneRules(torsionAtom, ruleList)) {
			        short[] torsion = TorsionDB.getTorsions(torsionID);
			        short[] frequency = TorsionDB.getTorsionFrequencies(torsionID);
				    short[][] range = TorsionDB.getTorsionRanges(torsionID);

					if (torsion == null) {
						TorsionPrediction prediction = new TorsionPrediction(mol, torsionAtom);
						torsion = prediction.getTorsions();
						frequency = prediction.getTorsionFrequencies();
						range = prediction.getTorsionRanges();
						}

                    if (torsion != null) {
/*
System.out.print("torsionID:"+torsionID+" torsions:"+torsion[0]);
for (int i=1; i<torsion.length; i++)
System.out.print(","+torsion[i]);
System.out.println();
*/
					    int[] atomToRotate = null;
					    int smallerSubstituentIndex = 0;
					    if (!mol.isRingBond(bond)) {
				    		boolean[][] isMemberAtom = new boolean[2][mol.getAllAtoms()];
				    		int[] count = new int[2];
					    	for (int i=0; i<2; i++)
					    		count[i] = mol.getSubstituent(torsionDetail.getRearAtom(i),
					    				torsionDetail.getCentralAtom(i), isMemberAtom[i], null, null);

					    	smallerSubstituentIndex = (count[0] < count[1]) ? 0 : 1;
					    	atomToRotate = new int[count[smallerSubstituentIndex]];
					    	int index = 0;
					    	for (int a=0; a<mol.getAllAtoms(); a++)
					    		if (isMemberAtom[smallerSubstituentIndex][a])
					    			atomToRotate[index++] = a;
					        }
					    ruleList.add(new TorsionRule(torsion, frequency, range, torsionAtom, atomToRotate, smallerSubstituentIndex));
                        }
				    }
				}
			}

		// In the CSD/COD torsion tables atom sequences cover non-H atoms only.
		// Often this is not a problem, because in sp3 chains distance rules cause
		// hydrogen gauche positioning. If we have an -OH or -NH2 connected to a pi-system,
		// the hydrogen needs to be in the pi-plane. We add artificial torsion rules
		// to rotate the hydrogen of ?=?-X-H into the pi-plane.
		for (int bond=0; bond<mol.getBonds(); bond++) {
			if (mol.getBondType(bond) == Molecule.cBondTypeSingle	// does not include up/down-bonds
			 && !mol.isRingBond(bond)) {
				for (int i=0; i<2; i++) {
					int heteroAtom = mol.getBondAtom(i, bond);
					int atomicNo = mol.getAtomicNo(heteroAtom);
					if (atomicNo > 6
					 && mol.getConnAtoms(heteroAtom) == 1) {
						int hCount = mol.getAllConnAtoms(heteroAtom) - 1;
						if (hCount != 0) {	// hetero atom with hydrogen(s) and one non-H neighbor
							int piAtom = mol.getBondAtom(1-i, bond);
							if (mol.getAtomPi(piAtom) == 1) {
								int piNeighbour = -1;
								for (int j=0; j<mol.getConnAtoms(piAtom); j++)
									if (mol.getConnBondOrder(piAtom, j) == 2)
										piNeighbour = mol.getConnAtom(piAtom, j);

								int[] torsionAtom = new int[4];
								torsionAtom[0] = piNeighbour;
								torsionAtom[1] = piAtom;
								torsionAtom[2] = heteroAtom;
								torsionAtom[3] = mol.getConnAtom(heteroAtom, 1);	// first H-neighbor

								int[] atomToRotate = new int[hCount];
								for (int j=0; j<hCount; j++)
									atomToRotate[j] = mol.getConnAtom(heteroAtom, mol.getConnAtoms(heteroAtom) + j);

								// carboxylic acid or similar are only in Z-conformation
								int stateCount = (!mol.isAromaticAtom(piAtom) && hCount == 1) ? 1 : 2;

								short[] torsion = new short[stateCount];
								short[] frequency = new short[stateCount];
								short[][] range = new short[stateCount][2];

								if (!mol.isAromaticAtom(piAtom) && hCount == 1) {
									torsion[0] = 0;

									frequency[0] = 100;

									range[0][0] = -15;
									range[0][1] = 15;
									}
								else {
									torsion[0] = 0;
									torsion[1] = 180;

									frequency[0] = 50;
									frequency[1] = 50;

									range[0][0] = -15;
									range[0][1] = 15;
									range[1][0] = 165;
									range[1][1] = 195;
									}

								ruleList.add(new TorsionRule(torsion, frequency, range, torsionAtom, atomToRotate, 1));
								}
							}
						}
					}
				}
			}
    	}

	private static boolean conflictWithPlaneRules(int[] torsionAtom, ArrayList<ConformationRule> ruleList) {
		for (ConformationRule rule:ruleList) {
			if (rule.getRuleType() == ConformationRule.RULE_TYPE_PLANE) {
				int[] planeAtom = rule.getAtomList();
				for (int ta:torsionAtom) {
					if (ta == -1)
						return false;

					boolean found = false;
					for (int pa : planeAtom) {
						if (ta == pa) {
							found = true;
							break;
							}
						}
					if (!found)
						return false;
					}

				return true;
				}
			}
		return false;
		}

	@Override
	public boolean apply(Conformer conformer, double cycleFactor) {
	    double currentTorsion = TorsionDB.calculateTorsionExtended(conformer, mAtom);
	    if (Double.isNaN(currentTorsion))
	    	return false;

	    if (currentTorsion < 0.0)
	        currentTorsion += 2 * Math.PI;

		int index = findApplicableTorsionIndex(currentTorsion);

		double severity = getSeverity(currentTorsion, index);
		if (severity == 0.0)
			return false;

	    double angleCorrection = Math.PI * mTorsion[index] / 180.0 - currentTorsion;
	    if (Math.abs(angleCorrection) > Math.PI)
	        angleCorrection = (angleCorrection < 0) ? 2*Math.PI + angleCorrection : angleCorrection - 2*Math.PI;

	    if (Math.abs(angleCorrection) < 0.001 * Math.PI)
	    	return false;

		Coordinates unit = conformer.getCoordinates(mAtom[2]).subC(conformer.getCoordinates(mAtom[1]));
		unit.unit();

		angleCorrection *= cycleFactor * severity;

	    StereoMolecule mol = conformer.getMolecule();

	    if (mAtomToRotate != null) {	// rotate smaller side of the molecule
    	    double rotation = (mSmallerSubstituentIndex == 0) ? -angleCorrection : angleCorrection;
            for (int a:mAtomToRotate)
            	rotateAtom(conformer, a, mAtom[1], unit, rotation);
	        }
	    else {	// rotate first and second atom shell from bond atoms; reduce angle for second shell atoms and if one side is more rigid
	    	int bond = mol.getBond(mAtom[1], mAtom[2]);
	    	boolean isFiveMemberedRing = (bond != -1 && mol.getBondRingSize(bond) <= 5);
	    	for (int i=1; i<=2; i++) {
			    double factor = (i==1 ? -2.0 : 2.0) * mol.getAtomRingBondCount(mAtom[i]);
	    	    for (int j=0; j<mol.getAllConnAtoms(mAtom[i]); j++) {
	    	    	int firstConn = mol.getConnAtom(mAtom[i], j);
	    	        if (firstConn != mAtom[3-i]) {
	    	            rotateGroup(conformer, firstConn, mAtom[i], unit, angleCorrection / factor);
	    	            if (!isFiveMemberedRing) {
		    	    	    for (int k=0; k<mol.getConnAtoms(firstConn); k++) {
		    	    	    	int secondConn = mol.getConnAtom(firstConn, k);
		    	    	        if (secondConn != firstConn && mol.getConnAtoms(secondConn) != 1) {
		    	    	            rotateGroup(conformer, secondConn, mAtom[i], unit, angleCorrection / (4f*factor));
		    	    	        	}
		    	    	    	}
	    	    	    	}
	    	        	}
	    	    	}
	    		}
	        }

/*
if (mol.getConnAtoms(mAtom[1])+mol.getConnAtoms(mAtom[2])==4) {
double after = TorsionDB.calculateTorsion(conformer, mAtom);
if (after < 0) after += 2 * Math.PI;
double improvement = Math.abs(Molecule.getAngleDif(currentTorsion, optTorsion))-Math.abs(Molecule.getAngleDif(after, optTorsion));
System.out.println((mAtomToRotate==null?"ring":"!ring")+" before:"+currentTorsion+" after:"+after+" wanted:"+optTorsion+" correction:"+angleCorrection+" factor:"+cycleFactor+" atoms:"+mAtom[0]+","+mAtom[1]+","+mAtom[2]+","+mAtom[3]+" gain:"+improvement+(improvement<0f?" WORSE!!!":""));
}
*/
	    return true;
		}

	private int findApplicableTorsionIndex(double torsion) {
		int index = -1;
		double minDif = Double.MAX_VALUE;
		for (int i=0; i<mTorsion.length; i++) {
			double t = Math.PI * mTorsion[i] / 180.0;
			double dif = Math.abs(torsion - t);
			if (dif > Math.PI)
				dif = 2*Math.PI - dif;
			dif /= (10+Math.sqrt(mFrequency[i]));	// normalize somewhat by frequency
			if (minDif > dif) {
				minDif = dif;
				index = i;
				}
			}
		return index;
		}

	/**
	 * Calculates a severity factor for the torsion depending on how
	 * far it is from the optimum and whether it is still in the range.
	 * @param angle from 0 to 2pi
	 * @param index of the torsion
	 * @return severity from 0 to 1
	 */
	private double getSeverity(double angle, int index) {
		double range0 = mRange[index][0] * Math.PI / 180;
		double range1 = mRange[index][1] * Math.PI / 180;
		double torsion = mTorsion[index] * Math.PI / 180;
		if (angle < torsion - Math.PI)
			angle += 2 * Math.PI;
		else if (angle > torsion + Math.PI)
			angle -= 2 * Math.PI;

		double dif = (torsion - angle) / (torsion - ((angle<torsion) ? range0 : range1));
		if (dif < 1.0) // in range
			return 0;

		if (dif < 2.0) {
			dif -= 1.0;
			return dif * dif;
			}

		return 1.0;
		}

	/**
	 * Rotate atom and all of its exclusive neighbor atoms.
	 * @param conformer
	 * @param atom
	 * @param refAtom
	 * @param unit
	 * @param theta
	 */
	private void rotateGroup(Conformer conformer, int atom, int refAtom, Coordinates unit, double theta) {
		rotateAtom(conformer, atom, refAtom, unit, theta);
	    StereoMolecule mol = conformer.getMolecule();
        for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
        	int connAtom = mol.getConnAtom(atom, i);
        	if (mol.getAllConnAtoms(connAtom) == 1)
        		rotateAtom(conformer, connAtom, refAtom, unit, theta);
        	}
		}

	private void rotateAtom(Conformer conformer, int atom, int refAtom, Coordinates unit, double theta) {
        double x = unit.x;
        double y = unit.y;
        double z = unit.z;
        double c = Math.cos(theta);
        double s = Math.sin(theta);
        double t = 1-c;
        double mx = conformer.getX(atom) - conformer.getX(refAtom);
        double my = conformer.getY(atom) - conformer.getY(refAtom);
        double mz = conformer.getZ(atom) - conformer.getZ(refAtom);
        conformer.setX(atom, conformer.getX(refAtom) + (t*x*x+c)*mx + (t*x*y+s*z)*my + (t*x*z-s*y)*mz);
        conformer.setY(atom, conformer.getY(refAtom) + (t*x*y-s*z)*mx + (t*y*y+c)*my + (t*y*z+s*x)*mz);
        conformer.setZ(atom, conformer.getZ(refAtom) + (t*x*z+s*y)*mx + (t*z*y-s*x)*my + (t*z*z+c)*mz);
	    }

	@Override
	public double addStrain(Conformer conformer, double[] atomStrain) {
		double angle = TorsionDB.calculateTorsionExtended(conformer, mAtom);
		if (Double.isNaN(angle))
			return 0.0;

		if (angle < 0.0)
			angle += 2 * Math.PI;

		double torsion = 180 * angle / Math.PI;

		int index = findApplicableTorsionIndex(angle);

		double severity = getSeverity(angle, index);
		if (severity == 0.0)
			return 0.0;

        double dif = Math.abs(torsion - mTorsion[index]);
        if (dif > 180)
            dif = 360 - dif;
        if (dif > 60)
            dif = 60;

	    StereoMolecule mol = conformer.getMolecule();

	    double penalty = severity*dif*dif/14400;	// make a 60 degree dif penalty the same as 0.5 Angstrom distance penalty

		double totalStrain = 0;
	    for (int i=0; i<mol.getAllConnAtoms(mAtom[1]); i++) {
	        if (mol.getConnAtom(mAtom[1], i) != mAtom[2]) {
	        	atomStrain[mol.getConnAtom(mAtom[1], i)] += penalty;
				totalStrain += penalty;
	        	}
	    	}
        for (int i=0; i<mol.getAllConnAtoms(mAtom[2]); i++) {
            if (mol.getConnAtom(mAtom[2], i) != mAtom[1]) {
            	atomStrain[mol.getConnAtom(mAtom[2], i)] += penalty;
				totalStrain += penalty;
		    	}
			}
		return totalStrain;
		}

	public boolean disableIfColliding(SelfOrganizedConformer conformer) {
	    StereoMolecule mol = conformer.getMolecule();
		double maxAtomStrain = 0;
		for (int i=1; i<=2; i++) {
		    for (int j=0; j<mol.getAllConnAtoms(mAtom[i]); j++) {
		    	int connAtom = mol.getConnAtom(mAtom[i], j);
		        if (connAtom != mAtom[3-i] && conformer.getAtomStrain(connAtom) > maxAtomStrain)
		        	maxAtomStrain = conformer.getAtomStrain(connAtom);
		    	}
			}

System.out.println(((maxAtomStrain < COLLIDING_ATOM_STRAIN) ? "kept:" : "disabled:")+toString());

		if (maxAtomStrain < COLLIDING_ATOM_STRAIN)
			return false;

		mIsEnabled = false;
		return true;
		}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder("torsion rule:");
		super.addAtomList(sb);
        sb.append(" torsions:");
        for (int i=0; i<mTorsion.length; i++) {
	        if (i != 0)
		        sb.append(",");
	        sb.append(mTorsion[i] + "(" + mFrequency[i] + ")");
            }
		return sb.toString();
		}
	}
