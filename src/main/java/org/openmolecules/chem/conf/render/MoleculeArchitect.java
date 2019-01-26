package org.openmolecules.chem.conf.render;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.VDWRadii;

/**
 * A MoleculeArchitect is a class that constructs a 3D molecule model
 * from cylinders and spheres. It talks to a MoleculeBuilder class that
 * implements createCylinder() and createSphere() to build the molecule
 * from primitives for a particular environment, e.g. for 3D-JavaFX or
 * the SunFlow ray-tracer.
 */
public class MoleculeArchitect {
	public final static String[] MODE_TEXT = { "Ball & Sticks", "Sticks", "Balls", "Wires" };
	public final static int CONSTRUCTION_MODE_BALL_AND_STICKS = 0;
	public final static int CONSTRUCTION_MODE_STICKS = 1;
	public final static int CONSTRUCTION_MODE_BALLS = 2;
	public final static int CONSTRUCTION_MODE_WIRES = 3;
	public static final int CONSTRUCTION_MODE_DEFAULT = CONSTRUCTION_MODE_BALL_AND_STICKS;

	public final static int HYDROGEN_MODE_ALL = 0;
	public final static int HYDROGEN_MODE_POLAR = 1;
	public final static int HYDROGEN_MODE_NONE = 2;
	public static final int HYDROGEN_MODE_DEFAULT = HYDROGEN_MODE_ALL;

	public final static int COLOR_MODE_ATOMIC_NO = 0;
	public final static int COLOR_MODE_WIRES = 1;
	public static final int COLOR_MODE_DEFAULT = COLOR_MODE_ATOMIC_NO;

	private static final double BALL_AND_STICK_SBOND_RADIUS = 0.18;
	private static final double BALL_AND_STICK_DBOND_RADIUS = 0.10;
	private static final double BALL_AND_STICK_TBOND_RADIUS = 0.08;
	private static final double BALL_AND_STICK_DBOND_SHIFT = 0.15;
	private static final double BALL_AND_STICK_TBOND_SHIFT = 0.26;

	public static final double STICK_SBOND_RADIUS = 0.16;
	private static final double STICK_DBOND_RADIUS = 0.08;
	private static final double STICK_PI_BOND_SHIFT = 0.40;

	// Radii of metal ligand bond dots
	private static final double BALL_AND_STICK_DOT_RADIUS = 0.12;
	private static final double STICK_DOT_RADIUS = 0.10;
	private static final double WIRE_DOT_RADIUS = 0.05;

	public static final double WIRE_SBOND_RADIUS = 0.04;
	private static final double WIRE_DBOND_RADIUS = 0.04;
	private static final double WIRE_PI_BOND_SHIFT = 0.18;

	public static final int BALL_AND_STICK_STICK_COLOR = 0xFFE0E0E0;

	public final static int[] ATOM_ARGB = {
			0xFFFF1493, 0xFFFFFFFF, 0xFFD9FFFF, 0xFFCC80FF, 0xFFC2FF00, 0xFFFFB5B5, //  ?, H,He,Li, Be,B
			0xFF909090, 0xFF3050F8, 0xFFFF0D0D, 0xFF90E050, 0xFFB3E3F5, 0xFFAB5CF2, //  C, N, O, F,Ne,Na
			0xFF8AFF00, 0xFFBFA6A6, 0xFFF0C8A0, 0xFFFF8000, 0xFFFFFF30, 0xFF1FF01F, // Mg,Al,Si, P, S,Cl
			0xFF80D1E3, 0xFF8F40D4, 0xFF3DFF00, 0xFFE6E6E6, 0xFFBFC2C7, 0xFFA6A6AB, // Ar, K,Ca,Sc,Ti, V
			0xFF8A99C7, 0xFF9C7AC7, 0xFFE06633, 0xFFF090A0, 0xFF50D050, 0xFFC88033, // Cr,Mn,Fe,Co,Ni,Cu
			0xFF7D80B0, 0xFFC28F8F, 0xFF668F8F, 0xFFBD80E3, 0xFFFFA100, 0xFFA62929, // Zn,Ga,Ge,As,Se,Br
			0xFF5CB8D1, 0xFF702EB0, 0xFF00FF00, 0xFF94FFFF, 0xFF94E0E0, 0xFF73C2C9, // Kr,Rb,Sr, Y,Zr,Nb
			0xFF54B5B5, 0xFF3B9E9E, 0xFF248F8F, 0xFF0A7D8C, 0xFF006985, 0xFFC0C0C0, // Mo,Tc,Ru,Rh,Pd,Ag
			0xFFFFD98F, 0xFFA67573, 0xFF668080, 0xFF9E63B5, 0xFFD47A00, 0xFF940094, // Cd,In,Sn,Sb,Te, I
			0xFF429EB0, 0xFF57178F, 0xFF00C900, 0xFF70D4FF, 0xFFFFFFC7, 0xFFD9FFC7, // Xe,Cs,Ba,La,Ce,Pr
			0xFFC7FFC7, 0xFFA3FFC7, 0xFF8FFFC7, 0xFF61FFC7, 0xFF45FFC7, 0xFF30FFC7, // Nd,Pm,Sm,Eu,Gd,Tb
			0xFF1FFFC7, 0xFF00FF9C, 0xFF00E675, 0xFF00D452, 0xFF00BF38, 0xFF00AB24, // Dy,Ho,Er,Tm,Yb,Lu
			0xFF4DC2FF, 0xFF4DA6FF, 0xFF2194D6, 0xFF267DAB, 0xFF266696, 0xFF175487, // Hf,Ta, W,Re,Os,Ir
			0xFFD0D0E0, 0xFFFFD123, 0xFFB8B8D0, 0xFFA6544D, 0xFF575961, 0xFF9E4FB5, // Pt,Au,Hg,Tl,Pb,Bi
			0xFFAB5C00, 0xFF754F45, 0xFF428296, 0xFF420066, 0xFF007D00, 0xFF70ABFA, // Po,At,Rn,Fr,Ra,Ac
			0xFF00BAFF, 0xFF00A1FF, 0xFF008FFF, 0xFF0080FF, 0xFF006BFF, 0xFF545CF2, // Th,Pa, U,Np,Pu,Am
			0xFF785CE3, 0xFF8A4FE3, 0xFFA136D4, 0xFFB31FD4, 0xFFB31FBA, 0xFFB30DA6, // Cm,Bk,Cf,Es,Fm,Md
			0xFFBD0D87, 0xFFC70066, 0xFFCC0059, 0xFFD1004F, 0xFFD90045, 0xFFE00038, // No,Lr,Rf,Db,Sg,Bh
			0xFFE6002E, 0xFFEB0026,													// Hs,Mt
	};


	private MoleculeBuilder mBuilder;
	private Coordinates center,delta,point1,point2;    // reused Coordinate objects
	private int mColorMode,mConstructionMode,mHydrogenMode;

	public MoleculeArchitect(MoleculeBuilder builder) {
		mBuilder = builder;
		center = new Coordinates();
		delta = new Coordinates();
		point1 = new Coordinates();
		point2 = new Coordinates();
		mConstructionMode = CONSTRUCTION_MODE_DEFAULT;
		mHydrogenMode = HYDROGEN_MODE_ALL;
		mColorMode = COLOR_MODE_DEFAULT;
		}

	public int getColorMode() {
		return mColorMode;
	}

	public int getConstructionMode() {
		return mConstructionMode;
	}

	public int getHydrogenMode() {
		return mHydrogenMode;
	}

	public void setColorMode(int mode) {
		mColorMode = mode;
	}

	public void setConstructionMode(int mode) {
		mConstructionMode = mode;
	}

	public void setHydrogenMode(int mode) {
		mHydrogenMode = mode;
	}

	public MoleculeBuilder getRenderer() {
		return mBuilder;
		}

/*	public void centerMolecule(Conformer conformer) {
		Coordinates cog = new Coordinates();
		double atomicNoSum = 0.0;
		StereoMolecule mol = conformer.getMolecule();
		for (int atom=0; atom<mol.getAllAtoms(); atom++) {
			Coordinates c = conformer.getCoordinates(atom);
			cog.add(c.x * mol.getAtomicNo(atom),
					c.y * mol.getAtomicNo(atom),
					c.z * mol.getAtomicNo(atom));
			atomicNoSum += mol.getAtomicNo(atom);
			}
		cog.scale(1.0 / atomicNoSum);

		// Here we move the conformer's internal coordinates to center it around its COG!!!!!!!
		for (int atom=0; atom<mol.getAllAtoms(); atom++)
			conformer.getCoordinates(atom).sub(cog);
		}*/

	public void buildMolecule(Conformer conformer) {
		StereoMolecule mol = conformer.getMolecule();
		mol.ensureHelperArrays(Molecule.cHelperRings);

		mBuilder.init();

		if (mConstructionMode == CONSTRUCTION_MODE_STICKS
		 || mConstructionMode == CONSTRUCTION_MODE_BALL_AND_STICKS
		 || mConstructionMode == CONSTRUCTION_MODE_WIRES)
			for (int bond=0; bond<mol.getAllBonds(); bond++)
				if (!skipAtom(mol, mol.getBondAtom(0, bond))
				 && !skipAtom(mol, mol.getBondAtom(1, bond)))
					buildBond(conformer, bond);

		if (mConstructionMode == CONSTRUCTION_MODE_STICKS
		 || mConstructionMode == CONSTRUCTION_MODE_BALL_AND_STICKS
		 || mConstructionMode == CONSTRUCTION_MODE_BALLS) {
			for (int atom=0; atom<mol.getAllAtoms(); atom++) {
				if (!skipAtom(mol, atom)) {
					int atomicNo = mol.getAtomicNo(atom);
					double radius = mol.isMarkedAtom(atom) ? VDWRadii.VDW_RADIUS[atomicNo]/4
								  : (mConstructionMode == CONSTRUCTION_MODE_STICKS) ? STICK_SBOND_RADIUS
								  : (mConstructionMode == CONSTRUCTION_MODE_BALLS) ? VDWRadii.VDW_RADIUS[atomicNo]*0.95  // to avoid collision with vdw-radii based surface
								  :							VDWRadii.VDW_RADIUS[atomicNo]/4;
					mBuilder.addSphere(atom, -1, conformer.getCoordinates(atom), radius, getAtomColor(mol, atom));
					}
				}
			}

		mBuilder.done();
		}

	private int getAtomColor(StereoMolecule mol, int atom) {
		return ATOM_ARGB[mol.getAtomicNo(atom)];
		}

	private boolean skipAtom(StereoMolecule mol, int atom) {
		return (mHydrogenMode == HYDROGEN_MODE_ALL) ? false
			 : (mHydrogenMode == HYDROGEN_MODE_NONE) ? true
			 : (mol.getConnAtoms(atom) != 1) ? false
			 : (mol.getAtomicNo(mol.getConnAtom(atom, 0)) == 6);
		}

	private void buildBond(Conformer conformer, int bond) {
		StereoMolecule mol = conformer.getMolecule();
		int atom1 = mol.getBondAtom(0, bond);
		int atom2 = mol.getBondAtom(1, bond);
		Coordinates c1 = conformer.getCoordinates(atom1);
		Coordinates c2 = conformer.getCoordinates(atom2);
		center.center(c1, c2);
		delta.set(c2).sub(c1);

		double d = delta.getLength();
		double dxy = Math.sqrt(delta.x * delta.x + delta.y * delta.y);
		double b = Math.asin(c2.z > c1.z ? dxy / d : -dxy / d);
		double c = (delta.x < 0.0) ? Math.atan(delta.y / delta.x) + Math.PI
				: (delta.x > 0.0) ? Math.atan(delta.y / delta.x)
				: (delta.y > 0.0) ? Math.PI / 2 : -Math.PI / 2;

		switch (mConstructionMode) {
		case CONSTRUCTION_MODE_BALL_AND_STICKS:
			buildBallAndStickBond(conformer, bond, d, b, c);
			break;
		case CONSTRUCTION_MODE_STICKS:
			buildStickBond(conformer, bond, d, b, c);
			break;
		case CONSTRUCTION_MODE_WIRES:
			buildStickBond(conformer, bond, d, b, c);
			break;
			}
		}

	private void buildBallAndStickBond(Conformer conformer, int bond, double d, double b, double c) {
		StereoMolecule mol = conformer.getMolecule();
		int color = BALL_AND_STICK_STICK_COLOR;
		int order = mol.getBondOrder(bond);
		if (order == 1) {
			double dd = 2*calculateBondReduction(mol, bond, 0.2);
			if (dd < d)
				mBuilder.addCylinder(bond, BALL_AND_STICK_SBOND_RADIUS, d-dd, center, b, c, color);
			return;
			}

		if (order == 2) {
			Coordinates ds = calculateDoubleBondShift(conformer, bond).scale(BALL_AND_STICK_DBOND_SHIFT);
			double dd = calculateBondReduction(mol, bond, 0.20+0.10);
			if (dd != 0f) {
				mBuilder.addCylinder(bond, BALL_AND_STICK_DBOND_RADIUS, d-dd, point1.set(center).add(ds), b, c, color);
				mBuilder.addCylinder(bond, BALL_AND_STICK_DBOND_RADIUS, d-dd, point1.set(center).sub(ds), b, c, color);
				}
			return;
			}

		if (order == 3) {
			Coordinates ds = calculateRandomOrthogonalShift(conformer, bond).scale(BALL_AND_STICK_TBOND_SHIFT);
			double dd1 = 2*calculateBondReduction(mol, bond, 0.11);
			double dd2 = 2*calculateBondReduction(mol, bond, 0.22+0.07);
			if (dd2 < d)
				mBuilder.addCylinder(bond, BALL_AND_STICK_TBOND_RADIUS, d-dd2, point1.set(center).add(ds), b, c, color);
			if (dd1 < d)
				mBuilder.addCylinder(bond, BALL_AND_STICK_TBOND_RADIUS, d-dd1, center, b, c, color);
			if (dd2 < d)
				mBuilder.addCylinder(bond, BALL_AND_STICK_TBOND_RADIUS, d-dd2, point1.set(center).sub(ds), b, c, color);
			return;
			}

		if (order == 0) {
			double dd = calculateBondReduction(mol, bond, 0.2);
			if (2*dd < d) {
				Coordinates p1 = conformer.getCoordinates(mol.getBondAtom(0, bond));
				Coordinates p2 = conformer.getCoordinates(mol.getBondAtom(1, bond));
				point1.between(p1, p2, dd/d);
				point2.between(p2, p1, dd/d);
				buildDottedBond(bond, color, color, point1, point2, BALL_AND_STICK_DOT_RADIUS, d-2*dd);
				}
			return;
			}
		}

	private void buildStickBond(Conformer conformer, int bond, double d, double b, double c) {
		StereoMolecule mol = conformer.getMolecule();

		int color1 = getAtomColor(mol, mol.getBondAtom(0, bond));
		int color2 = getAtomColor(mol, mol.getBondAtom(1, bond));

		int atom1 = mol.getBondAtom(0, bond);
		int atom2 = mol.getBondAtom(1, bond);

		Coordinates p1 = conformer.getCoordinates(atom1);
		Coordinates p2 = conformer.getCoordinates(atom2);

		int order = mol.getBondOrder(bond);

		if (order == 0) {
			double r = (mConstructionMode == CONSTRUCTION_MODE_WIRES) ? WIRE_DOT_RADIUS : STICK_DOT_RADIUS;
			buildDottedBond(bond, color1, color2, p1, p2, r, d);
			return;
			}

		double r1 = (mConstructionMode == CONSTRUCTION_MODE_WIRES) ? WIRE_SBOND_RADIUS : STICK_SBOND_RADIUS;

		if (order == 1) {
			buildStickBond(bond, color1, color2, p1, p2, r1, d, b, c);
			return;
			}

		double r2 = (mConstructionMode == CONSTRUCTION_MODE_WIRES) ? WIRE_DBOND_RADIUS : STICK_DBOND_RADIUS;
		double piShift = (mConstructionMode == CONSTRUCTION_MODE_WIRES) ? WIRE_PI_BOND_SHIFT : STICK_PI_BOND_SHIFT;

		if (order == 2) {
			Coordinates ds = calculateDoubleBondShift(conformer, bond).scale(piShift);
			buildStickBond(bond, color1, color2, p1, p2, r1, d, b, c);
			buildPiStickBond(mol, bond, color1, color2, point1.set(p1).add(ds),
														point2.set(p2).add(ds), r2, piShift, d, b, c);
			return;
			}

		if (order == 3) {
			Coordinates ds = calculateRandomOrthogonalShift(conformer, bond).scale(piShift);
			buildStickBond(bond, color1, color2, p1, p2, r1, d, b, c);
			buildPiStickBond(mol, bond, color1, color2, point1.set(p1).sub(ds),
														point2.set(p2).sub(ds), r2, piShift, d, b, c);
			buildPiStickBond(mol, bond, color1, color2, point1.set(p1).sub(ds),
														point2.set(p2).sub(ds), r2, piShift, d, b, c);
			return;
			}
		}

	private void buildStickBond(int bond, int color1, int color2, Coordinates p1, Coordinates p2,
	                            double r, double d, double b, double c) {
		if (color1 == color2) {
			center.center(p1, p2);
			mBuilder.addCylinder(bond, r, d, center, b, c, color1);
			}
		else {
			center.between(p1, p2, 0.25);
			mBuilder.addCylinder(bond, r, d / 2, center, b, c, color1);
			center.between(p1, p2, 0.75);
			mBuilder.addCylinder(bond, r, d / 2, center, b, c, color2);
			}
		}

	private void buildDottedBond(int bond, int color1, int color2, Coordinates p1, Coordinates p2,
								 double r, double d) {
		int dots = 2 * Math.max(2, (int)Math.round(d / (r * 5)));
		double dd = r / (dots-1);
		for (int i=1; i<dots-1; i++) {
			center.between(p1, p2, (double)i/(dots-1));
			mBuilder.addSphere(-1, bond, center, r, i*2<dots ? color1 : color2);
			}
		}

	/**
	 * Builds the second (and third in case of triple bond) shortened and thinner stick
	 * that indicates the unsaturation in STICK mode. Both stick ends are rounded by adding
	 * a sphere and the stick is created from two differently colored cylinders, if both
	 * and atoms have different atomic numbers.
	 * @param mol
	 * @param bond
	 * @param color1 color of bond atom 0
	 * @param color2 color of bond atom 1
	 * @param c1 bond atom 0 position parallel shifted
	 * @param c2 bond atom 1 position parallel shifted
	 * @param r stick radius
	 * @param piShift distance between sigma bond and thinner pi bond
	 * @param d distance between both bond atoms
	 * @param b rotation around y axis
	 * @param c rotation around z axis
	 */
	private void buildPiStickBond(StereoMolecule mol, int bond, int color1, int color2,
	                              Coordinates c1, Coordinates c2, double r, double piShift, double d, double b, double c) {
		center.center(c1, c2);

		double centerShift = piShift / (2 * d);

		double l1 = d / 2;
		Coordinates p1 = c1;
		if (mol.getAllConnAtoms(mol.getBondAtom(0, bond)) != 1) {
			p1 = point1.between(c1, c2, centerShift);
			l1 -= piShift / 2;
			}

		double l2 = d / 2;
		Coordinates p2 = c2;
		if (mol.getAllConnAtoms(mol.getBondAtom(1, bond)) != 1) {
			p2 = point2.between(c1, c2, 1.0 - centerShift);
			l2 -= piShift / 2;
			}

		if (mConstructionMode != CONSTRUCTION_MODE_WIRES) {
			mBuilder.addSphere(-1, -1, p1, r, color1);
			mBuilder.addSphere(-1, -1, p2, r, color2);
			}

		if (color1 == color2) {
			center.center(p1, p2);
			mBuilder.addCylinder(bond, r, l1 + l2, center, b, c, color1);
			}
		else {
			p1.center(center);
			mBuilder.addCylinder(bond, r, l1, p1, b, c, color1);
			p2.center(center);
			mBuilder.addCylinder(bond, r, l2, p2, b, c, color2);
			}
		}

	private Coordinates calculateDoubleBondShift(Conformer conformer, int bond) {
		StereoMolecule mol = conformer.getMolecule();
		int maxScore = 0;
		int maxNeighbor = -1;
		int maxSide = -1;
		for (int i=0; i<2; i++) {
			int atom = mol.getBondAtom(i, bond);
			int otherAtom = mol.getBondAtom(1-i, bond);
			for (int j=0; j<mol.getConnAtoms(atom); j++) {
				int connAtom = mol.getConnAtom(atom, j);
				if (connAtom != otherAtom) {
					int score = getReferenceNeighborScore(mol, atom, bond, connAtom, mol.getConnBond(atom, j));
					if (maxScore < score) {
						maxScore = score;
						maxNeighbor = connAtom;
						maxSide = i;
						}
					}
				}
			}

		if (maxNeighbor != -1) {
			Coordinates ca = conformer.getCoordinates(mol.getBondAtom(maxSide, bond));
			Coordinates vb = conformer.getCoordinates(mol.getBondAtom(1-maxSide, bond)).subC(ca).unit();
			Coordinates vc = conformer.getCoordinates(maxNeighbor).subC(ca).unit();
			vb.scale(-vb.dot(vc));
			return vc.addC(vb).unit();
			}

		return calculateRandomOrthogonalShift(conformer, bond);
		}

	/**
	 * Returns the preferred double bond neighbor. Double bonds are constructed as one
	 * stick connecting the involved atoms and another parallel one moved to that side,
	 * where the reference bond neighbor points to.
	 * @param mol
	 * @param atom
	 * @param bond
	 * @param connAtom
	 * @param connBond
	 * @return
	 */
	private int getReferenceNeighborScore(StereoMolecule mol, int atom, int bond, int connAtom, int connBond) {
		if (mol.isRingBond(bond)) {
			RingCollection ringSet = mol.getRingSet();
			int score = 512;
			for (int i=0; i<ringSet.getSize(); i++) {
				if (ringSet.isBondMember(i, bond)) {
					if (ringSet.isBondMember(i, connBond)) {
						if (ringSet.isAromatic(i))
							score += 256;
						score += 256 - ringSet.getSize();
						}
					}
				}
			return score;
			}
		else {
			return (mol.getConnAtoms(atom) == 2 ? 256 : 0) + mol.getAtomicNo(connAtom);
			}
		}

	private Coordinates calculateRandomOrthogonalShift(Conformer conformer, int bond) {
		StereoMolecule mol = conformer.getMolecule();
		Coordinates c1 = conformer.getCoordinates(mol.getBondAtom(0, bond));
		Coordinates c2 = conformer.getCoordinates(mol.getBondAtom(1, bond));
		Coordinates v = c2.subC(c1).unit();
		double lxy = Math.sqrt(v.x*v.x+v.y*v.y);
		return (lxy == 0.0) ? new Coordinates(1.0, 0.0, 0.0)
							: new Coordinates(-v.x*v.z/lxy, -v.y*v.z/lxy, lxy);
	}

	private double calculateBondReduction(StereoMolecule mol, int bond, double sideShift) {
		double atomRadius = Math.min(VDWRadii.VDW_RADIUS[mol.getAtomicNo(mol.getBondAtom(0, bond))],
									 VDWRadii.VDW_RADIUS[mol.getAtomicNo(mol.getBondAtom(1, bond))]) / 4;
		return (sideShift >= atomRadius) ? 0.0 : Math.sqrt(atomRadius*atomRadius - sideShift*sideShift);
		}
	}
