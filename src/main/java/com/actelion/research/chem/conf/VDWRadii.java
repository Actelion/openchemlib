package com.actelion.research.chem.conf;

public interface VDWRadii {
	
	/** 
	 * VDW Radii indexed by atomic numbers, taken from DOI: 10.1039/c3dt50599e
	 * Santiago Alvarez: A cartography of the van der Waals territories
	 * Published on 01 May 2013 on http://pubs.rsc.org | doi:10.1039/C3DT50599E
	 * 
	 * (in a few cases no value was given (marked ??), where we use estimated values)
	 */
	public static final float[] VDW_RADIUS = new float[]{
		1.00f,	1.20f,	1.43f,	2.12f,	// ??,H,He,Li
		1.98f,	1.91f,	1.77f,	1.66f,	// Be,B,C,N
		1.50f,	1.46f,	1.58f,	2.50f,	// O,F,Ne,Na
		2.51f,	2.25f,	2.19f,	1.90f,	// Mg,Al,Si,P
		1.89f,	1.82f,	1.83f,	2.73f,	// S,Cl,Ar,K
		2.62f,	2.58f,	2.46f,	2.42f,	// Ca,Sc,Ti,V
		2.45f,	2.45f,	2.44f,	2.40f,	// Cr,Mn,Fe,Co
		2.40f,	2.38f,	2.39f,	2.32f,	// Ni,Cu,Zn,Ga
		2.29f,	1.88f,	1.82f,	1.86f,	// Ge,As,Se,Br
		2.25f,	3.21f,	2.84f,	2.75f,	// Kr,Rb,Sr,Y
		2.52f,	2.56f,	2.45f,	2.44f,	// Zr,Nb,Mo,Tc
		2.46f,	2.44f,	2.15f,	2.53f,	// Ru,Rh,Pd,Ag
		2.49f,	2.43f,	2.42f,	2.47f,	// Cd,In,Sn,Sb
		1.99f,	2.04f,	2.06f,	3.48f,	// Te,I,Xe,Cs
		3.03f,	2.98f,	2.88f,	2.92f,	// Ba,La,Ce,Pr
		2.95f,	2.92f,	2.90f,	2.87f,	// Nd,??,Sm,Eu
		2.83f,	2.79f,	2.87f,	2.81f,	// Gd,Tb,Dy,Ho
		2.83f,	2.79f,	2.80f,	2.74f,	// Er,Tm,Yb,Lu
		2.63f,	2.53f,	2.57f,	2.49f,	// Hf,Ta,W,Re
		2.48f,	2.41f,	2.29f,	2.32f,	// Os,Ir,Pt,Au
		2.45f,	2.47f,	2.60f,	2.54f, 	// Hg,Tl,Pb,Bi
		2.5f,	2.5f,	2.5f,	2.5f,	// ??,??,??,??
		2.5f,	2.80f,	2.93f,	2.88f,	// ??,Ac,Th,Pa
		2.71f,	2.82f,	2.81f,	2.83f,	// U,Np,Pu,Am
		3.05f,	3.40f,	3.05f,	2.70f	// Cm,Bk,Cf,Es
	};

	/**
	 * These covalent single bond radii values were taken from the
	 * following source:
	 * 'Molecular Double-Bond Covalent Radii for Elements Li–E112', 2009,
	 * Pyykkö and Atsumi, doi: 10.1002/chem.200901472
	 */
	public static final float[] COVALENT_RADIUS = new float[] {
			0.25f, 0.32f, 0.46f, 1.33f, // ?,H,He,Li
			1.02f, 0.85f, 0.75f, 0.71f, // Be,B,C,N
			0.63f, 0.64f, 0.96f, 1.60f, // O,F,Ne,Na
			1.39f, 1.26f, 1.16f, 1.11f, // Mg,Al,Si,P
			1.03f, 0.99f, 1.07f, 1.96f, // S,Cl,Ar,K
			1.71f, 1.48f, 1.36f, 1.34f, // Ca,Sc,Ti,V
			1.22f, 1.19f, 1.16f, 1.11f, // Cr,Mn,Fe,Co
			1.10f, 1.20f, 1.20f, 1.24f, // Ni,Cu,Zn,Ga
			1.21f, 1.21f, 1.16f, 1.14f, // Ge,As,Se,Br
			1.21f, 2.10f, 1.85f, 1.63f, // Kr,Rb,Sr,Y
			1.54f, 1.47f, 1.38f, 1.28f, // Zr,Nb,Mo,Tc
			1.25f, 1.25f, 1.20f, 1.39f, // Ru,Rh,Pd,Ag
			1.44f, 1.46f, 1.40f, 1.40f, // Cd,In,Sn,Sb
			1.36f, 1.33f, 1.35f, 2.32f, // Te,I,Xe,Cs
			1.96f, 1.80f, 1.63f, 1.76f, // Ba,La,Ce,Pr
			1.74f, 1.73f, 1.72f, 1.68f, // Nd,Pm,Sm,Eu
			1.69f, 1.68f, 1.67f, 1.66f, // Gd,Tb,Dy,Ho
			1.65f, 1.64f, 1.70f, 1.62f, // Er,Tm,Yb,Lu
			1.52f, 1.46f, 1.37f, 1.31f, // Hf,Ta,W,Re
			1.29f, 1.22f, 1.23f, 1.24f, // Os,Ir,Pt,Au
			1.42f, 1.50f, 1.44f, 1.51f, // Hg,Tl,Pb,Bi
			1.45f, 1.47f, 1.45f, 2.23f, // Po,At,Rn,Fr
			2.01f, 1.86f, 1.75f, 1.69f, // Ra,Ac,Th,Pa
			1.70f, 1.71f, 1.72f, 1.66f, // U,Np,Pu,Am
			1.66f, 1.68f, 1.68f, 1.65f, // Cm,Bk,Cf,Es
			1.67f, 1.73f, 1.76f, 1.61f, // Fm,Md,No,Lr
			1.57f, 1.49f, 1.43f, 1.41f, // Rf,Db,Sg,Bh
			1.34f, 1.29f, 1.28f, 1.21f, // Hs,Mt,Ds,Rg
			1.37f, 1.36f, 1.43f, 1.62f, // Cn,Uut,Fl,Uup
			1.75f, 1.65f, 1.57f         // Lv, Uus, Uuo
	};

	public static float getVDWRadius(int atomicNo) {
		return VDW_RADIUS[atomicNo < VDW_RADIUS.length ? atomicNo : 6]; // we assume some kind of carbon from the MDL special types
	}

	public static float getCovalentRadius(int atomicNo) {
		return COVALENT_RADIUS[atomicNo < COVALENT_RADIUS.length ? atomicNo : 6]; // we assume some kind of carbon from the MDL special types
	}
}
