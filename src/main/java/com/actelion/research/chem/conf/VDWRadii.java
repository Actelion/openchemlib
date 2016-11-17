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
		0.0f,	1.20f,	1.43f,	2.12f,	// ??,H,He,Li
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
	 * Covalent Radii taken from the JMol program 
	 * Originally taken from http://www.daylight.com/meetings/mug01/Sayle/m4xbondage.html, 
	 * I found them to be inacurrate and incomplete
	 */
	public static final float[] COVALENT_RADIUS = new float[]{
		0.00f,	0.23f,	0.93f,	0.68f,	//0
		0.35f,	0.83f,	0.77f,	0.75f,	//4
		0.73f,	0.64f,	1.12f,	0.97f,	//8
		1.10f,	1.35f,	1.12f,	0.75f,	//12
		1.02f,	0.99f,	1.57f,	1.33f,	//16
		0.99f,	1.44f,	1.47f,	1.33f,	//20
		1.35f,	1.35f,	1.34f,	1.33f,	//24
		1.50f,	1.52f,	1.45f,	1.22f,	//28
		1.17f,	1.12f,	1.22f,	1.21f,	//32
		1.91f,	1.47f,	1.12f,	1.78f,	//36
		1.56f,	1.48f,	1.47f,	1.35f,	//40
		1.40f,	1.45f,	1.50f,	1.59f,	//44
		1.69f,	1.63f,	1.46f,	1.46f,	//48
		1.47f,	1.40f,	1.98f,	1.67f,	//52
		1.34f,	1.87f,	1.83f,	1.82f	//56		
	};
}
