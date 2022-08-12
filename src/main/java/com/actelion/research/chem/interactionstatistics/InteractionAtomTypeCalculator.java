package com.actelion.research.chem.interactionstatistics;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.PeriodicTable;
import com.actelion.research.chem.StereoMolecule;


public class InteractionAtomTypeCalculator {
	
	
	
	public enum FunctionalGroup {NITRO("NO2",1),ESTER("COOR",2),CARBOXYL("COO",3),
			SULFONAMIDE("HNSO2R",4),SULFONATE("SO3",5), PHOSPHONATE("PO3",6),SULFOXIDE("SO",7),
			SULFONE("SO2",8),AMIDE("HNCO",9),AMIDINE("N=C-N",10),GUANIDINE("N-C(-N)=N",11),N_SP2_TAUT("N=C-N(Ar)",12),
			SP2_AMINE("C=CNX2",13),ENOL("C=COH",14);
	
			private String s;
			private int id;
			
			FunctionalGroup(String s, int id){
				this.s = s;
				this.id = id;
			}
			
			public String getString() {return s;}
			
			public int getId() {return id;}

		@Override
		public String toString() {
			final StringBuilder sb = new StringBuilder("FunctionalGroup{");
			sb.append("s='").append(s).append('\'');
			sb.append(", id=").append(id);
			sb.append('}');
			return sb.toString();
		}
	}
	

	
	public enum AtomPropertyMask { ATOMIC_NO(0x0000007F),HYBRID(0x00000180),NONHNEIGHBOURS(0x00000E00),AROM(0x00001000),
		STABILIZED(0x00002000), FUNCTIONAL_GROUP(0x000FC000), BASIC(0x000001FF), EXTENDED(0x00003FFF),SPECIFIC(0x000FCFFF);

		private final int mask;
		AtomPropertyMask(int mask) { this.mask = mask;}
		public int getMask() { return mask;}
	}
	


	
	public enum AtomPropertyShift{HYBRID_SHIFT(7), NEIGHBOURS_SHIFT(9), FUNCTIONAL_GROUP_SHIFT(16);
		private final int shift;
		AtomPropertyShift(int shift) { this.shift = shift;}
		public int getShift() { return shift;}
	}

	
	public enum AtomFlagCount{FUNC_GROUP_FLAG_COUNT(20),BASIC_ATOM_FLAG_COUNT(12), EXTENDED_ATOM_FLAG_COUNT(14);
		private final int count;
		AtomFlagCount(int count) { this.count = count;}
		public int getCount() { return count;}
	}
	
	
	private static int getFunctionalGroupValue(StereoMolecule complex, int atom) {
		return getFunctionalGroup(complex,atom)<<AtomPropertyShift.FUNCTIONAL_GROUP_SHIFT.getShift();
		
	}
	
	private static int getFunctionalGroup(StereoMolecule complex,int atom) {
		complex.ensureHelperArrays(Molecule.cHelperNeighbours);
		int atomNo = complex.getAtomicNo(atom);
		int neighbours = 0;
		int doubleBonds = 0;
		int tripleBonds = 0;
		

		for(int i=0; i<complex.getConnAtoms(atom); i++) {
			int order = complex.getConnBondOrder(atom, i);
			if(order==2) doubleBonds++;
			else if(order==3) tripleBonds++;
			neighbours++;
		}

		if(atomNo==8) { //oxygen
			if(neighbours==1) {
				if(getPaths(complex,atom,7,-1,8,-1)>0)return FunctionalGroup.NITRO.id;
				//int nPOs = getPaths(complex,a,15,-1,8,-1);
				int nPs = getNeighbours(complex,atom,15,-1);
				if(nPs>0) {
					int p = complex.getConnAtom(atom, 0);
					int nPOXs = getNeighbours(complex,p,8,-1);
					int nPORs = getPaths(complex,p,8,-1,-1,-1);
					int nPOs = nPOXs-nPORs;
					if(nPOs==1) return 0;
					else if(nPOs>=2) return FunctionalGroup.PHOSPHONATE.id;
				}
				int nSs = getNeighbours(complex,atom,16,-1);
				if(nSs>0) {
					int s = complex.getConnAtom(atom, 0);
					int nSNs = getPaths(complex,atom,16,-1,7,1);
					int nSOXs = getNeighbours(complex,s,8,-1);
					int nSORs = getPaths(complex,s,8,-1,-1,-1);
					int nSOs = nSOXs-nSORs;
					if (nSOs>0 && nSNs>0) return FunctionalGroup.SULFONAMIDE.id;
					else if (nSOs==1) return FunctionalGroup.SULFOXIDE.id;
					else if (nSOs==2) return FunctionalGroup.SULFONE.id;
					else if (nSOs==3) return FunctionalGroup.SULFONATE.id;
				}
			}
			int nCOs = getPaths(complex,atom,6,-1,8,-1);
				
			if(nCOs>0) {
				for(int i=0;i<complex.getConnAtoms(atom);i++) {
					int aa = complex.getConnAtom(atom,i);
					if(complex.getAtomicNo(aa)!=6)continue;
					for(int j=0;j<complex.getConnAtoms(aa);j++) {
						int aaa = complex.getConnAtom(aa, j);
						if(aaa==atom)continue;
						else if (complex.getAtomicNo(aaa)==8) { // C(=O)O
							int dBds = complex.getBondOrder(complex.getBond(aa,aaa))> 1 ? 1:0;
							if(doubleBonds+dBds>0) {
								if(complex.getNonHydrogenNeighbourCount(aaa)==1 && neighbours==1) return FunctionalGroup.CARBOXYL.id;
								else if (complex.getAllConnAtoms(aaa)==1 && neighbours>1) return FunctionalGroup.ESTER.id;
								else if (complex.getAllConnAtoms(aaa)>1 && neighbours==1) return FunctionalGroup.ESTER.id;
							}
							
						}
				
			}
				}
			}
			
			if(doubleBonds<1 && complex.getNonHydrogenNeighbourCount(atom)==1) {
				if(getPaths(complex, atom, 6, 1, 6, 2)>0) {
					return FunctionalGroup.ENOL.id;
			}
			}
				
					
			

				
			if(doubleBonds==1) { //sp2 oxygen
				if(getPaths(complex,atom,6,-1,7,1)>0) { //amide
					return FunctionalGroup.AMIDE.id;
				}

							
			}
			return 0;
			}
		else if(atomNo==7) { //nitrogen
			if(neighbours==4) return 0;
			if(complex.isAromaticAtom(atom)) { //aromatic nitrogen
				/**
				//check for imidazole 
				if(neighbours==2) { //only secondary nitrogens 
					
					imidazole:for(int i=0; i<complex.getAllConnAtoms(a); i++) {
						int aa = complex.getConnAtom(a, i);
						if(!complex.isRingBond(complex.getBond(a, aa)) || complex.getAtomicNo(aa)!=6 || !complex.isAromaticAtom(aa))continue imidazole;
						for(int j=0; j<complex.getAllConnAtoms(a); j++) {
							int aaa = complex.getConnAtom(aa, j);
							if(aaa==a) continue;
							if(complex.getAtomicNo(aaa)==7 && complex.isAromaticAtom(aaa) && complex.getConnAtoms(aaa)==2) return N_IMIDAZOLE; 
						}
					}
				}
				*/
				if(neighbours<=2) { //check for generic tautomeric  nitrogens in aromatic rings
					for(int i=0; i<complex.getAllConnAtoms(atom); i++) {
						int aa = complex.getConnAtom(atom, i);
						if(!complex.isAromaticAtom(aa))continue;
						for(int j=0; j<complex.getConnAtoms(aa); j++) {
							int aaa = complex.getConnAtom(aa, j);
							if(!complex.isAromaticAtom(aaa))continue;
							if(aaa==atom) continue;
							if(complex.getAtomicNo(aaa)==7 && complex.getConnAtoms(aaa)<=2) return FunctionalGroup.N_SP2_TAUT.id;
					}
				}
				
				}

				return 0;
			}
			else { //nitrogen is not aromatic
				if(getNeighbours(complex,atom,8,-1)>=2) return FunctionalGroup.NITRO.id;
				if(neighbours==1 && tripleBonds>0) return 0;
				if(getPaths(complex,atom,6,1,8,2)>0) { //amide
					return FunctionalGroup.AMIDE.id;
				}
				if(getPaths(complex,atom,16,1,8,2)>0) {
					return FunctionalGroup.SULFONAMIDE.id;
				}
			    
				boolean isAmidine = false; 
				
				guanidine:for(int i=0; i<complex.getAllConnAtoms(atom); i++) {
					int aa = complex.getConnAtom(atom, i);
					if(complex.getAtomicNo(aa)==6 && complex.getConnAtoms(aa)==3 && getNeighbours(complex, aa, 7, 2)>0) {
						isAmidine = true;
						for(int j=0; j<complex.getAllConnAtoms(aa); j++) {
							int aaa = complex.getConnAtom(aa, j);
							if(atom==aaa)continue;
							if(complex.getAtomicNo(aaa)!=7 || complex.isAromaticAtom(aaa)) continue guanidine;
							//else isAmidine = true;
						}
						return FunctionalGroup.GUANIDINE.id;

					}
				}
				if(isAmidine) {
					return FunctionalGroup.AMIDINE.id;
				}
				if(doubleBonds==0 && complex.isFlatNitrogen(atom)) {
					return FunctionalGroup.SP2_AMINE.id;
				}
		}
		}
		return 0;
	}
	
	private static int getHybridizationValue(StereoMolecule mol, int atom) {
		if(mol.getAtomicNo(atom) == 6 || mol.getAtomicNo(atom) == 7 || mol.getAtomicNo(atom) == 8) {
		int pi = mol.getAtomPi(atom);
		int sp = (pi == 2) ? 1 : (pi == 1 || mol.getAtomicNo(atom) == 5)  ? 2 : 3;
		
		return sp << AtomPropertyShift.HYBRID_SHIFT.getShift();
		}
		else return 0;

	}
	
	
	
	private static int getAtomicNoValue(StereoMolecule mol, int atom) {
		int atomicNo = mol.getAtomicNo(atom);
		return atomicNo;
		
	}
	

	
	private static int getAromValue(StereoMolecule mol, int atom) {
		int arom = mol.isAromaticAtom(atom) ? 1 : 0;
		int type = 0;
		if(arom>0)
			type = AtomPropertyMask.AROM.getMask();
		
		return type;
	}
	

	private static int getStabilizationValue(StereoMolecule mol, int atom) {
		int stab = mol.isStabilizedAtom(atom) ? 1 : 0;
		int type = 0;
		if(stab>0) 
			type = AtomPropertyMask.STABILIZED.getMask();
		
		return type;
	}
	

	  
	/**
	* get connected atoms of specified element that are bonded to the atom of interest with the specific bond order
	* -1 for any bond order
	* @param mol
	* @param a
	* @param aaAtomicNo
	* @param aaBondOrder
	* @return
	*/
	private static int getNeighbours(StereoMolecule mol,int a, int aaAtomicNo, int aaBondOrder) { 
		int n = 0;
		for(int i=0;i<mol.getAllConnAtoms(a);i++) {
			int aa = mol.getConnAtom(a, i);
			int bondOrder = mol.getBondOrder(mol.getBond(a, aa));
			if(aaAtomicNo>=0 && mol.getAtomicNo(aa)!=aaAtomicNo) continue;
			if(aaBondOrder>0 && bondOrder!=aaBondOrder) continue;
			n++;
		}
		return n;
		
		}



	private static int getPaths(StereoMolecule mol,int a, int aaAtomicNo, int aaBondOrder,int aaaAtomicNo, int aaaBondOrder) { 
		int n = 0;
		for(int i=0;i<mol.getAllConnAtoms(a);i++) {
			int aa = mol.getConnAtom(a, i);
			int bondOrder1 = mol.getBondOrder(mol.getBond(a, aa));
			if(aaAtomicNo>0 &&mol.getAtomicNo(aa)!=aaAtomicNo) continue;
			if(aaBondOrder>0 && bondOrder1!=aaBondOrder) continue;
			for(int j=0;j<mol.getAllConnAtoms(aa);j++) {
				int aaa = mol.getConnAtom(aa, j);
				if(aaa==a)continue;
				int bondOrder2 = mol.getBondOrder(mol.getBond(aa, aaa));
				if(aaaAtomicNo>0 && mol.getAtomicNo(aaa)!=aaaAtomicNo) continue;
				if(aaaBondOrder>0 && bondOrder2!=aaaBondOrder) continue;
				n++;
				
			}
		}
		return n;
		
		}
	
	
	
	public static int getAtomType(StereoMolecule mol, int atom) {
		int type = 0;
		type += getFunctionalGroupValue(mol,atom);
		type += getHybridizationValue(mol,atom);
		type += getStabilizationValue(mol,atom);
		type += getAromValue(mol,atom);
		type += getAtomicNoValue(mol,atom);

		return type;	
	}
	
	public static int getAtomType(FunctionalGroup fg, int atomicNo, boolean isAromatic, int hybridization,
			boolean isStabilized) {
		int type = 0;
		if(fg!=null) {
			type+=fg.id<<AtomPropertyShift.FUNCTIONAL_GROUP_SHIFT.getShift();
		}
		type+=hybridization << AtomPropertyShift.HYBRID_SHIFT.getShift();
		if(isStabilized)
				type += AtomPropertyMask.STABILIZED.getMask();
		if(isAromatic)
			type += AtomPropertyMask.AROM.getMask();
		
		type+= atomicNo;
		
		return type;
	}
	
	public static String getString(int atomType) {
		StringBuilder sb = new StringBuilder();
		sb.append(PeriodicTable.symbol(atomType & AtomPropertyMask.ATOMIC_NO.getMask()));
		int hybrid = (atomType & AtomPropertyMask.HYBRID.getMask())>>AtomPropertyShift.HYBRID_SHIFT.getShift();
		if(hybrid>0) {
			sb.append("." + hybrid);
		}
		int neighbours = (atomType & AtomPropertyMask.NONHNEIGHBOURS.getMask())>>AtomPropertyShift.NEIGHBOURS_SHIFT.getShift();
		sb.append("-" + neighbours + "-");
		String stab = (atomType & AtomPropertyMask.STABILIZED.getMask())>0 ? ".St" : "" ;
		sb.append(stab);
		String arom = (atomType & AtomPropertyMask.AROM.getMask())>0 ? ".Ar" : "" ;
		sb.append(arom);
		int fgValue = (atomType & AtomPropertyMask.FUNCTIONAL_GROUP.getMask())>>AtomPropertyShift.FUNCTIONAL_GROUP_SHIFT.getShift();
		for(FunctionalGroup fg : FunctionalGroup.values()) {
			if (fg.getId()==fgValue) {
				sb.append("(");
				sb.append(fg.getString());
				sb.append(")");
			}
		}
		return sb.toString();
		
	}

	public static int getAtomicNumber(int atomType){
		return atomType & AtomPropertyMask.ATOMIC_NO.getMask();
	}

	public static boolean isCarbonInteraction(int atomType){
		int atNo = getAtomicNumber(atomType);
		return (atNo==6)?true:false;
	}

	public static boolean isAromatic(int atomType){
		return ((atomType & AtomPropertyMask.AROM.getMask()) > 0) ? true : false;
	}


	public static void setInteractionTypes(Molecule3D mol) {
		for (int i = 0; i < mol.getAtoms(); i++) {
			int atomType = getAtomType(mol, i);
			mol.setInteractionAtomType(i, atomType);
		}
	}
	
	public static int getGenericDonor() {//represents a hydroxyl group
		return getAtomType(null, 8, false, 3,
			false);
	}
	
	public static int getGenericAcceptor() {//represents a hydroxyl group
		return getAtomType(null, 8, false, 3,
			false);
	}
	
	public static int getGenericPosCharge() {//represents sp3 primary amine
		return getAtomType(null, 7, false, 3,
			false);
	}
	
	public static int getGenericNegCharge() {//represents a carboxylate oxygen
		return getAtomType(FunctionalGroup.CARBOXYL, 8, false, 2,
			false);
	}
	
}
	
