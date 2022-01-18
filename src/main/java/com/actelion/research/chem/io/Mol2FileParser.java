/*
 * Created on Dec 20, 2004
 *
 */
package com.actelion.research.chem.io;

import com.actelion.research.chem.*;

import java.io.LineNumberReader;
import java.io.Reader;
import java.io.Writer;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

/**
 * 
 * @author freyssj
 */
public class Mol2FileParser extends AbstractParser {
	
	public static final int iNO_CHARGES = 0; 
	public static final int iDEL_RE = 1;
	public static final int iGASTEIGER = 2;
	public static final int iGAST_HUCK = 3;
	public static final int iHUCKEL = 4;
	public static final int iPULLMAN = 5;
	public static final int iGAUSS80_CHARGES = 6; 
	public static final int iAMPAC_CHARGES = 7;
	public static final int iMULLIKEN_CHARGES = 8;
	public static final int iDICT_CHARGES = 9;
	public static final int iMMFF94_CHARGES = 10;
	public static final int iUSER_CHARGES = 11;
	
	private static final NumberFormat NF_PARTIAL_CHARGES = new DecimalFormat(" 0.0000;-0.0000");

	private static final String sNO_CHARGES = "NO_CHARGES"; 
	private static final String sDEL_RE = "DEL_RE";
	private static final String sGASTEIGER = "GASTEIGER";
	private static final String sGAST_HUCK = "GAST_HUCK";
	private static final String sHUCKEL = "HUCKEL";
	private static final String sPULLMAN = "PULLMAN";
	private static final String sGAUSS80_CHARGES = "GAUSS80_CHARGES"; 
	private static final String sAMPAC_CHARGES = "AMPAC_CHARGES";
	private static final String sMULLIKEN_CHARGES = "MULLIKEN_CHARGES";
	private static final String sDICT_CHARGES = "DICT_CHARGES";
	private static final String sMMFF94_CHARGES = "MMFF94_CHARGES";
	private static final String sUSER_CHARGES = "USER_CHARGES";
	
	private static HashMap<Integer, String> hmIndex_CHARGETYPE;
	
	private boolean isLoadHydrogen = true;

	private void addMol(List<Molecule3D> res, Molecule3D m, Set<Integer> aromaticAtoms, Set<Integer> aromaticBonds) {
		if(m==null || m.getAllAtoms()==0) return;
		
		if(true) {
			//Use AromaticityResolver		
			ExtendedMolecule mol = new StereoMolecule(m);
			int[] bondMap = mol.getHandleHydrogenBondMap();
			int j =0;
			//for (int b : aromaticBonds) {
			//	if(mol.getAtomicNo(mol.getBondAtom(0, b))!=8 && mol.getAtomicNo(mol.getBondAtom(1, b))!=8) {
			//		mol.setBondType(b, ExtendedMolecule.cBondTypeDelocalized);
			//	}
			//}
			new AromaticityResolver(mol).locateDelocalizedDoubleBonds(null,true,true);
			
	
			for (int i=0; i<mol.getBonds(); i++) {
				m.setBondOrder(i, mol.getBondOrder(bondMap[i]));
			}
		}
		/*
		else {
			//Use BondsCalculator.aromatize
			boolean suc = BondsCalculator.aromatize(m, aromaticAtoms, aromaticBonds);
			if(!suc) {
				System.err.println("Could not aromatize the molecule");
			}
		//}		
		 * 
		 */

		assignCharges(m);

		m.setAllAtomFlag(Molecule3D.LIGAND, true);

		res.add(m);
	}
	
	private void assignCharges(Molecule3D mol) {
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (mol.getAtomicNo(atom) == 7)
				if(mol.getOccupiedValence(atom)==4) {
					//sulfonium ion
					mol.setAtomCharge(atom, 1);
				}
				else if(mol.getOccupiedValence(atom)==2) {
					mol.setAtomCharge(atom, -1);
				}
			
			if (mol.getAtomicNo(atom) == 8)
				if(mol.getOccupiedValence(atom)==3) {
					//sulfonium ion
					mol.setAtomCharge(atom, 1);
				}
				else if(mol.getOccupiedValence(atom)==1) {
					mol.setAtomCharge(atom, -1);
				}
			if (mol.getAtomicNo(atom) == 16)
				if(mol.getOccupiedValence(atom)==3) {
					//sulfonium ion
					mol.setAtomCharge(atom, 1);
				}
				else if(mol.getOccupiedValence(atom)==1)
					mol.setAtomCharge(atom, -1);
			
		}
	}
	
	
	public static String getChargeType(int type) {
		if(hmIndex_CHARGETYPE==null){
			hmIndex_CHARGETYPE = new HashMap<Integer, String>();
			hmIndex_CHARGETYPE.put(iNO_CHARGES, sNO_CHARGES);
			hmIndex_CHARGETYPE.put(iDEL_RE, sDEL_RE);
			hmIndex_CHARGETYPE.put(iGASTEIGER, sGASTEIGER);
			hmIndex_CHARGETYPE.put(iGAST_HUCK, sGAST_HUCK);
			hmIndex_CHARGETYPE.put(iHUCKEL, sHUCKEL);
			hmIndex_CHARGETYPE.put(iPULLMAN, sPULLMAN);
			hmIndex_CHARGETYPE.put(iGAUSS80_CHARGES, sGAUSS80_CHARGES);
			hmIndex_CHARGETYPE.put(iAMPAC_CHARGES, sAMPAC_CHARGES);
			hmIndex_CHARGETYPE.put(iMULLIKEN_CHARGES, sMULLIKEN_CHARGES);
			hmIndex_CHARGETYPE.put(iDICT_CHARGES, sDICT_CHARGES);
			hmIndex_CHARGETYPE.put(iMMFF94_CHARGES, sMMFF94_CHARGES);
			hmIndex_CHARGETYPE.put(iUSER_CHARGES, sUSER_CHARGES);
		}
		
		return hmIndex_CHARGETYPE.get(type);
	}
	
	
	/**
	 * @see com.actelion.research.chem.io.AbstractParser#loadGroup(String, java.io.Reader, int, int)
	 */
	@Override
	public List<Molecule3D> loadGroup(String fileName, Reader in, int from, int to) throws Exception {
		List<Molecule3D> res = new ArrayList<Molecule3D>();
		String line;

		Molecule3D m = new Molecule3D();
		Set<Integer> aromaticAtoms = new HashSet<Integer>();
		Set<Integer> aromaticBonds = new HashSet<Integer>();
		Map<Integer, Integer> nToA = new TreeMap<Integer, Integer>();
		LineNumberReader reader = new LineNumberReader(in);
		
		int state = 0;
		int lineNo = 0;
		int count = -1;
		while ((line = reader.readLine()) != null) {
			if (line.startsWith("@<TRIPOS>")) lineNo = 0;
			
			if (line.startsWith("@<TRIPOS>MOLECULE")) state = 0;
			else if (line.startsWith("@<TRIPOS>ATOM")) state = 1;
			else if (line.startsWith("@<TRIPOS>BOND")) state = 2;
			else if (line.startsWith("@<TRIPOS>")) state = -1;
			else {
				lineNo++;
				switch (state) {
				case 0: {
					if(lineNo==1) {
						//We have a new molecule, first process the last molecule and add it to the list
						count++;
						if(from>count) continue;
						if(to>=0 && to<count) break;
						addMol(res, m, aromaticAtoms, aromaticBonds);
						//Now create a new molecule
						m = new Molecule3D();
						aromaticAtoms.clear();
						aromaticBonds.clear();
						String molName = line.trim();
						m.setName(molName);
					}
					break;
				}
				case 1: {
					try {

						StringTokenizer st = new StringTokenizer(line, "\t ");
						int n = Integer.parseInt(st.nextToken().trim());
						String atomName = st.nextToken().trim();
						double x = Double.parseDouble(st.nextToken().trim());
						double y = Double.parseDouble(st.nextToken().trim());
						double z = Double.parseDouble(st.nextToken().trim());
						String atomClass = st.hasMoreTokens()? st.nextToken().trim(): "";
						String chainId = st.hasMoreTokens()? st.nextToken().trim(): "";
						String amino = st.hasMoreTokens()? st.nextToken().trim(): "";
						String charge = st.hasMoreTokens()? st.nextToken().trim(): "";
						/*
						while(line.length()<90) line += "                       ";
						int n = Integer.parseInt(line.substring(0, 8).trim());
						String atomName = line.substring(8, 12).trim();
						double x = Double.parseDouble(line.substring(17, 26).trim());
						double y = Double.parseDouble(line.substring(27, 36).trim());
						double z = Double.parseDouble(line.substring(37, 46).trim());
						String atomClass = line.substring(46, 53).trim();
						String chainId = line.substring(54, 58).trim();
						String amino = line.substring(59, 69).trim();
						String charge = line.substring(69, 81).trim();
						*/
						
						String elt = new StringTokenizer(atomClass, ".").nextToken();
						boolean aromatic = atomClass.endsWith(".ar");
						int atomicNo = Molecule.getAtomicNoFromLabel(elt);
						if (atomicNo < 0)
							throw new Exception("Invalid Atomic Number for " + n + ": " + elt);
	
						if (!isLoadHydrogen && atomicNo <= 1)
							continue;
	
						int a = m.addAtom(atomicNo);
						m.setAtomX(a, x);
						//invert y and z coordinates for compatibility with Java coordinate system (analogously to Molfileparser)
						m.setAtomY(a, -y);
						m.setAtomZ(a, -z);
						m.setAtomName(a, atomName);
						m.setAtomChainId(a, chainId);
						m.setAtomAmino(a, amino);
						if(aromatic) aromaticAtoms.add(a);
						try {
							m.setPartialCharge(a, Double.parseDouble(charge));
						} catch (Exception e) {
						}
						
	
						nToA.put(n, a);
					} catch (NumberFormatException e) {
						throw new Exception("Invalid number at line "+reader.getLineNumber()+": "+line);
					}
					break;
				}
				case 2: {
					StringTokenizer st = new StringTokenizer(line);
					if(st.countTokens()<3) continue;
					st.nextToken(); // bondNo
					int n1 = Integer.parseInt(st.nextToken());
					int n2 = Integer.parseInt(st.nextToken());
					String o = st.nextToken();
					
					
					Integer i1 = nToA.get(n1);
					Integer i2 = nToA.get(n2);

					if (i1 == null || i2 == null) continue;
					
					int order;
					if (o.equals("ar")) {
						//order = 1;
						//aromatic.add(i1);
						//aromatic.add(i2);
						order = Molecule.cBondTypeDelocalized;
					} else if (o.equals("am")) {
						order = Molecule.cBondTypeSingle;
						//order = 1;
					} else if(o.equals("un")) {
						continue;
					} else if(o.equals("nc")) {
						continue;
					} else if(o.equals("du")) {
						continue;
					} else {
						order = Integer.parseInt(o);
						
						if (order== 1) {
							order = Molecule.cBondTypeSingle;
						} else if (order == 2) { 
							order = Molecule.cBondTypeDouble;
						} else if (order == 3) {
							order = Molecule.cBondTypeTriple;
						} else {
							throw new RuntimeException("Unknown bond type " + order + ".");
						}
					}

					int a1 = i1.intValue();
					int a2 = i2.intValue();
					int b = m.addBond(a1, a2,order);

					//m.setBondOrder(b, order);
					//if(o.equals("ar")) {
					//	aromaticBonds.add(b);
					//}
					
					break;
				}
				}
			}
		}
		count++;
		int[] bondMap = m.getHandleHydrogenBondMap();
		int[] atomMap = m.getHandleHydrogenMap();
		Set<Integer> aromaticAtomsMapped = new HashSet<Integer>();
		Set<Integer> aromaticBondsMapped = new HashSet<Integer>();
		aromaticAtoms.stream().forEach(aa -> aromaticAtomsMapped.add(atomMap[aa]));
		aromaticBonds.stream().forEach(ab -> aromaticBondsMapped.add(bondMap[ab]));
		if(from>count) {}
		else if(to>=0 && to<count) {}
		else {
			m.ensureHelperArrays(Molecule.cHelperRings);
			addMol(res, m, aromaticAtoms, aromaticBonds);
		}
		return res;
	}

	@Override
	public void save(Molecule3D mol, Writer writer) throws Exception {
		save(Collections.singletonList(mol), writer);
	}
	
	@Override
	public void save(List<Molecule3D> mols, Writer writer) throws Exception {
		save(mols, iNO_CHARGES, writer);
	}

	/**
	 * Writes also the partial charges
	 * @param mols
	 * @param chargeType
	 * @param writer
	 * @throws Exception
	 */
	public void save(List<Molecule3D> mols, int chargeType, Writer writer) throws Exception {
		DecimalFormat df = new DecimalFormat("0.0000");

		for (Molecule3D mol : mols) {
			
			//Count the number of bonds and atoms (except lone pairs)
			int nAtoms = 0;
			int nBonds = 0;
			for (int i = 0; i < mol.getAllAtoms(); i++) {
				if(mol.getAtomicNo(i)<=0) continue;
				nAtoms++;
			}
			for (int i = 0; i < mol.getAllBonds(); i++) {
				if(mol.getAtomicNo(mol.getBondAtom(0, i))<=0) continue;
				if(mol.getAtomicNo(mol.getBondAtom(1, i))<=0) continue;
				nBonds++;
			}

			// Write the molecule
			writer.write("@<TRIPOS>MOLECULE" + NEWLINE);
			writer.write(mol.getName() + NEWLINE);
			writeR(writer, "" + nAtoms, 5);
			writeR(writer, "" + nBonds, 6);
			writer.write(NEWLINE);
			writer.write("SMALL" + NEWLINE);
			writer.write(getChargeType(chargeType) + NEWLINE);
	
			// Write the atoms
			writer.write("@<TRIPOS>ATOM" + NEWLINE); 			
			for (int i = 0; i < mol.getAllAtoms(); i++) {
				if(mol.getAtomicNo(i)<=0) continue;
				writeR(writer, "" + (i+1), 7);
				writer.write(" ");
				
				String desc = mol.getAtomName(i);
				if(desc==null || desc.length()==0) {
					desc = "" + ExtendedMolecule.cAtomLabel[mol.getAtomicNo(i)];
				}
				
				writeL(writer, desc, 8);
				writeR(writer, df.format(mol.getAtomX(i)), 10);
				writeR(writer, df.format(mol.getAtomY(i)), 10);
				writeR(writer, df.format(mol.getAtomZ(i)), 10);
				writer.write(" ");
				writeL(writer, ExtendedMolecule.cAtomLabel[mol.getAtomicNo(i)] + (mol.isAromaticAtom(i)?".ar":""), 8);
				writeR(writer, mol.getAtomChainId(i)!=null && mol.getAtomChainId(i).length()>0? mol.getAtomChainId(i): "1", 5);
				writer.write("  ");
				writeL(writer, mol.getAtomAmino(i)!=null && mol.getAtomAmino(i).length()>0? mol.getAtomAmino(i): "1", 11);
				writer.write(NF_PARTIAL_CHARGES.format(mol.getPartialCharge(i)));
	
				writer.write(NEWLINE);
			}
	
			// Write the bonds
			writer.write("@<TRIPOS>BOND" + NEWLINE);
			for (int i = 0; i < mol.getAllBonds(); i++) {
				if(mol.getAtomicNo(mol.getBondAtom(0, i))<=0) continue;
				if(mol.getAtomicNo(mol.getBondAtom(1, i))<=0) continue;
				writeR(writer, "" + (i+1), 6);
				writeR(writer, "" + (mol.getBondAtom(0, i)+1), 5);
				writeR(writer, "" + (mol.getBondAtom(1, i)+1), 5);
				writer.write(" ");
				writeL(writer, "" + mol.getBondOrder(i), 5);
				writer.write(NEWLINE);
			}
		}
	}

	
	public static void main(String[] args) throws Exception {
		List<Molecule3D> res = new Mol2FileParser().loadGroup("c:/mopac11971.mol2");
		//List<Molecule3D> res = new Mol2FileParser().loadGroup("c:/orig_pose_min.mol2");
		System.out.println("models=" +res.size()+" atm="+res.get(0).getAllAtoms()+" "+res.get(0).getAtoms());
		//new Mol2FileParser().save(res.get(0), "c:/t.mol2");
	}
}
