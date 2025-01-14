package com.actelion.research.chem.io;

import java.io.BufferedReader;
import java.io.StringReader;
import java.util.BitSet;
import java.util.Map;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.chem.io.CDXMLParser.CDBond;
import com.actelion.research.chem.io.CDXMLParser.CDNode;

/**
 * A parser for binary CDX or XML CDXML files.
 * 
 * For the specification of CDX/CDXML, see
 * https://iupac.github.io/IUPAC-FAIRSpec
 * 
 * All the binary parser does is convert the CDX to CDXML and then feed that to
 * XML parser. (Seemed like a useful thing to have anyway, and it was simpler
 * this way.)
 * 
 * See https://iupac.github.io/IUPAC-FAIRSpec for the full, detailed
 * specification.
 * 
 * revvity site:
 * 
 * https://support.revvitysignals.com/hc/en-us/articles/4408233129748-Where-is-the-ChemDraw-SDK-located
 * 
 * Their link:
 * 
 * https://web.archive.org/web/20221209095323/https://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/
 * 
 * WayBack machine Overview:
 * 
 * https://web.archive.org/web/20240000000000*
 * /https://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx
 * 
 * Partial archives:
 * 
 * https://web.archive.org/web/20160911235313/http://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/index.htm
 * 
 * https://web.archive.org/web/20160310081515/http://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/
 * 
 * https://web.archive.org/web/20100503174209/http://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/
 * 
 * Unfortunately, there appears to be no single archive that has all the images,
 * and so some of those are missing.
 * 
 * Here we are just looking for simple aspects that could be converted to valid
 * 2D MOL files, SMILES, and InChI.
 * 
 * Fragments (such as CH2CH2OH) and "Nickname"-type fragments such as Ac and Ph,
 * are processed correctly. But their 2D representations are pretty nuts.
 * ChemDraw does not make any attempt to place these in reasonable locations.
 * That said, Jmol's 3D minimization does a pretty fair job, and the default is
 * to do that minimization.
 * 
 * 2024.12.19 The CDXMLParser now supports square bracket multiple x-[-y-]n-c notation.
 * 
 * @author hansonr@stolaf.edu
 * 
 *
 */
public class CDXParser extends XmlReader implements CDXMLParser.CDXReaderI {

	private static final int BOND_ORDER_NULL = 0;
	private static final int BOND_ORDER_STEREO_EITHER = -1;

	private CDXMLParser parser;

	/**
	 * settable to true to force coordinate cleaning
	 */
	private boolean doCleanCoordinates;

	public void setCleanCoordinates() {
		doCleanCoordinates = true;
	}
	
	public CDXParser() {
		parser = new CDXMLParser(this);
	}

	public static StereoMolecule parseFile(String path) {
		byte[] cdx = ParserUtils.getURLContentsAsBytes(path);
		StereoMolecule mol = new StereoMolecule();
		return (new CDXParser().parse(mol, cdx) ? mol : null);
	}

	private StereoMolecule mol;

	public boolean parse(StereoMolecule mol, String cdxml) {
		this.mol = mol;

		reader = new BufferedReader(new StringReader(cdxml));
		err = parseXML();
		if (err != null)
			return false;
		parser.finalizeParsing();
		createMolecule();
		return true;
	}

	public boolean parse(StereoMolecule mol, byte[] cdx) {
		return (get(mol, cdx) != null);
	}

	private StereoMolecule get(StereoMolecule mol, byte[] cdx) {
		if (cdx == null || cdx.length == 0)
			return null;
		try {
			String cdxml;
			if (cdx[0] == 0x56) {
				cdxml = CDX2CDXML.fromCDX(cdx);
				// dump(cdxml);
			} else {
				cdxml = new String(cdx, "utf-8");
			}
			return (new CDXParser().parse(mol, cdxml) 
					&& mol.getAllAtoms() > 0 ? mol : null);
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	@Override
	public void processStartElement(String localName, String nodeName) {
		parser.processStartElement(localName, atts);
	}

	@Override
	void processEndElement(String localName) {
		parser.processEndElement(localName, chars.toString());
	}

	@Override
	public void setKeepChars(boolean TF) {
		super.setKeepChars(TF);
	}

	public void handleCoordinates(Map<String, String> atts) {
		parser.setAtom("p", atts);
	}

	public int getBondOrder(String key) {
		switch (key) {
		case "1":
		case "single":
			return StereoMolecule.cBondTypeSingle;
		case "2":
		case "double":
			return StereoMolecule.cBondTypeDouble;
		case "3":
		case "triple":
			return StereoMolecule.cBondTypeTriple;
		case "up":
			return StereoMolecule.cBondTypeUp;
		case "down":
			return StereoMolecule.cBondTypeDown;
		case "either":
			return BOND_ORDER_STEREO_EITHER;
		case "null":
			return BOND_ORDER_NULL;
		case "delocalized":
			return StereoMolecule.cBondTypeDelocalized;
		default:
		case "partial":
			return StereoMolecule.cBondTypeMetalLigand;
		}
	}

	private void createMolecule() {
		BitSet bs = parser.bsAtoms;
		for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i + 1)) {
			CDNode a = parser.getAtom(i);
			int ia = mol.addAtom(a.x, -a.y);
			a.index = ia;
			mol.setAtomCharge(ia, a.formalCharge);
			mol.setAtomicNo(ia, a.elementNumber);
			if (a.isotope != null)
				mol.setAtomMass(ia, parser.parseInt(a.isotope));
		}
		bs = parser.bsBonds;
		for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i + 1)) {
			CDBond bond = parser.getBond(i);
			int ib = mol.addBond(bond.atom1.index, bond.atom2.index);
			mol.setBondType(ib, bond.order);
		}
		
		if (doCleanCoordinates || parser.hasConnectedFragments()) {
			mol.ensureHelperArrays(Molecule.cHelperParities);
			new CoordinateInventor(0).invent(mol);
		}
	}

	@Override
	public void warn(String warning) {
		System.err.println(warning);
	}

}
