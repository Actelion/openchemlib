/* $RCSfile$
 * $Author: hansonr $
 * $Date: 2006-08-02 11:48:43 -0500 (Wed, 02 Aug 2006) $
 * $Revision: 5364 $
 *
 * Copyright (C) 2003-2005  Miguel, Jmol Development, www.jmol.org
 *
 * Contact: jmol-developers@lists.sf.net
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package com.actelion.research.chem.io;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import com.actelion.research.chem.io.XmlReader.XmlHandler;

/**
 * A reader for CambridgeSoft CDXML files.
 * 
 * Copied from org.jmol.adapter.readers.xml.CDXMLParser.java
 * 
 * For the full specification of CDX/CDXML, see
 * https://iupac.github.io/IUPAC-FAIRSpec
 * 
 * revvity site:
 * 
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
 * If minimization and addition of H is not desired, use FILTER "NOH" or FILTER
 * "NO3D"
 * 
 * XmlChemDrawReader also serves as the reader for binary CDX files, as
 * CDXReader subclasses this class. See that class for details.
 * 
 * @author hansonr@stolaf.edu
 * 
 */

class CDXMLParser extends XmlHandler  {

  public interface CDXReaderI {

    int getBondOrder(String string);

    void handleCoordinates(Map<String, String> atts);

    void setKeepChars(boolean b);

    void warn(String string);

  }

  private double minX = Double.MAX_VALUE;
  private double minY = Double.MAX_VALUE;
  private double minZ = Double.MAX_VALUE;
  private double maxZ = -Double.MAX_VALUE;
  private double maxY = -Double.MAX_VALUE;
  private double maxX = -Double.MAX_VALUE;

  private int idnext = 100000;
  
  protected BitSet bsAtoms = new BitSet(), bsBonds = new BitSet();

  /**
   * main list
   */
  List<CDNode> atoms = new ArrayList<CDNode>();

  /**
   * main list
   */
  List<CDBond> bonds = new ArrayList<CDBond>();

  Map<String, CDBond> bondIDMap = new HashMap<String, CDBond>();

  private Stack<BracketedGroup> bracketedGroups;

  protected CDXReaderI rdr;
  public Map<String, CDNode> mapCloned;

  /**
   * CDNode extends Atom in order to maintain information about fragments,
   * connectivity, and validity
   * 
   */
  class CDNode implements Cloneable {

    public int index;
    public String id;
    public int intID;
    public String isotope;
    public String element;
    public int elementNumber;
    public double x, y, z;
    public int formalCharge;
    
    String nodeType;
    String warning;

    boolean isValid = true;
    boolean isConnected;
    boolean isExternalPt;
    boolean isFragment; // could also be a Nickname
    /**
     * fragment ID for the fragment containing this node
     */
    String outerFragmentID;

    /**
     * fragment ID of this fragment node
     */
    String innerFragmentID;
    
    /**
     * text associated with this node
     */
    public String text;
    
    /**
     * fragment atom parent node
     */
    CDNode parentNode;
    /**
     * list of connection bonds, ordered by ID
     */
    List<Object[]> orderedConnectionBonds;
    /**
     * for an external point, the actual atom associated with it in the fragment  
     */
    CDNode internalAtom;
    /**
     * for a fragment, the list of external points for a fragment, ordered by sequence in the label
     */
    List<CDNode> orderedExternalPoints;

    /**
     * 0x0432 For multicenter attachment nodes or variable attachment nodes a
     * list of IDs of the nodes which are multiply or variably attached to this
     * node array of attachment id values;
     * 
     * for example, in ferrocene, we are attaching to all of the carbon atoms if
     * this node is the special point that indicates that attachment
     */
    private String[] attachments;
    /**
     * 0x0431 BondOrdering An ordering of the bonds to this node used for
     * stereocenters fragments and named alternative groups with more than one
     * attachment.
     * 
     */
    private String[] bondOrdering;
    /**
     * 0x0505 ConnectionOrder An ordered list of attachment points within a
     * fragment.
     * 
     */
    private String[] connectionOrder;

    public boolean hasMultipleAttachments;
    CDNode attachedAtom;

    /**
     * atom-atom edges, including special points
     */
    public BitSet bsConnections;

    /**
     * atoms of a nickname fragment
     */
    BitSet bsFragmentAtoms;
    int clonedIndex= -1;

    CDNode(int index, String id, String nodeType, String fragmentID,
        CDNode parent) {
      this.id = id;
      this.index = index;
      this.outerFragmentID = fragmentID;
      this.intID = Integer.parseInt(id);
      this.nodeType = nodeType;
      this.parentNode = parent;
      this.bsConnections = new BitSet();
      isFragment = "Fragment".equals(nodeType) || "Nickname".equals(nodeType);
      isExternalPt = "ExternalConnectionPoint".equals(nodeType);
      if (isFragment) {
        bsFragmentAtoms = new BitSet();
      } else if (parent != null && !isExternalPt) {
        parent.bsFragmentAtoms.set(index);
      }
      // isGeneric = "GenericNickname".equals(nodeType);
    }

		public void set(double x, double y, double z) {
			this.x = x;
			this.y = y;
			this.z = z;
		}

    public void setInnerFragmentID(String id) {
      innerFragmentID = id;
    }

    void setBondOrdering(String[] bondOrdering) {
      this.bondOrdering = bondOrdering;
    }

    void setConnectionOrder(String[] connectionOrder) {
      this.connectionOrder = connectionOrder;
    }

    void setMultipleAttachments(String[] attachments) {
      this.attachments = attachments;
      hasMultipleAttachments = true;
    }

    /**
     * keep these in order
     * 
     * @param externalPoint
     */
    void addExternalPoint(CDNode externalPoint) {
      if (orderedExternalPoints == null)
        orderedExternalPoints = new ArrayList<CDNode>();
      int i = orderedExternalPoints.size();
      while (--i >= 0 && orderedExternalPoints.get(i).intID >= externalPoint.internalAtom.intID) {
        // continue;
      }
      orderedExternalPoints.add(++i, externalPoint);
    }

    public void setInternalAtom(CDNode a) {
      internalAtom = a;
			if (parentNode != null) {
        parentNode.addExternalPoint(this);
      }
    }

    void addAttachedAtom(CDBond bond, int pt) {
      if (orderedConnectionBonds == null)
        orderedConnectionBonds = new ArrayList<Object[]>();
      int i = orderedConnectionBonds.size();
      while (--i >= 0
          && ((Integer) orderedConnectionBonds.get(i)[0]).intValue() > pt) {
        // continue;
      }
      orderedConnectionBonds.add(++i, new Object[] { Integer.valueOf(pt), bond });
    }

    void fixAttachments() {
      if (hasMultipleAttachments && attachedAtom != null) {
        // something like Ferrocene
        int order = rdr.getBondOrder("partial");
        CDNode a1 = attachedAtom;
        for (int i = attachments.length; --i >= 0;) {
          CDNode a = (CDNode) objectsByID.get(attachments[i]);
          if (a != null) {
            addBond(new CDBond(null, a1.id, a.id, order));
          }
        }
      }

      if (orderedExternalPoints == null || text == null)
        return;
      // fragments and Nicknames
      int n = orderedExternalPoints.size();
      if (n != orderedConnectionBonds.size()) {
        System.err.println(
            "XmlCdxReader cannot fix attachments for fragment " + text);
        return;
      }
      if (bondOrdering == null) {
        bondOrdering = new String[n];
        for (int i = 0; i < n; i++) {
          bondOrdering[i] = ((CDBond)orderedConnectionBonds.get(i)[1]).id;
        }
      }
      if (connectionOrder == null) {
        connectionOrder = new String[n];      
        for (int i = 0; i < n; i++) {
          connectionOrder[i] = orderedExternalPoints.get(i).id;
        }
      }
      
        for (int i = 0; i < n; i++) {
          CDBond b = (CDBond) objectsByID.get(bondOrdering[i]);
          CDNode a = ((CDNode) objectsByID.get(connectionOrder[i])).internalAtom;
          updateExternalBond(b, a);
        }
    }

    /**
     * Replace the fragment connection (to this fragment node) in bond b with
     * the internal atom a.
     * 
     * @param bond2f
     * @param intAtom
     */
    private void updateExternalBond(CDBond bond2f, CDNode intAtom) {
			bsBonds.set(bond2f.index);
			bond2f.disconnect();
      if (bond2f.atomIndex2 == index) {
				bond2f.connect(bond2f.atom1, intAtom);
      } else if (bond2f.atomIndex1 == index) {
				bond2f.connect(intAtom, bond2f.atom2);
      } else {
				System.err.println("CDXMLParser attachment failed! " + intAtom + " " + bond2f);
      }
      
    }

		@Override
    public CDNode clone() {
			CDNode a;
			try {
				a = (CDNode) super.clone();
				a.index = atoms.size();
				a.id = nextID();
				mapCloned.put(id, a);
				a.clonedIndex = index;
        a.bsConnections = new BitSet();
				objectsByID.put(a.id, a);
				// for minimization, don't put these on top of each other.

//				a.x += Math.random()* 0.1 - 0.05;
//        a.y += Math.random()* 0.1 - 0.05;
//        a.z += Math.random()* 0.1 - 0.05;
				addAtom(a); 
				return a;
			} catch (CloneNotSupportedException e) {
				return null;
			}
		}
    @Override
    public String toString() {
      return "[CDNode " + id + " " + elementNumber + " index=" + index
          + " ext=" + isExternalPt + " frag=" + isFragment + " "
          + " " + x + " " + y +"]";
    }

    public double distance(CDNode a2) {
      double dx = (x - a2.x);
      double dy = (y - a2.y);
      return Math.sqrt(dx * dx + dy * dy);
    }

    public void addBracketedAtoms(BitSet bsBracketed) {
      if (isFragment)
        bsBracketed.or(bsFragmentAtoms);
      else if (!isExternalPt)
        bsBracketed.set(index);
     }

  }

  class CDBond implements Cloneable {
    public int atomIndex1, atomIndex2;
    public int order;
    String id, id1, id2;

    CDNode atom1;
    CDNode atom2; 

    int index;

    boolean invalidated;

    CDBond(String id, String id1, String id2, int order) {
      this.id = (id == null ? nextID() : id);
      this.id1 = id1;
      this.id2 = id2;
      this.order = order;
      atom1 = (CDNode) objectsByID.get(id1);
      atom2 = (CDNode) objectsByID.get(id2);
      atomIndex1 = atom1.index;
      atomIndex2 = atom2.index;
			atom1.bsConnections.set(atomIndex2);
			atom2.bsConnections.set(atomIndex1);
			bondIDMap.put(getBondKey(atomIndex1, atomIndex2), this);
    }

		public void connect(CDNode atomA, CDNode atomB) {
			atom1 = atomA;
			atomIndex1 = atomA.index;
			atom2 = atomB;
			atomIndex2 = atomB.index;
			atomA.bsConnections.set(atomB.index);
			atomB.bsConnections.set(atomA.index);
      bondIDMap.put(getBondKey(atomIndex1, atomIndex2), this);
		}

		public void disconnect() {
			atom1.bsConnections.clear(atomIndex2);
			atom2.bsConnections.clear(atomIndex1);
			bondIDMap.remove(getBondKey(atomIndex1, atomIndex2));
		}

    CDNode getOtherNode(CDNode a) {
      return getAtom(atomIndex1 == a.index ? atomIndex2
          : atomIndex1);
    }

    public void invalidate() {
      invalidated = true;
      bsBonds.clear(index);
      atomIndex1 = atomIndex2 = -1;
    }
    
    public boolean isValid() {
      return (!invalidated && atom1.intID >= 0 && atom2.intID >= 0);
    }

    @Override
    public String toString() {
      return "[CDBond index " + index +" id " + id 
          + " v=" + isValid()
          + " id1=" + id1 + "/" + atom1.index +"/" + atom1.clonedIndex 
          + " id2=" + id2 + "/" + atom2.index +"/" + atom2.clonedIndex + "]";
    }

  }

  class BracketedGroup {
	    String id;
	    String[] ids;
	    String[] bondIDs = new String[2];
	    String[] innerAtomIDs = new String[2];
	    int repeatCount;
	    int pt;

	    BracketedGroup(String id, String[] ids, int repeatCount) {
	      this.id = id;
	      this.ids = ids;
	      this.repeatCount = repeatCount;
	    }

    public void addCrossing(String innerAtomID, String bondID) {
	      if (pt == 2) {
	        System.err.println("BracketedGroup has more than two crossings");
	        return;
	      }
	      bondIDs[pt] = bondID;
	      innerAtomIDs[pt] = innerAtomID;
	      pt++;
	    }

	    public void process() {
	      if (pt != 2 || repeatCount < 2)
	        return;
	      // a1o-|-a1i...a2i-|-a2o...[trailing atoms]
	      CDBond b1 = (CDBond) objectsByID.get(bondIDs[1]);
        CDNode a1i = (CDNode) objectsByID.get(innerAtomIDs[1]);
        CDNode a2i = (CDNode) objectsByID.get(innerAtomIDs[0]);
	      CDBond b2 = (CDBond) objectsByID.get(bondIDs[0]);
	      CDNode a2o = b2.getOtherNode(a2i);
	      
	      BitSet bsTrailingAtoms = new BitSet();
	      buildBranch(a2i, a2o, null, null, bsTrailingAtoms, null, null);
	      double[] offset = new double[] {a2o.x - a1i.x, a2o.y - a1i.y, a2o.z - a1i.z }; 
	      // a1o-b1-a1i...a2i-b2-a2o
	      // ....|------------|
	      // becomes
	      // a1o-b1-a1i...a2i-b1-a1i'...a2i'-b2-a2o
        // ....|------------|--------------|
	      BitSet bsBracketed = new BitSet();
	      for (int i = ids.length; --i >= 0;) {
	        CDNode node = (CDNode) objectsByID.get(ids[i]);
	        node.addBracketedAtoms(bsBracketed);
	      }
	      CDNode a1 = a1i;
	      CDNode a2iLast = a2i;
	      for (int i = 1; i < repeatCount; i++) {
	        int nAtoms = getAtomCount();
	        CDNode[] newNodes = duplicateBracketAtoms(a1i, a2i, bsBracketed);
          a1 = newNodes[0];
	        patch(a2iLast, b1, a1);
	        a2iLast = newNodes[1];
	        shiftAtoms(nAtoms, offset, i, null);
	      }
	      patch(a2iLast, b2, a2o);
	      shiftAtoms(0, offset, repeatCount - 1, bsTrailingAtoms);
        b2.invalidate();     
	    }

    private void shiftAtoms(int firstAtom, double[] d, int multiplier,
                            BitSet bs) {
      if (bs == null) {
        for (int i = getAtomCount(); --i >= firstAtom;) {
          CDNode a = getAtom(i);
          a.x += d[0] * multiplier;
          a.y += d[1] * multiplier;
          a.z += d[2] * multiplier;
        }
      } else {
        for (int i = bs.nextSetBit(0); i >=0; i = bs.nextSetBit(i + 1) ) {
          CDNode a = getAtom(i);
          a.x += d[0] * multiplier;
          a.y += d[1] * multiplier;
          a.z += d[2] * multiplier;         
        }
      }
    }


      private CDNode[] duplicateBracketAtoms(CDNode a1i, CDNode a2i, BitSet bsBracketed) {
	      mapCloned = new HashMap<>();     CDNode aNew = a1i.clone();
	      CDNode[] a12i = new CDNode[] {aNew, a1i == a2i ? aNew : null};
        BitSet bsDone = new BitSet();
	      buildBranch(null, a1i, a2i, aNew, bsDone, bsBracketed, a12i);
	      return a12i;
	    }

    private void buildBranch(CDNode prev, CDNode aRoot, CDNode aEnd, CDNode a,
                             BitSet bsDone, BitSet bsBracketed, CDNode[] a12i) {
      bsDone.set(aRoot.index);
      CDNode aNext = null;
      for (int i = aRoot.bsConnections
          .nextSetBit(0); i >= 0; i = aRoot.bsConnections.nextSetBit(i + 1)) {
        CDNode aBranch = getAtom(i);
        if (aBranch == prev)
          continue;
        boolean isNew = !bsDone.get(i);
        if (bsBracketed != null) {
          if (!bsBracketed.get(i) || aBranch.isExternalPt)
            continue;
          CDBond bBranch = getBond(aRoot, aBranch);
          aNext = (isNew ? aBranch.clone() : mapCloned.get(aBranch.id));
          if (aBranch == aEnd && a12i[1] == null) {
            a12i[1] = aNext;
          }
          CDBond b = (bBranch.atom1 == aBranch
               ? new CDBond(null, aNext.id, a.id, bBranch.order)
               : new CDBond(null, a.id, aNext.id, bBranch.order)
          );
          addBond(b);
        }
        if (isNew) {
          buildBranch(aRoot, aBranch, aEnd, aNext, bsDone, bsBracketed, a12i);
        }
      }
    }

	    private void patch(CDNode a1, CDBond b, CDNode a2) {
	      CDBond b1 = addBond(new CDBond(null, a1.id, a2.id, b.order));
	      b1.disconnect();
	      b1.connect(a1, a2);
	    }
	  }

  /**
   * not public
   * @param reader
   */
  CDXMLParser(CDXReaderI reader) {
    this.rdr = reader;
  }

  private Stack<String> fragments = new Stack<String>();
  private String thisFragmentID;
  private CDNode thisNode;
  private CDNode thisAtom;
  private boolean ignoreText;
  private Stack<CDNode> nodes = new Stack<CDNode>();
  private List<CDNode> nostereo = new ArrayList<CDNode>();
  Map<String, Object> objectsByID = new HashMap<String, Object>();
  
  /**
   * set true if connected fragments are found, which may need cleaning
   */
  private boolean haveConnectedFragments;


  /**
   * temporary holder of style chunks within text objects
   */
  private String textBuffer;

  public void processStartElement(String localName, Map<String, String> atts) {
    String id = atts.get("id");
    switch (localName) {
    case "n":
      objectsByID.put(id, setNode(id, atts));
      break;
    case "b":
      objectsByID.put(id, setBond(id, atts));
      break;
    case "t":
      textBuffer = "";
      break;
    case "s":
      rdr.setKeepChars(true);
      break;
    case "fragment":
      objectsByID.put(id, setFragment(id, atts));
      break;
    case "objecttag":
      // avoid messages about unidentified text for square brackets
      switch(atts.get("name")) {
      case "parameterizedBracketLabel":
      case "bracketusage":
        ignoreText = true;
        break;
      }
      break;
    case "bracketedgroup":
      setBracketedGroup(id, atts);
      break;
    case "crossingbond":
      BracketedGroup bg = (bracketedGroups == null || bracketedGroups.isEmpty() ? null
          : bracketedGroups.get(bracketedGroups.size() - 1));
      if (bg != null && bg.repeatCount > 0) {
        bg.addCrossing(atts.get("inneratomid"), atts.get("bondid"));
      }
      break;
    }
  }
  
  public String nextID() {
    return "" + (idnext++);
  }

  public static String getBondKey(int atomIndex1, int atomIndex2) {
    return Math.min(atomIndex1,  atomIndex2) + "_" + Math.max(atomIndex1, atomIndex2);
  }

  public CDBond getBond(CDNode a, CDNode b) {
    return bondIDMap.get(getBondKey(a.index, b.index));
  }

  void processEndElement(String localName, String chars) {
    switch (localName) {
    case "fragment":
      thisFragmentID = fragments.pop();
      return;
    case "objecttag":
      ignoreText = false;
      return;
    case "n":
      thisNode = (nodes.size() == 0 ? null : nodes.pop());
      return;
    case "bracketedgroup":
      break;
    case "s":
      textBuffer += chars.toString();
      break;
    case "t":
      if (ignoreText) {
      } else if (thisNode == null) {
        System.out.println("CDXReader unassigned text: " + textBuffer);
      } else {
        thisNode.text = textBuffer;
        if (thisAtom.elementNumber == 0) {
          System.err.println(
              "XmlChemDrawReader: Problem with \"" + textBuffer + "\"");
        }
        if (thisNode.warning != null)
          rdr.warn("Warning: " + textBuffer + " " + thisNode.warning);
      }
      textBuffer = "";
      break;
    }

    rdr.setKeepChars(false);
  }


  private String[] getTokens(String s) {
    return s.split("\\s");
  }

  private String[] split(String s, String p) {
    return s.split(p);
  }

  int parseInt(String s) {
    try {
      return Integer.parseInt(s);
    } catch (Exception e) {
      return Integer.MIN_VALUE;
    }
  }

  private double parseDouble(String s) {
    try {
      return Double.parseDouble(s);
    } catch (Exception e) {
      return Double.NaN;
    }
  }

	/**
	 * Set the atom information. Reading:
	 * 
	 * NodeType, Warning. Element, Isotope, Charge, xyz, p
	 * 
	 * 3D coordinates xyz is only used if there are no 2D p coordinates. This may
	 * not be possible. I don't know. These aren't real 3D coordinates, just
	 * enhanced z values.
	 * 
	 * @param id
	 * @param atts
	 * @return thisNode
	 */
	private CDNode setNode(String id, Map<String, String> atts) {
		String nodeType = atts.get("nodetype");
		if (thisNode != null)
			nodes.push(thisNode);
		if (nodeType != null) {
			switch (nodeType) {
			case "_":
				// internal Jmol code for ignored node
				thisAtom = thisNode = null;
				return null;
			case "ExternalConnectionPoint":
				haveConnectedFragments = true;
				break;
			}
		}
		thisAtom = thisNode = new CDNode(atoms.size(), id, nodeType, thisFragmentID, thisNode);
		addAtom(thisNode);

		String w = atts.get("warning");
		if (w != null) {
			thisNode.warning = w.replace("&apos;", "'");
			thisNode.isValid = (w.indexOf("ChemDraw can't interpret") < 0);
		}

		String element = atts.get("element");
		String s = atts.get("genericnickname");
		if (s != null) {
			element = s;
		}
		thisAtom.element = element;
		thisAtom.elementNumber = (short) Math.max(0,
				(!checkWarningOK(w) ? 0 : element == null ? 6 : parseInt(element)));
		thisAtom.isotope = atts.get("isotope");
		s = atts.get("charge");
		if (s != null) {
			thisAtom.formalCharge = parseInt(s);
		}

		rdr.handleCoordinates(atts);

		s = atts.get("attachments");
		if (s != null) {
			thisNode.setMultipleAttachments(split(s.trim(), " "));
		}

		s = atts.get("bondordering");
		if (s != null) {
			thisNode.setBondOrdering(split(s.trim(), " "));
		}
		return thisNode;
	}

  private boolean checkWarningOK(String warning) {
    return (warning == null || warning.indexOf("valence") >= 0
        || warning.indexOf("very close") >= 0
        || warning.indexOf("two identical colinear bonds") >= 0);
  }
  
  private CDNode setFragment(String id, Map<String, String> atts) {
	    fragments.push(thisFragmentID = id);
	    CDNode fragmentNode = (thisNode == null || !thisNode.isFragment ? null : thisNode);
	    if (fragmentNode != null) {
	      fragmentNode.setInnerFragmentID(id);
	    }
	    String s = atts.get("connectionorder");
	    if (s != null) {
	      thisNode.setConnectionOrder(s.trim().split(" "));
	    }
	    return fragmentNode;
	  }

  /**
   * Process the bond tags. We only look at the following attributes:
   * 
   * B beginning atom (atom1)
   * 
   * E ending atom (atom2)
   * 
   * BeginAttach associates atom1 with a fragment
   * 
   * EndAttach associates atom2 with a fragment
   * 
   * Order -- the bond order
   * 
   * Display -- wedges and such
   * 
   * Display2 -- only important here for partial bonds
   * 
   * bonds to multiple attachments are not actually made.
   * 
   * @param id
   * @param atts 
   * @return the bond
   * 
   */
  private CDBond setBond(String id, Map<String, String> atts) {
    String atom1 = atts.get("b");
    String atom2 = atts.get("e");
    String a = atts.get("beginattach");
    int beginAttach = (a == null ? 0 : parseInt(a));
    a = atts.get("endattach");
    int endAttach = (a == null ? 0 : parseInt(a));
    String s = atts.get("order");
    String disp = atts.get("display");
    String disp2 = atts.get("display2");
    int order = rdr.getBondOrder("null");
    boolean invertEnds = false;
    if (disp == null) {
      if (s == null) {
        order = 1;
      } else if (s.equals("1.5")) {
        order = rdr.getBondOrder("delocalized");
      } else {
        if (s.indexOf(".") > 0 && !"Dash".equals(disp2)) {
          // partial only works with "dash" setting for second line
          s = s.substring(0, s.indexOf("."));
        }
        order = rdr.getBondOrder(s);
      }
    } else if (disp.equals("WedgeBegin")) {
      order = rdr.getBondOrder("up"); // near
    } else if (disp.equals("Hash") || disp.equals("WedgedHashBegin")) {
      order = rdr.getBondOrder("down");
    } else if (disp.equals("WedgeEnd")) {
      invertEnds = true;
      order = rdr.getBondOrder("up");
    } else if (disp.equals("WedgedHashEnd")) {
      invertEnds = true;
      order = rdr.getBondOrder("down");
    } else if (disp.equals("Bold")) {
      order = rdr.getBondOrder("single");
    } else if (disp.equals("Wavy")) {
      order = rdr.getBondOrder("either");
    }
    if (order == rdr.getBondOrder("null")) {
      // dative, ionic, hydrogen, threecenter
      System.err.println("XmlChemDrawReader ignoring bond type " + s);
      return null;
    }
    CDBond b = (invertEnds ? new CDBond(id, atom2, atom1, order)
        : new CDBond(id, atom1, atom2, order));

    CDNode node1 = getAtom(b.atomIndex1);
    CDNode node2 = getAtom(b.atomIndex2);

    if (order == rdr.getBondOrder("either")) {
      if (!nostereo.contains(node1))
        nostereo.add(node1);
      if (!nostereo.contains(node2))
        nostereo.add(node2);
    }

    if (node1.hasMultipleAttachments) {
      node1.attachedAtom = node2;
      return b;
    } else if (node2.hasMultipleAttachments) {
      node2.attachedAtom = node1;
      return b;
    }

    if (node1.isFragment && beginAttach == 0)
      beginAttach = 1;
    if (node2.isFragment && endAttach == 0)
      endAttach = 1;
    if (beginAttach > 0) {
      (invertEnds ? node2 : node1).addAttachedAtom(b, beginAttach);
    }
    if (endAttach > 0) {
      (invertEnds ? node1 : node2).addAttachedAtom(b, endAttach);
    }
    if (node1.isExternalPt) {
      node1.setInternalAtom(node2);
    }
    if (node2.isExternalPt) {
      node2.setInternalAtom(node1);
    }

    addBond(b);

    return b;
  }

  private void setBracketedGroup(String id, Map<String, String> atts) {
	    String usage = atts.get("bracketusage");
	    if (bracketedGroups == null)
	      bracketedGroups = new Stack<>();
	    if ("MultipleGroup".equals(usage)) {
	      String[] ids = getTokens(atts.get("bracketedobjectids"));
	      int repeatCount = parseInt(atts.get("repeatcount"));
	      bracketedGroups.add(new BracketedGroup(id, ids, repeatCount));
		  haveConnectedFragments = true;
	    }
	  }

  /**
   * Set the 2D or pseudo-3D coordinates of the atoms. ChemDraw pseudo-3D is
   * just a z-layering of chunks of the molecule. Nothing really useful. These
   * coordinates are ignored if there are any atoms also with 2D coordinates or
   * for FILTER "NO3D". So, pretty much, the z coordinates are never used.
   * 
   * @param key
   * @param atts 
   */
  public void setAtom(String key, Map<String, String> atts) {
    String xyz = atts.get(key);
    String[] tokens = getTokens(xyz);
    double x = parseDouble(tokens[0]);
    double y = -parseDouble(tokens[1]);
    double z = (key == "xyz" ? parseDouble(tokens[2]) : 0);
    if (x < minX)
      minX = x;
    if (x > maxX)
      maxX = x;
    if (y < minY)
      minY = y;
    if (y > maxY)
      maxY = y;
    if (z < minZ)
      minZ = z;
    if (z > maxZ)
      maxZ = z;
    thisAtom.set(x, y, z);
  }


  /**
   * Remove fragment, external point, or invalid unconnected nodes (including
   * unconnected carbon nodes, which can arise from deletions (in my experience)
   * and are then not noticed because they have no associated text.
   */
  private void fixInvalidAtoms() {
    for (int i = getAtomCount(); --i >= 0;) {
      CDNode a = getAtom(i);
      a.intID = Integer.MIN_VALUE;
      if (a.isFragment || a.isExternalPt || !a.isConnected
          && (!a.isValid || a.elementNumber < 10)) {
        bsAtoms.clear(a.index);
      }
    }
    reserializeAtoms();

    bsBonds.clear();
    for (int i = getBondCount(); --i >= 0;) {
        CDBond b = getBond(i);
        if (b.isValid()) {
          bsBonds.set(i);
        } else {
          // fragment or invalidated
        }
      }
  }

  private void reserializeAtoms() {
    for (int p = 0, i = bsAtoms.nextSetBit(0); i >= 0; i = bsAtoms.nextSetBit(i + 1)) {
      getAtom(i).intID = ++p;
    }
  }

  private void reindexAtomsAndBonds() {
    reserializeAtoms();
    for (int p = 0, i = bsAtoms.nextSetBit(0); i >= 0; i = bsAtoms.nextSetBit(i + 1)) {
      getAtom(i).index = p++;
    }
    for (int i = bsBonds.nextSetBit(0); i >= 0; i = bsBonds.nextSetBit(i + 1)) {
      CDBond b = getBond(i);
      b.atomIndex1 = b.atom1.index;
      b.atomIndex2 = b.atom2.index;
    }
  }
  /**
   * First fix all the attachments, tying together the atoms identified as
   * ExternalConnectionPoints with atoms of bonds indicating "BeginAttach" or
   * "EndAttach".
   * 
   * Then flag all unconnected atoms and also remove any wedges or hashes that
   * are associated with bonds to atoms that also have wavy bonds.
   */
  private void fixConnections() {

    // fix attachments for fragments

    for (int i = getAtomCount(); --i >= 0;) {
      CDNode a = getAtom(i);
      if (a.isFragment || a.hasMultipleAttachments)
        a.fixAttachments();
    }

    // indicate all atoms that are connected

    for (int i = 0, n = getBondCount(); i < n; i++) {
      CDBond b = getBond(i);
      if (b == null) {
        continue; // bond to nickname
      }
      CDNode a1 = b.atom1;
      CDNode a2 = b.atom2;
      
      a1.isConnected = true;
      a2.isConnected = true;
      if (nostereo.contains(a1) != nostereo.contains(a2)) {
        // wavy line, so no stereo bonds here
        b.order = 1;
      }
    }
  }

  private void fixBracketedGroups() {
    if (bracketedGroups == null)
      return;
    //dumpGraph();
    
    for (int i = bracketedGroups.size(); --i >= 0;) {
      bracketedGroups.remove(i).process();
    }
    
    //dumpGraph();

  }

  void dumpGraph() {
    for (int i = 0, n = getAtomCount(); i < n; i++) {
      CDNode a =  getAtom(i);
      System.out.println("CDXMLP " + i
          + " id="  + a.id + a.bsConnections
          + " cid="  + a.clonedIndex
          + " fa=" + a.bsFragmentAtoms
          + " xp=" + a.isExternalPt
          + " ifd=" + a.innerFragmentID
          + " ofd= " + a.outerFragmentID
          
          );
    }
    for (int i = 0, n = getBondCount(); i < n; i++) {
      CDBond b =  getBond(i);
      System.out.println("CDXMLP bond " + i + " " 
          + b.atomIndex1 + " " + b.atomIndex2 + b);
    }
    System.out.println(bondIDMap);
    return;
  }

  /**
   * Adjust the scale to have an average bond length of 1.45 Angstroms. This is
   * just to get the structure in the range of other structures rather than
   * being huge.
   * 
   */
  private void centerAndScale() {
    if (minX > maxX)
      return;
    double sum = 0;
    int n = 0;
    double lenH = 1;
    for (int i = getBondCount(); --i >= 0;) {
      CDBond b = getBond(i);
      CDNode a1 = b.atom1;
      CDNode a2 = b.atom2;
      double d = a1.distance(a2);
      if (a1.elementNumber > 1 && a2.elementNumber > 1) {
        sum += d;
        n++;
      } else {
        lenH = d;
      }
    }
    double f = (sum > 0 ? 1.45d * n / sum : lenH > 0 ? 1 / lenH : 1);
    // in case somehow ChemDraw uses Cartesians.
    if (f > 0.5)
      f = 1;

    double cx = (maxX + minX) / 2;
    double cy = (maxY + minY) / 2;
    double cz = (maxZ + minZ) / 2;
    for (int i = getAtomCount(); --i >= 0;) {
      CDNode a = getAtom(i);
      a.x = (a.x - cx) * f;
      a.y = (a.y - cy) * f;
      a.z = (a.z - cz) * f;
    }
  }

  CDNode getAtom(int i) {
    return atoms.get(i);
  }

  CDNode addAtom(CDNode atom) {
    atoms.add(atom);
    bsAtoms.set(atom.index);
    return atom;
  }
  
  int getAtomCount() {
    return atoms.size();
  }

  CDBond addBond(CDBond b) {
    bsBonds.set(b.index = getBondCount());
    bonds.add(b);
    return b;
  }

  CDBond getBond(int i) {
    return bonds.get(i);
  }
  
  int getBondCount() {
    return bonds.size();
  }

  public void finalizeParsing() {
    fixConnections();
    fixInvalidAtoms();
    fixBracketedGroups();
    centerAndScale();
    reindexAtomsAndBonds();    
  }

	public boolean hasConnectedFragments() {
		return haveConnectedFragments;
	}
}
