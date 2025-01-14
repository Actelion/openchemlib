package com.actelion.research.chem.io;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Stack;

/**
 * Minimal ChemDraw CDX to CDXML converter.
 * 
 * Just the basic necessities for CDXParser.
 * 
 * @author hansonr@stolaf.edu
 */
public class CDX2CDXML {

	private Stack<String> objects = new Stack<String>();
	private StringBuffer sb = new StringBuffer();
	private int sbpt;
	private ByteBuffer buf;
	private final static String[] cdxAttributes = { "HashSpacing", "2.50", "MarginWidth", "1.60", "LineWidth", "0.60",
			"BoldWidth", "2", "BondLength", "14.40", "BondSpacing", "18" };

	public CDX2CDXML() {
		// for JavaScript dynamic loading
	}

	public static String fromCDX(byte[] cdx) throws Exception {
		return new CDX2CDXML().cdxToCdxml(cdx);
	}

	private String cdxToCdxml(byte[] cdx) throws Exception {
		buf = ByteBuffer.wrap(cdx);
		buf.order(ByteOrder.LITTLE_ENDIAN);
		try {
			openDocument(sb);
			appendHeader(sb);
			buf.position(22); // header
			processObject(buf.getShort());
			sb.append("</CDXML>\n");
		} catch (Exception e) {
			System.out.println(sb + "\n" + objects);
			e.printStackTrace();
			return null;
		}
		return sb.toString();
	}

	private static void appendHeader(StringBuffer sb) {
		sb.append("<!DOCTYPE CDXML SYSTEM \"http://www.cambridgesoft.com/xml/cdxml.dtd\" >\n");
		startOpenTag(sb, "CDXML");
		addAttributes(sb, cdxAttributes);
		terminateTag(sb);
	}

	/*
	 * Process objects recursively.
	 * 
	 */
	private void processObject(int type) throws Exception {
		int id = buf.getInt();
		type = type & 0xFFFF;
		boolean terminated = false;
		String name = null;
		// boolean hasAttributes = false;
		// System.out.println("type " + Integer.toHexString(type) + " id " + id);
		switch (type) {
		case kCDXObj_Document: // 0x8000
		case kCDXObj_Group: // 0x8002
		default:
			id = Integer.MIN_VALUE;
			terminated = true;
			break;
		case kCDXObj_Page: // 0x8001
			name = "page";
			id = Integer.MIN_VALUE;
			// hasAttributes = true;
			break;
		case kCDXObj_Fragment: // 0x8003
			name = "fragment";
			// hasAttributes = true;
			break;
		case kCDXObj_Node: // 0x8004
			name = "n";
			// hasAttributes = true;
			break;
		case kCDXObj_Bond: // 0x8005
			name = "b";
			// hasAttributes = true;
			break;
		case kCDXObj_Text: // 0x8006
			name = "t";
			id = Integer.MIN_VALUE;
			// hasAttributes = false;
			break;
		case kCDXObj_BracketedGroup: // 0x8017
			name = "bracketedgroup";
			break;
		case kCDXObj_BracketAttachment: // 0x8018
			name = "bracketattachment";
			break;
		case kCDXObj_CrossingBond: // 0x8019
			name = "crossingbond";
			break;
		}
		sbpt = sb.length();
		objects.push(name);
		if (name != null) {
			startOpenTag(sb, name);
			if (id != Integer.MIN_VALUE) {
				addAttribute(sb, "id", "" + id);
			}
		}
		int prop;
		while ((prop = buf.getShort()) != 0) {
			if ((prop & 0x8000) != 0) {
				if (!terminated) {
					terminateTag(sb);
					terminated = true;
				}
				processObject(prop);
				continue;
			}
			int len = readLength();
			switch (type) {
			case kCDXObj_Node:
				writeNodeProperties(prop, len);
				break;
			case kCDXObj_Text:
				if (!terminated) {
					terminateTag(sb);
					terminated = true;
				}
				writeTextProperty(prop, len);
				break;
			case kCDXObj_Bond:
				writeBondProperties(prop, len);
				break;
			case kCDXObj_BracketedGroup:
				writeBracketedGroupProperties(prop, len);
				break;
			case kCDXObj_CrossingBond:
				writeCrossingBondProperties(prop, len);
				break;
			default:
				skip(len);
				break;
			}
		}
		if (name != null) {
			if (!terminated) {
				terminateEmptyTag(sb);
			} else {
				closeTag(sb, name);
			}
		}
	}

	private void writeBracketedGroupProperties(int prop, int len) throws Exception {
		switch (prop) {
		case kCDXProp_BracketedObjects:
			addAttribute(sb, "BracketedObjectIDs", readArray(len));
			break;
		case kCDXProp_Bracket_RepeatCount:
			addAttribute(sb, "RepeatCount", "" + (int) readFloat64());
			break;
		case kCDXProp_Bracket_Usage:
			int usage = readInt(len);
			String sval = null;
			switch (usage) {
			case 16:
				sval = "MultipleGroup";
				break;
			default:
				removeObject();
				return;
			}
			addAttribute(sb, "BracketUsage", sval);
			break;
		default:
			skip(len);
		}

	}

	private void writeCrossingBondProperties(int prop, int len) throws Exception {
		switch (prop) {
		case kCDXProp_Bracket_BondID:
			addAttribute(sb, "BondID", "" + readInt(len));
			break;
		case kCDXProp_Bracket_InnerAtomID:
			addAttribute(sb, "InnerAtomID", "" + readInt(len));
			break;
		default:
			skip(len);
		}
	}

	private void writeNodeProperties(int prop, int len) throws Exception {
		switch (prop) {
		case kCDXProp_2DPosition:
			double y = toPoint(readInt(len >> 1));
			double x = toPoint(readInt(len >> 1));
			addAttribute(sb, "p", x + " " + y);
			break;
		case kCDXProp_Node_Type:
			String nodeType = getNodeType(readInt(len));
			// Jmol will put "_" here for nodes it does not need
			addAttribute(sb, "NodeType", nodeType);
			break;
		case kCDXProp_Node_Element: // 0x0402 The atomic number of the atom representing this node. (INT16)
			addAttribute(sb, "Element", "" + readInt(len));
			break;
		case kCDXProp_Atom_Isotope: // 0x0420 The absolute isotopic mass of an atom (2 for deuterium 14 for
									// carbon-14). (INT16)
			addAttribute(sb, "Isotope", "" + readInt(len));
			break;
		case kCDXProp_Atom_Charge: // 0x0421 The atomic charge of an atom. (INT8)
			addAttribute(sb, "Charge", "" + readInt(len));
			break;
		case kCDXProp_ChemicalWarning:
			addAttribute(sb, "Warning", readString(len));
			break;
		case kCDXProp_Atom_BondOrdering: // 0x0431 An ordering of the bonds to this node used for stereocenters
											// fragments and named alternative groups with more than one attachment.
											// (CDXObjectIDArray)
			addAttribute(sb, "BondOrdering", readArray(len));
			break;
		case kCDXProp_Frag_ConnectionOrder: // 0x0505 An ordered list of attachment points within a fragment.
											// (CDXObjectIDArray)
			addAttribute(sb, "ConnectionOrder", readArray(len));
			break;
		case kCDXProp_Node_Attachments: // 0x0432 For multicenter attachment nodes or variable attachment nodes a list
										// of IDs of the nodes which are multiply or variably attached to this node.
										// (CDXObjectIDArrayWithCounts)
			addAttribute(sb, "Attachments", readArray(-1));
			break;
		case kCDXProp_Atom_GenericNickname: // 0x0433 The name of the generic nickname. (CDXString)
			addAttribute(sb, "GenericNickname", readString(len));
			break;
		default:
			skip(len);
		}
	}

	private void writeBondProperties(int prop, int len) throws Exception {
		switch (prop) {
		case kCDXProp_Bond_Order:
			String order = getBondOrder(readInt(len));
			if (order == null) {
				removeObject();
				return;
			}
			addAttribute(sb, "Order", order);
			break;
		case kCDXProp_Bond_Display:
			String d = getBondDisplay(readInt(len));
			if (d == null) {
				removeObject();
				return;
			}
			addAttribute(sb, "Display", d);
			break;
		case kCDXProp_Bond_Display2:
			String d2 = getBondDisplay(readInt(len));
			if (d2 != null)
				addAttribute(sb, "Display2", d2);
			break;
		case kCDXProp_Bond_Begin:
			addAttribute(sb, "B", "" + readInt(len));
			break;
		case kCDXProp_Bond_End:
			addAttribute(sb, "E", "" + readInt(len));
			break;
		case kCDXProp_Bond_BeginAttach:
			addAttribute(sb, "BeginAttach", "" + readInt(len));
			break;
		case kCDXProp_Bond_EndAttach:
			addAttribute(sb, "EndAttach", "" + readInt(len));
			break;
		default:
			skip(len);
		}
	}

	private void writeTextProperty(int prop, int len) throws Exception {
		switch (prop) {
		case kCDXProp_Text:
			String text = readString(len);
			System.out.println("CDXMLW text=" + text);
			openTag(sb, "s");
			// remove new line char
			sb.setLength(sb.length() - 1);
			sb.append(wrapCData(text));
			closeTag(sb, "s");
			break;
		default:
			skip(len);
		}
	}

	/**
	 * wrap the string as character data, with replacements for [ noted as a list
	 * starting with * after the CDATA termination
	 * 
	 * @param s
	 * @return wrapped text
	 */
	public static String wrapCData(String s) {
		return (s.indexOf("&") < 0 && s.indexOf("<") < 0 ? s
				: "<![CDATA[" + s.replace("]]>", "]]]]><![CDATA[>") + "]]>");
	}

	private static String getNodeType(int n) {
		String name = null;
		switch (n) {
		case kCDXNodeType_Unspecified:
			return "Unspecified";
		case kCDXNodeType_Element:
			return "Element";
		case kCDXNodeType_Nickname:
			return "Nickname";
		case kCDXNodeType_Fragment:
			return "Fragment";
		case kCDXNodeType_GenericNickname:
			return "GenericNickname";
		case kCDXNodeType_MultiAttachment:
			return "MultiAttachment";
		case kCDXNodeType_VariableAttachment:
			return "VariableAttachment";
		case kCDXNodeType_ExternalConnectionPoint:
			return "ExternalConnectionPoint";
		case kCDXNodeType_ElementList:
			name = "ElementList";
			break;
		case kCDXNodeType_ElementListNickname:
			name = "ElementListNickname";
			break;
		case kCDXNodeType_Formula:
			name = "Formula";
			break;
		case kCDXNodeType_AnonymousAlternativeGroup:
			name = "AnonymousAlternativeGroup";
			break;
		case kCDXNodeType_NamedAlternativeGroup:
			name = "NamedAnonymousGroup";
			break;
		case kCDXNodeType_LinkNode:
			name = "LinkNode";
			break;
		}
		System.err.println("CDXMLWriter Node type " + name + " not identified");
		return "_";
	}

	private static String getBondDisplay(int i) {
		switch (i) {
		case kCDXBondDisplay_Solid:
			return "Solid";
		case kCDXBondDisplay_Dash:
			return "Dash";
		case kCDXBondDisplay_Hash:
			return "Hash";
		case kCDXBondDisplay_WedgedHashBegin:
			return "WedgedHashBegin";
		case kCDXBondDisplay_WedgedHashEnd:
			return "WedgedHashEnd";
		case kCDXBondDisplay_Bold:
			return "Bold";
		case kCDXBondDisplay_WedgeBegin:
			return "WedgeBegin";
		case kCDXBondDisplay_WedgeEnd:
			return "WedgeEnd";
		case kCDXBondDisplay_Wavy:
			return "Wavy";
		case kCDXBondDisplay_HollowWedgeBegin:
			return "HollowWedgeBegin";
		case kCDXBondDisplay_HollowWedgeEnd:
			return "HollowWedgeEnd";
		case kCDXBondDisplay_WavyWedgeBegin:
			return "WavyWedgeBegin";
		case kCDXBondDisplay_WavyWedgeEnd:
			return "WavyWedgeEnd";
		case kCDXBondDisplay_Dot:
			return "Dot";
		case kCDXBondDisplay_DashDot:
			return "DashDot";
		}
		return null;
	}

	private static String getBondOrder(int i) {
		switch (i) {
		case kCDXBondOrder_Single:
			return "1";
		case kCDXBondOrder_Double:
			return "2";
		case kCDXBondOrder_Triple:
			return "3";
		case kCDXBondOrder_Quadruple:
			return "4";
		case kCDXBondOrder_Quintuple:
			return "5";
		case kCDXBondOrder_Sextuple:
			return "6";
		case kCDXBondOrder_Half:
			return "0.5";
		case kCDXBondOrder_OneHalf:
			return "1.5";
		case kCDXBondOrder_TwoHalf:
			return "2.5";
		case kCDXBondOrder_ThreeHalf:
			return "3.5";
		case kCDXBondOrder_FourHalf:
			return "4.5";
		case kCDXBondOrder_FiveHalf:
			return "5.5";
		case kCDXBondOrder_Dative:
			return "dative";
		case kCDXBondOrder_Ionic:
			return "ionic";
		case kCDXBondOrder_Hydrogen:
			return "hydrogen";
		case kCDXBondOrder_ThreeCenter:
			return "threecenter";
		}
		return null;
	}

	private void removeObject() {
		sb.setLength(sbpt);
	}

	private void skip(int len) {
		buf.position(buf.position() + len);
	}

	private double readFloat64() throws Exception {
		return buf.getDouble();
	}

	private int readInt(int len) throws Exception {
		switch (len) {
		case 1:
			return (256 + buf.get()) % 256;
		case 2:
			return buf.getShort();
		case 4:
			return buf.getInt();
		case 8:
			return (int) buf.getLong();
		}
		System.err.println("CDXMLWriter.readInt len " + len);
		return 0;
	}

	private String readString(int len) throws Exception {
		int nStyles = buf.getShort();
		len -= 2;
		switch (nStyles) {
		case 0:
			break;
		default:
			skip(nStyles * 10);
			len -= nStyles * 10;
			break;
		}
	    byte[] temp = new byte[len];
	    buf.get(temp, 0, len);
	    return new String(temp, 0, len, "UTF-8");
	}

	private String readArray(int len) throws Exception {
		int n = (len < 0 ? buf.getShort() : len / 4);
		String s = "";
		for (int i = n; --i >= 0;) {
			s += " " + buf.getInt();
		}
		return s.trim();
	}

	private int readLength() throws Exception {
		int len = buf.getShort();
		if (len == -1) {
			// large 4-byte length
			len = buf.getInt();
		}
		return len;
	}

	private static double toPoint(int i) {
		return Math.round(i / 655.36) / 100.0;
	}

	public static void main(String[] args) {
//    try {
//      FileInputStream fis = new FileInputStream("c:/temp/t2.cdx");
//      BinaryDocument doc = new BinaryDocument();
//      doc.setStream(new BufferedInputStream(fis), false);
//      System.out.println(fromCDX(doc));
//    } catch (Exception e) {
//      e.printStackTrace();
//    }
	}

	// from
	// https://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/CDXConstants.h
	// via Excel

	// General properties.
	// private final static int kCDXProp_EndObject = 0x0000; // 0x0000 Marks end of
	// object.
	// private final static int kCDXProp_CreationUserName = 0x0001; // 0x0001 The
	// name of the creator (program user's name) of the document. (CDXString)
	// private final static int kCDXProp_CreationDate = 0x0002; // 0x0002 The time
	// of object creation. (CDXDate)
	// private final static int kCDXProp_CreationProgram = 0x0003; // 0x0003 The
	// name of the program including version and platform that created the
	// associated CDX object. ChemDraw 4.0 uses ChemDraw 4.0 as the value of
	// CreationProgram. (CDXString)
	// private final static int kCDXProp_ModificationUserName = 0x0004; // 0x0004
	// The name of the last modifier (program user's name) of the document.
	// (CDXString)
	// private final static int kCDXProp_ModificationDate = 0x0005; // 0x0005 Time
	// of the last modification. (CDXDate)
	// private final static int kCDXProp_ModificationProgram = 0x0006; // 0x0006 The
	// name of the program including version and platform of the last program to
	// perform a modification. ChemDraw 4.0 uses ChemDraw 4.0 as the value of
	// CreationProgram. (CDXString)
	// private final static int kCDXProp_Unused1 = 0x0007; // 0x0007 Table of
	// contents. (obsolete)
	// private final static int kCDXProp_Name = 0x0008; // 0x0008 Name of an object.
	// (CDXString)
	// private final static int kCDXProp_Comment = 0x0009; // 0x0009 An arbitrary
	// string intended to be meaningful to a user. (CDXString)
	// private final static int kCDXProp_ZOrder = 0x000A; // 0x000A Back-to-front
	// ordering index in 2D drawing. (INT16)
	// private final static int kCDXProp_RegistryNumber = 0x000B; // 0x000B A
	// registry or catalog number of a molecule object. (CDXString)
	// private final static int kCDXProp_RegistryAuthority = 0x000C; // 0x000C A
	// string that specifies the authority which issued a registry or catalog
	// number. Some examples of registry authorities are CAS Beilstein Aldrich and
	// Merck. (CDXString)
	// private final static int kCDXProp_Unused2 = 0x000D; // 0x000D Indicates that
	// this object (the reference object) is an alias to an object elsewhere in the
	// document (the target object). The attributes and contained objects should be
	// taken from the target object. (obsolete)
	// private final static int kCDXProp_RepresentsProperty = 0x000E; // 0x000E
	// Indicates that this object represents some property in some other object.
	// (CDXRepresentsProperty)
	// private final static int kCDXProp_IgnoreWarnings = 0x000F; // 0x000F
	// Signifies whether chemical warnings should be suppressed on this object.
	// (CDXBooleanImplied)
	private final static int kCDXProp_ChemicalWarning = 0x0010; // 0x0010 A warning concerning possible chemical
																// problems with this object. (CDXString)
	// private final static int kCDXProp_Visible = 0x0011; // 0x0011 The object is
	// visible if non-zero. (CDXBoolean)

	// Fonts.
	// private final static int kCDXProp_FontTable = 0x0100; // 0x0100 A list of
	// fonts used in the document. (CDXFontTable)

	// Coordinates.
	private final static int kCDXProp_2DPosition = 0x0200; // 0x0200 The 2D location (in the order of vertical and
															// horizontal locations) of an object. (CDXPoint2D)
	// private final static int kCDXProp_3DPosition = 0x0201; // 0x0201 The 3D
	// location (in the order of X- Y- and Z-locations in right-handed coordinate
	// system) of an object in CDX coordinate units. The precise meaning of this
	// attribute varies depending on the type of object. (CDXPoint3D)
	// private final static int kCDXProp_2DExtent = 0x0202; // 0x0202 The width and
	// height of an object in CDX coordinate units. The precise meaning of this
	// attribute varies depending on the type of object. (CDXPoint2D)
	// private final static int kCDXProp_3DExtent = 0x0203; // 0x0203 The width
	// height and depth of an object in CDX coordinate units (right-handed
	// coordinate system). The precise meaning of this attribute varies depending on
	// the type of object. (CDXPoint3D)
	// private final static int kCDXProp_BoundingBox = 0x0204; // 0x0204 The
	// smallest rectangle that encloses the graphical representation of the object.
	// (CDXRectangle)
	// private final static int kCDXProp_RotationAngle = 0x0205; // 0x0205 The
	// angular orientation of an object in degrees * 65536. (INT32)
	// private final static int kCDXProp_BoundsInParent = 0x0206; // 0x0206 The
	// bounds of this object in the coordinate system of its parent (used for pages
	// within tables). (CDXRectangle)
	// private final static int kCDXProp_3DHead = 0x0207; // 0x0207 The 3D location
	// (in the order of X- Y- and Z-locations in right-handed coordinate system) of
	// the head of an object in CDX coordinate units. The precise meaning of this
	// attribute varies depending on the type of object. (CDXPoint3D)
	// private final static int kCDXProp_3DTail = 0x0208; // 0x0208 The 3D location
	// (in the order of X- Y- and Z-locations in right-handed coordinate system) of
	// the tail of an object in CDX coordinate units. The precise meaning of this
	// attribute varies depending on the type of object. (CDXPoint3D)
	// private final static int kCDXProp_TopLeft = 0x0209; // 0x0209 The location of
	// the top-left corner of a quadrilateral object possibly in a rotated or skewed
	// frame. (CDXPoint2D)
	// private final static int kCDXProp_TopRight = 0x020A; // 0x020A The location
	// of the top-right corner of a quadrilateral object possibly in a rotated or
	// skewed frame. (CDXPoint2D)
	// private final static int kCDXProp_BottomRight = 0x020B; // 0x020B The
	// location of the bottom-right corner of a quadrilateral object possibly in a
	// rotated or skewed frame. (CDXPoint2D)
	// private final static int kCDXProp_BottomLeft = 0x020C; // 0x020C The location
	// of the bottom-left corner of a quadrilateral object possibly in a rotated or
	// skewed frame. (CDXPoint2D)

	// Colors.
	// private final static int kCDXProp_ColorTable = 0x0300; // 0x0300 The color
	// palette used throughout the document. (CDXColorTable)
	// private final static int kCDXProp_ForegroundColor = 0x0301; // 0x0301 The
	// foreground color of an object represented as the two-based index into the
	// object's color table. (UINT16)
	// private final static int kCDXProp_BackgroundColor = 0x0302; // 0x0302 The
	// background color of an object represented as the two-based index into the
	// object's color table. (INT16)

	// Atom properties.
	private final static int kCDXProp_Node_Type = 0x0400; // 0x0400 The type of a node object. (INT16)
	// private final static int kCDXProp_Node_LabelDisplay = 0x0401; // 0x0401 The
	// characteristics of node label display. (INT8)
	private final static int kCDXProp_Node_Element = 0x0402; // 0x0402 The atomic number of the atom representing this
																// node. (INT16)
	// private final static int kCDXProp_Atom_ElementList = 0x0403; // 0x0403 A list
	// of atomic numbers. (CDXElementList)
	// private final static int kCDXProp_Atom_Formula = 0x0404; // 0x0404 The
	// composition of a node representing a fragment whose composition is known but
	// whose connectivity is not. For example C<sub>4</sub>H<sub>9</sub> represents
	// a mixture of the 4 butyl isomers. (CDXFormula)
	private final static int kCDXProp_Atom_Isotope = 0x0420; // 0x0420 The absolute isotopic mass of an atom (2 for
																// deuterium 14 for carbon-14). (INT16)
	private final static int kCDXProp_Atom_Charge = 0x0421; // 0x0421 The atomic charge of an atom. (INT8)
	// private final static int kCDXProp_Atom_Radical = 0x0422; // 0x0422 The atomic
	// radical attribute of an atom. (UINT8)
	// private final static int kCDXProp_Atom_RestrictFreeSites = 0x0423; // 0x0423
	// Indicates that up to the specified number of additional substituents are
	// permitted on this atom. (UINT8)
	// private final static int kCDXProp_Atom_RestrictImplicitHydrogens = 0x0424; //
	// 0x0424 Signifies that implicit hydrogens are not allowed on this atom.
	// (CDXBooleanImplied)
	// private final static int kCDXProp_Atom_RestrictRingBondCount = 0x0425; //
	// 0x0425 The number of ring bonds attached to an atom. (INT8)
	// private final static int kCDXProp_Atom_RestrictUnsaturatedBonds = 0x0426; //
	// 0x0426 Indicates whether unsaturation should be present or absent. (INT8)
	// private final static int kCDXProp_Atom_RestrictRxnChange = 0x0427; // 0x0427
	// If present signifies that the reaction change of an atom must be as
	// specified. (CDXBooleanImplied)
	// private final static int kCDXProp_Atom_RestrictRxnStereo = 0x0428; // 0x0428
	// The change of stereochemistry of an atom during a reaction. (INT8)
	// private final static int kCDXProp_Atom_AbnormalValence = 0x0429; // 0x0429
	// Signifies that an abnormal valence for an atom is permitted.
	// (CDXBooleanImplied)
	// private final static int kCDXProp_Unused3 = 0x042A; // 0x042A
	// private final static int kCDXProp_Atom_NumHydrogens = 0x042B; // 0x042B The
	// number of (explicit) hydrogens in a labeled atom consisting of one heavy atom
	// and (optionally) the symbol H (e.g. CH<sub>3</sub>). (UINT16)
	// private final static int kCDXProp_Unused4 = 0x042C; // 0x042C
	// private final static int kCDXProp_Unused5 = 0x042D; // 0x042D
	// private final static int kCDXProp_Atom_HDot = 0x042E; // 0x042E Signifies the
	// presence of an implicit hydrogen with stereochemistry specified equivalent to
	// an explicit H atom with a wedged bond. (CDXBooleanImplied)
	// private final static int kCDXProp_Atom_HDash = 0x042F; // 0x042F Signifies
	// the presence of an implicit hydrogen with stereochemistry specified
	// equivalent to an explicit H atom with a hashed bond. (CDXBooleanImplied)
	// private final static int kCDXProp_Atom_Geometry = 0x0430; // 0x0430 The
	// geometry of the bonds about this atom. (INT8)
	private final static int kCDXProp_Atom_BondOrdering = 0x0431; // 0x0431 An ordering of the bonds to this node used
																	// for stereocenters fragments and named alternative
																	// groups with more than one attachment.
																	// (CDXObjectIDArray)
	private final static int kCDXProp_Node_Attachments = 0x0432; // 0x0432 For multicenter attachment nodes or variable
																	// attachment nodes a list of IDs of the nodes which
																	// are multiply or variably attached to this node.
																	// (CDXObjectIDArrayWithCounts)
	private final static int kCDXProp_Atom_GenericNickname = 0x0433; // 0x0433 The name of the generic nickname.
																		// (CDXString)
	// private final static int kCDXProp_Atom_AltGroupID = 0x0434; // 0x0434 The ID
	// of the alternative group object that describes this node. (CDXObjectID)
	// private final static int kCDXProp_Atom_RestrictSubstituentsUpTo = 0x0435; //
	// 0x0435 Indicates that substitution is restricted to no more than the
	// specified value. (UINT8)
	// private final static int kCDXProp_Atom_RestrictSubstituentsExactly = 0x0436;
	// // 0x0436 Indicates that exactly the specified number of substituents must be
	// present. (UINT8)
	// private final static int kCDXProp_Atom_CIPStereochemistry = 0x0437; // 0x0437
	// The node's absolute stereochemistry according to the Cahn-Ingold-Prelog
	// system. (INT8)
	// private final static int kCDXProp_Atom_Translation = 0x0438; // 0x0438
	// Provides for restrictions on whether a given node may match other more- or
	// less-general nodes. (INT8)
	// private final static int kCDXProp_Atom_AtomNumber = 0x0439; // 0x0439 Atom
	// number as text. (CDXString)
	// private final static int kCDXProp_Atom_ShowQuery = 0x043A; // 0x043A Show the
	// query indicator if non-zero. (CDXBoolean)
	// private final static int kCDXProp_Atom_ShowStereo = 0x043B; // 0x043B Show
	// the stereochemistry indicator if non-zero. (CDXBoolean)
	// private final static int kCDXProp_Atom_ShowAtomNumber = 0x043C; // 0x043C
	// Show the atom number if non-zero. (CDXBoolean)
	// private final static int kCDXProp_Atom_LinkCountLow = 0x043D; // 0x043D Low
	// end of repeat count for link nodes. (FLOAT64)
	// private final static int kCDXProp_Atom_LinkCountHigh = 0x043E; // 0x043E High
	// end of repeat count for link nodes. (INT16)
	// private final static int kCDXProp_Atom_IsotopicAbundance = 0x043F; // 0x043F
	// Isotopic abundance of this atom's isotope. (INT8)
	// private final static int kCDXProp_Atom_ExternalConnectionType = 0x0440; //
	// 0x0440 Type of external connection for atoms of type
	// kCDXNodeType_ExternalConnectionPoint. (INT8)

	// Molecule properties.
	// private final static int kCDXProp_Mole_Racemic = 0x0500; // 0x0500 Indicates
	// that the molecule is a racemic mixture. (CDXBoolean)
	// private final static int kCDXProp_Mole_Absolute = 0x0501; // 0x0501 Indicates
	// that the molecule has known absolute configuration. (CDXBoolean)
	// private final static int kCDXProp_Mole_Relative = 0x0502; // 0x0502 Indicates
	// that the molecule has known relative stereochemistry but unknown absolute
	// configuration. (CDXBoolean)
	// private final static int kCDXProp_Mole_Formula = 0x0503; // 0x0503 The
	// molecular formula representation of a molecule object. (CDXFormula)
	// private final static int kCDXProp_Mole_Weight = 0x0504; // 0x0504 The average
	// molecular weight of a molecule object. (FLOAT64)
	private final static int kCDXProp_Frag_ConnectionOrder = 0x0505; // 0x0505 An ordered list of attachment points
																		// within a fragment. (CDXObjectIDArray)

	// Bond properties.
	private final static int kCDXProp_Bond_Order = 0x0600; // 0x0600 The order of a bond object. (INT16)
	private final static int kCDXProp_Bond_Display = 0x0601; // 0x0601 The display type of a bond object. (INT16)
	private final static int kCDXProp_Bond_Display2 = 0x0602; // 0x0602 The display type for the second line of a double
																// bond. (INT16)
	// private final static int kCDXProp_Bond_DoublePosition = 0x0603; // 0x0603 The
	// position of the second line of a double bond. (INT16)
	private final static int kCDXProp_Bond_Begin = 0x0604; // 0x0604 The ID of the CDX node object at the first end of a
															// bond. (CDXObjectID)
	private final static int kCDXProp_Bond_End = 0x0605; // 0x0605 The ID of the CDX node object at the second end of a
															// bond. (CDXObjectID)
	// private final static int kCDXProp_Bond_RestrictTopology = 0x0606; // 0x0606
	// Indicates the desired topology of a bond in a query. (INT8)
	// private final static int kCDXProp_Bond_RestrictRxnParticipation = 0x0607; //
	// 0x0607 Specifies that a bond is affected by a reaction. (INT8)
	private final static int kCDXProp_Bond_BeginAttach = 0x0608; // 0x0608 Indicates where within the Bond_Begin node a
																	// bond is attached. (UINT8)
	private final static int kCDXProp_Bond_EndAttach = 0x0609; // 0x0609 Indicates where within the Bond_End node a bond
																// is attached. (UINT8)
	// private final static int kCDXProp_Bond_CIPStereochemistry = 0x060A; // 0x060A
	// The bond's absolute stereochemistry according to the Cahn-Ingold-Prelog
	// system. (INT8)
	// private final static int kCDXProp_Bond_BondOrdering = 0x060B; // 0x060B
	// Ordered list of attached bond IDs. (CDXObjectIDArray)
	// private final static int kCDXProp_Bond_ShowQuery = 0x060C; // 0x060C Show the
	// query indicator if non-zero. (CDXBoolean)
	// private final static int kCDXProp_Bond_ShowStereo = 0x060D; // 0x060D Show
	// the stereochemistry indicator if non-zero. (CDXBoolean)
	// private final static int kCDXProp_Bond_CrossingBonds = 0x060E; // 0x060E
	// Unordered list of IDs of bonds that cross this one (either above or below).
	// (CDXObjectIDArray)
	// private final static int kCDXProp_Bond_ShowRxn = 0x060F; // 0x060F Show the
	// reaction-change indicator if non-zero. (CDXBoolean)

	// Text properties.
	private final static int kCDXProp_Text = 0x0700; // 0x0700 The text of a text object. (CDXString)
	// private final static int kCDXProp_Justification = 0x0701; // 0x0701 The
	// horizontal justification of a text object. (INT8)
	// private final static int kCDXProp_LineHeight = 0x0702; // 0x0702 The line
	// height of a text object. (UINT16)
	// private final static int kCDXProp_WordWrapWidth = 0x0703; // 0x0703 The
	// word-wrap width of a text object. (INT16)
	// private final static int kCDXProp_LineStarts = 0x0704; // 0x0704 The number
	// of lines of a text object followed by that many values indicating the
	// zero-based text position of each line start. (INT16ListWithCounts)
	// private final static int kCDXProp_LabelAlignment = 0x0705; // 0x0705 The
	// alignment of the text with respect to the node position. (INT8)
	// private final static int kCDXProp_LabelLineHeight = 0x0706; // 0x0706 Text
	// line height for atom labels (INT16)
	// private final static int kCDXProp_CaptionLineHeight = 0x0707; // 0x0707 Text
	// line height for non-atomlabel text objects (INT16)
	// private final static int kCDXProp_InterpretChemically = 0x0708; // 0x0708
	// Signifies whether to the text label should be interpreted chemically (if
	// possible). (CDXBooleanImplied)

	// Document properties.
	// private final static int kCDXProp_MacPrintInfo = 0x0800; // 0x0800 The 120
	// byte Macintosh TPrint data associated with the CDX document object. Refer to
	// Macintosh Toolbox manual for detailed description. (Unformatted)
	// private final static int kCDXProp_WinPrintInfo = 0x0801; // 0x0801 The
	// Windows DEVMODE structure associated with the CDX document object.
	// (Unformatted)
	// private final static int kCDXProp_PrintMargins = 0x0802; // 0x0802 The outer
	// margins of the Document. (CDXRectangle)
	// private final static int kCDXProp_ChainAngle = 0x0803; // 0x0803 The default
	// chain angle setting in degrees * 65536. (INT32)
	// private final static int kCDXProp_BondSpacing = 0x0804; // 0x0804 The spacing
	// between segments of a multiple bond measured relative to bond length. (INT16)
	// private final static int kCDXProp_BondLength = 0x0805; // 0x0805 The default
	// bond length. (CDXCoordinate)
	// private final static int kCDXProp_BoldWidth = 0x0806; // 0x0806 The default
	// bold bond width. (CDXCoordinate)
	// private final static int kCDXProp_LineWidth = 0x0807; // 0x0807 The default
	// line width. (CDXCoordinate)
	// private final static int kCDXProp_MarginWidth = 0x0808; // 0x0808 The default
	// amount of space surrounding atom labels. (CDXCoordinate)
	// private final static int kCDXProp_HashSpacing = 0x0809; // 0x0809 The default
	// spacing between hashed lines used in wedged hashed bonds. (CDXCoordinate)
	// private final static int kCDXProp_LabelStyle = 0x080A; // 0x080A The default
	// style for atom labels. (CDXFontStyle)
	// private final static int kCDXProp_CaptionStyle = 0x080B; // 0x080B The
	// default style for non-atomlabel text objects. (CDXFontStyle)
	// private final static int kCDXProp_CaptionJustification = 0x080C; // 0x080C
	// The horizontal justification of a caption (non-atomlabel text object) (INT8)
	// private final static int kCDXProp_FractionalWidths = 0x080D; // 0x080D
	// Signifies whether to use fractional width information when drawing text.
	// (CDXBooleanImplied)
	// private final static int kCDXProp_Magnification = 0x080E; // 0x080E The view
	// magnification factor (INT16)
	// private final static int kCDXProp_WidthPages = 0x080F; // 0x080F The width of
	// the document in pages. (INT16)
	// private final static int kCDXProp_HeightPages = 0x0810; // 0x0810 The height
	// of the document in pages. (INT16)
	// private final static int kCDXProp_DrawingSpaceType = 0x0811; // 0x0811 The
	// type of drawing space used for this document. (INT8)
	// private final static int kCDXProp_Width = 0x0812; // 0x0812 The width of an
	// object in CDX coordinate units possibly in a rotated or skewed frame.
	// (CDXCoordinate)
	// private final static int kCDXProp_Height = 0x0813; // 0x0813 The height of an
	// object in CDX coordinate units possibly in a rotated or skewed frame.
	// (CDXCoordinate)
	// private final static int kCDXProp_PageOverlap = 0x0814; // 0x0814 The amount
	// of overlap of pages when a poster is tiled. (CDXCoordinate)
	// private final static int kCDXProp_Header = 0x0815; // 0x0815 The text of the
	// header. (CDXString)
	// private final static int kCDXProp_HeaderPosition = 0x0816; // 0x0816 The
	// vertical offset of the header baseline from the top of the page.
	// (CDXCoordinate)
	// private final static int kCDXProp_Footer = 0x0817; // 0x0817 The text of the
	// footer. (CDXString)
	// private final static int kCDXProp_FooterPosition = 0x0818; // 0x0818 The
	// vertical offset of the footer baseline from the bottom of the page.
	// (CDXCoordinate)
	// private final static int kCDXProp_PrintTrimMarks = 0x0819; // 0x0819 If
	// present trim marks are to printed in the margins. (CDXBooleanImplied)
	// private final static int kCDXProp_LabelStyleFont = 0x081A; // 0x081A The
	// default font family for atom labels. (INT16)
	// private final static int kCDXProp_CaptionStyleFont = 0x081B; // 0x081B The
	// default font style for captions (non-atom-label text objects). (INT16)
	// private final static int kCDXProp_LabelStyleSize = 0x081C; // 0x081C The
	// default font size for atom labels. (INT16)
	// private final static int kCDXProp_CaptionStyleSize = 0x081D; // 0x081D The
	// default font size for captions (non-atom-label text objects). (INT16)
	// private final static int kCDXProp_LabelStyleFace = 0x081E; // 0x081E The
	// default font style for atom labels. (INT16)
	// private final static int kCDXProp_CaptionStyleFace = 0x081F; // 0x081F The
	// default font face for captions (non-atom-label text objects). (INT16)
	// private final static int kCDXProp_LabelStyleColor = 0x0820; // 0x0820 The
	// default color for atom labels (INT16)
	// private final static int kCDXProp_CaptionStyleColor = 0x0821; // 0x0821 The
	// default color for captions (non-atom-label text objects). (INT16)
	// private final static int kCDXProp_BondSpacingAbs = 0x0822; // 0x0822 The
	// absolute distance between segments of a multiple bond. (CDXCoordinate)
	// private final static int kCDXProp_LabelJustification = 0x0823; // 0x0823 The
	// default justification for atom labels. (INT8)
	// private final static int kCDXProp_FixInplaceExtent = 0x0824; // 0x0824
	// Defines a size for OLE In-Place editing. (CDXPoint2D)
	// private final static int kCDXProp_Side = 0x0825; // 0x0825 A specific side of
	// an object (rectangle). (INT16)
	// private final static int kCDXProp_FixInplaceGap = 0x0826; // 0x0826 Defines a
	// padding for OLE In-Place editing. (CDXPoint2D)

	// Graphic object properties.
	// private final static int kCDXProp_Graphic_Type = 0x0A00; // 0x0A00 The type
	// of graphical object. (INT16)
	// private final static int kCDXProp_Line_Type = 0x0A01; // 0x0A01 The type of a
	// line object. (INT16)
	// private final static int kCDXProp_Arrow_Type = 0x0A02; // 0x0A02 The type of
	// arrow object which represents line arrow arc rectangle or orbital. (INT16)
	// private final static int kCDXProp_Rectangle_Type = 0x0A03; // 0x0A03 The type
	// of a rectangle object. (INT16)
	// private final static int kCDXProp_Oval_Type = 0x0A04; // 0x0A04 The type of
	// an arrow object that represents a circle or ellipse. (INT16)
	// private final static int kCDXProp_Orbital_Type = 0x0A05; // 0x0A05 The type
	// of orbital object. (INT16)
	// private final static int kCDXProp_Bracket_Type = 0x0A06; // 0x0A06 The type
	// of symbol object. (INT16)
	// private final static int kCDXProp_Symbol_Type = 0x0A07; // 0x0A07 The type of
	// symbol object. (INT16)
	// private final static int kCDXProp_Curve_Type = 0x0A08; // 0x0A08 The type of
	// curve object. (INT16)
	// private final static int kCDXProp_Arrow_HeadSize = 0x0A20; // 0x0A20 The size
	// of the arrow's head. (INT16)
	// private final static int kCDXProp_Arc_AngularSize = 0x0A21; // 0x0A21 The
	// size of an arc (in degrees * 10 so 90 degrees = 900). (INT16)
	// private final static int kCDXProp_Bracket_LipSize = 0x0A22; // 0x0A22 The
	// size of a bracket. (INT16)
	// private final static int kCDXProp_Curve_Points = 0x0A23; // 0x0A23 The
	// B&eacute;zier curve's control point locations. (CDXCurvePoints)
	private final static int kCDXProp_Bracket_Usage = 0x0A24; // 0x0A24 The syntactical chemical meaning of the bracket
																// (SRU mer mon xlink etc). (INT8)
	// private final static int kCDXProp_Polymer_RepeatPattern = 0x0A25; // 0x0A25
	// The head-to-tail connectivity of objects contained within the bracket. (INT8)
	// private final static int kCDXProp_Polymer_FlipType = 0x0A26; // 0x0A26 The
	// flip state of objects contained within the bracket. (INT8)
	private final static int kCDXProp_BracketedObjects = 0x0A27; // 0x0A27 The set of objects contained in a
																	// BracketedGroup. (CDXObjectIDArray)
	private final static int kCDXProp_Bracket_RepeatCount = 0x0A28; // 0x0A28 The number of times a multiple-group
																	// BracketedGroup is repeated. (INT16)
	// private final static int kCDXProp_Bracket_ComponentOrder = 0x0A29; // 0x0A29
	// The component order associated with a BracketedGroup. (INT16)
	// private final static int kCDXProp_Bracket_SRULabel = 0x0A2A; // 0x0A2A The
	// label associated with a BracketedGroup that represents an SRU. (CDXString)
	// private final static int kCDXProp_Bracket_GraphicID = 0x0A2B; // 0x0A2B The
	// ID of a graphical object (bracket brace or parenthesis) associated with a
	// Bracket Attachment. (CDXObjectID)
	private final static int kCDXProp_Bracket_BondID = 0x0A2C; // 0x0A2C The ID of a bond that crosses a Bracket
																// Attachment. (CDXObjectID)
	private final static int kCDXProp_Bracket_InnerAtomID = 0x0A2D; // 0x0A2D The ID of the node located within the
																	// Bracketed Group and attached to a bond that
																	// crosses a Bracket Attachment. (CDXObjectID)
	// private final static int kCDXProp_Curve_Points3D = 0x0A2E; // 0x0A2E The
	// B&eacute;zier curve's control point locations. (CDXCurvePoints3D)

	// Objects. = ;
	private final static int kCDXObj_Document = 0x8000; // 0x8000
	private final static int kCDXObj_Page = 0x8001; // 0x8001
	private final static int kCDXObj_Group = 0x8002; // 0x8002
	private final static int kCDXObj_Fragment = 0x8003; // 0x8003
	private final static int kCDXObj_Node = 0x8004; // 0x8004
	private final static int kCDXObj_Bond = 0x8005; // 0x8005
	private final static int kCDXObj_Text = 0x8006; // 0x8006
	// private final static int kCDXObj_Graphic = 0x8007; // 0x8007
	// private final static int kCDXObj_Curve = 0x8008; // 0x8008
	// private final static int kCDXObj_EmbeddedObject = 0x8009; // 0x8009
	// private final static int kCDXObj_NamedAlternativeGroup = 0x800a; // 0x800a
	// private final static int kCDXObj_TemplateGrid = 0x800b; // 0x800b
	// private final static int kCDXObj_RegistryNumber = 0x800c; // 0x800c
	// private final static int kCDXObj_ReactionScheme = 0x800d; // 0x800d
	// private final static int kCDXObj_ReactionStep = 0x800e; // 0x800e
	// private final static int kCDXObj_ObjectDefinition = 0x800f; // 0x800f
	// private final static int kCDXObj_Spectrum = 0x8010; // 0x8010
	// private final static int kCDXObj_ObjectTag = 0x8011; // 0x8011
	// private final static int kCDXObj_OleClientItem = 0x8012; // 0x8012
	// private final static int kCDXObj_Sequence = 0x8013; // 0x8013
	// private final static int kCDXObj_CrossReference = 0x8014; // 0x8014
	// private final static int kCDXObj_Splitter = 0x8015; // 0x8015
	// private final static int kCDXObj_Table = 0x8016; // 0x8016
	private final static int kCDXObj_BracketedGroup = 0x8017; // 0x8017
	private final static int kCDXObj_BracketAttachment = 0x8018; // 0x8018
	private final static int kCDXObj_CrossingBond = 0x8019; // 0x8019
	// private final static int kCDXObj_Border = 0x8020; // 0x8020
	// private final static int kCDXObj_Geometry = 0x8021; // 0x8021
	// private final static int kCDXObj_Constraint = 0x8022; // 0x8022
	// private final static int kCDXObj_TLCPlate = 0x8023; // 0x8023
	// private final static int kCDXObj_TLCLane = 0x8024; // 0x8024
	// private final static int kCDXObj_TLCSpot = 0x8025; // 0x8025
	// Add new objects here
	// private final static int kCDXObj_UnknownObject = 0x8FFF; // 0x8FFF

	private final static int kCDXNodeType_Unspecified = 0;
	private final static int kCDXNodeType_Element = 1;
	private final static int kCDXNodeType_ElementList = 2;
	private final static int kCDXNodeType_ElementListNickname = 3;
	private final static int kCDXNodeType_Nickname = 4;
	private final static int kCDXNodeType_Fragment = 5;
	private final static int kCDXNodeType_Formula = 6;
	private final static int kCDXNodeType_GenericNickname = 7;
	private final static int kCDXNodeType_AnonymousAlternativeGroup = 8;
	private final static int kCDXNodeType_NamedAlternativeGroup = 9;
	private final static int kCDXNodeType_MultiAttachment = 10;
	private final static int kCDXNodeType_VariableAttachment = 11;
	private final static int kCDXNodeType_ExternalConnectionPoint = 12;
	private final static int kCDXNodeType_LinkNode = 13;

	// private final static int kCDXLabelDisplay_Auto = 1;
	// private final static int kCDXLabelDisplay_Left = 2;
	// private final static int kCDXLabelDisplay_Center = 3;
	// private final static int kCDXLabelDisplay_Right = 4;
	// private final static int kCDXLabelDisplay_Above = 5;
	// private final static int kCDXLabelDisplay_Below = 6;
	// private final static int kCDXLabelDisplay_BestInitial = 7;

	// Same as MDL codes
	// private final static int kCDXRadical_None = 0;
	// private final static int kCDXRadical_Singlet = 1;
	// private final static int kCDXRadical_Doublet = 2;
	// private final static int kCDXRadical_Triplet = 3;

	// private final static int kCDXIsotope_Natural = 0;

	// private final static int kCDXRingBondCount_Unspecified = -1;
	// private final static int kCDXRingBondCount_NoRingBonds = 0;
	// private final static int kCDXRingBondCount_AsDrawn = 1;
	// private final static int kCDXRingBondCount_SimpleRing = 2;
	// private final static int kCDXRingBondCount_Fusion = 3;
	// private final static int kCDXRingBondCount_SpiroOrHigher = 4;

	// private final static int kCDXUnsaturation_Unspecified = 0;
	// private final static int kCDXUnsaturation_MustBeAbsent = 1;
	// private final static int kCDXUnsaturation_MustBePresent = 2;

	// private final static int kCDXReactionStereo_Unspecified = 0;
	// private final static int kCDXReactionStereo_Inversion = 1;
	// private final static int kCDXReactionStereo_Retention = 2;

	// private final static int kCDXTranslation_Equal = 0;
	// private final static int kCDXTranslation_Broad = 1;
	// private final static int kCDXTranslation_Narrow = 2;
	// private final static int kCDXTranslation_Any = 3;

	// private final static int kCDXAbundance_Unspecified = 0;
	// private final static int kCDXAbundance_Any = 1;
	// private final static int kCDXAbundance_Natural = 2;
	// private final static int kCDXAbundance_Enriched = 3;
	// private final static int kCDXAbundance_Deficient = 4;
	// private final static int kCDXAbundance_Nonnatural = 5;

	// private final static int kCDXExternalConnection_Unspecified = 0;
	// private final static int kCDXExternalConnection_Diamond = 1;
	// private final static int kCDXExternalConnection_Star = 2;
	// private final static int kCDXExternalConnection_PolymerBead = 3;
	// private final static int kCDXExternalConnection_Wavy = 4;

	// private final static int kCDXAtomGeometry_Unknown = 0;
	// private final static int kCDXAtomGeometry_1Ligand = 1;
	// private final static int kCDXAtomGeometry_Linear = 2;
	// private final static int kCDXAtomGeometry_Bent = 3;
	// private final static int kCDXAtomGeometry_TrigonalPlanar = 4;
	// private final static int kCDXAtomGeometry_TrigonalPyramidal = 5;
	// private final static int kCDXAtomGeometry_SquarePlanar = 6;
	// private final static int kCDXAtomGeometry_Tetrahedral = 7;
	// private final static int kCDXAtomGeometry_TrigonalBipyramidal = 8;
	// private final static int kCDXAtomGeometry_SquarePyramidal = 9;
	// private final static int kCDXAtomGeometry_5Ligand = 10;
	// private final static int kCDXAtomGeometry_Octahedral = 11;
	// private final static int kCDXAtomGeometry_6Ligand = 12;
	// private final static int kCDXAtomGeometry_7Ligand = 13;
	// private final static int kCDXAtomGeometry_8Ligand = 14;
	// private final static int kCDXAtomGeometry_9Ligand = 15;
	// private final static int kCDXAtomGeometry_10Ligand = 16;

	private final static int kCDXBondOrder_Single = 0x0001;
	private final static int kCDXBondOrder_Double = 0x0002;
	private final static int kCDXBondOrder_Triple = 0x0004;
	private final static int kCDXBondOrder_Quadruple = 0x0008;
	private final static int kCDXBondOrder_Quintuple = 0x0010;
	private final static int kCDXBondOrder_Sextuple = 0x0020;
	private final static int kCDXBondOrder_Half = 0x0040;
	private final static int kCDXBondOrder_OneHalf = 0x0080;
	private final static int kCDXBondOrder_TwoHalf = 0x0100;
	private final static int kCDXBondOrder_ThreeHalf = 0x0200;
	private final static int kCDXBondOrder_FourHalf = 0x0400;
	private final static int kCDXBondOrder_FiveHalf = 0x0800;
	private final static int kCDXBondOrder_Dative = 0x1000;
	private final static int kCDXBondOrder_Ionic = 0x2000;
	private final static int kCDXBondOrder_Hydrogen = 0x4000;
	private final static int kCDXBondOrder_ThreeCenter = 0x8000;
	// private final static int kCDXBondOrder_SingleOrDouble = kCDXBondOrder_Single
	// | kCDXBondOrder_Double
	// private final static int kCDXBondOrder_SingleOrAromatic =
	// kCDXBondOrder_Single | kCDXBondOrder_OneHalf
	// private final static int kCDXBondOrder_DoubleOrAromatic =
	// kCDXBondOrder_Double | kCDXBondOrder_OneHalf
	// private final static int kCDXBondOrder_Any = -1;

	private final static int kCDXBondDisplay_Solid = 0;
	private final static int kCDXBondDisplay_Dash = 1;
	private final static int kCDXBondDisplay_Hash = 2;
	private final static int kCDXBondDisplay_WedgedHashBegin = 3;
	private final static int kCDXBondDisplay_WedgedHashEnd = 4;
	private final static int kCDXBondDisplay_Bold = 5;
	private final static int kCDXBondDisplay_WedgeBegin = 6;
	private final static int kCDXBondDisplay_WedgeEnd = 7;
	private final static int kCDXBondDisplay_Wavy = 8;
	private final static int kCDXBondDisplay_HollowWedgeBegin = 9;
	private final static int kCDXBondDisplay_HollowWedgeEnd = 10;
	private final static int kCDXBondDisplay_WavyWedgeBegin = 11;
	private final static int kCDXBondDisplay_WavyWedgeEnd = 12;
	private final static int kCDXBondDisplay_Dot = 13;
	private final static int kCDXBondDisplay_DashDot = 14;

	// private final static int kCDXBondDoublePosition_AutoCenter = 0x0000;
	// private final static int kCDXBondDoublePosition_AutoRight = 0x0001;
	// private final static int kCDXBondDoublePosition_AutoLeft = 0x0002;
	// private final static int kCDXBondDoublePosition_UserCenter = 0x0100;
	// private final static int kCDXBondDoublePosition_UserRight = 0x0101;
	// private final static int kCDXBondDoublePosition_UserLeft = 0x0102;

	// private final static int kCDXBondTopology_Unspecified = 0;
	// private final static int kCDXBondTopology_Ring = 1;
	// private final static int kCDXBondTopology_Chain = 2;
	// private final static int kCDXBondTopology_RingOrChain = 3;

	// private final static int kCDXBondReactionParticipation_Unspecified = 0;
	// private final static int kCDXBondReactionParticipation_ReactionCenter = 1;
	// private final static int kCDXBondReactionParticipation_MakeOrBreak = 2;
	// private final static int kCDXBondReactionParticipation_ChangeType = 3;
	// private final static int kCDXBondReactionParticipation_MakeAndChange = 4;
	// private final static int kCDXBondReactionParticipation_NotReactionCenter = 5;
	// private final static int kCDXBondReactionParticipation_NoChange = 6;
	// private final static int kCDXBondReactionParticipation_Unmapped = 7;

	// private final static int kCDXTextJustification_Right = -1;
	// private final static int kCDXTextJustification_Left = 0;
	// private final static int kCDXTextJustification_Center = 1;
	// private final static int kCDXTextJustification_Full = 2;
	// private final static int kCDXTextJustification_Above = 3;
	// private final static int kCDXTextJustification_Below = 4;
	// private final static int kCDXTextJustification_Auto = 5;
	// private final static int kCDXTextJustification_BestInitial = 6;

	public static void openDocument(StringBuffer sb) {
		sb.append("<?xml version=\"1.0\"?>\n");
	}

	static protected void openTag(StringBuffer sb, String name) {
		sb.append("<").append(name).append(">\n");
	}

	static protected void startOpenTag(StringBuffer sb, String name) {
		sb.append("<").append(name);
	}

	static protected void terminateTag(StringBuffer sb) {
		sb.append(">\n");
	}

	static protected void terminateEmptyTag(StringBuffer sb) {
		sb.append("/>\n");
	}

	static protected void appendEmptyTag(StringBuffer sb, String name, String[] attributes) {
		startOpenTag(sb, name);
		addAttributes(sb, attributes);
		terminateEmptyTag(sb);
	}

	static protected void addAttributes(StringBuffer sb, String[] attributes) {
		for (int i = 0; i < attributes.length; i++) {
			addAttribute(sb, attributes[i], attributes[++i]);
		}
	}

	static protected void addAttribute(StringBuffer sb, String key, String val) {
		// System.out.println("CMLW addAttribute " + key + "=" + val);
		sb.append(" ").append(key).append("=").append(esc(val));
	}

	static protected void closeTag(StringBuffer sb, String name) {
		sb.append("</").append(name).append(">\n");
	}

	private final static String escapable = "\\\\\tt\rr\nn\"\"";

	public static String esc(String str) {
		if (str == null || str.length() == 0)
			return "\"\"";
		boolean haveEscape = false;
		int i = 0;
		for (; i < escapable.length(); i += 2)
			if (str.indexOf(escapable.charAt(i)) >= 0) {
				haveEscape = true;
				break;
			}
		if (haveEscape)
			while (i < escapable.length()) {
				int pt = -1;
				char ch = escapable.charAt(i++);
				char ch2 = escapable.charAt(i++);
				StringBuffer sb = new StringBuffer();
				int pt0 = 0;
				while ((pt = str.indexOf(ch, pt + 1)) >= 0) {
					sb.append(str.substring(pt0, pt)).append('\\').append(ch2);
					pt0 = pt + 1;
				}
				sb.append(str.substring(pt0, str.length()));
				str = sb.toString();
			}
		return "\"" + escUnicode(str) + "\"";
	}

	public static String escUnicode(String str) {
		for (int i = str.length(); --i >= 0;)
			if (str.charAt(i) > 0x7F) {
				String s = "0000" + Integer.toHexString(str.charAt(i));
				str = str.substring(0, i) + "\\u" + s.substring(s.length() - 4) + str.substring(i + 1);
			}
		return str;
	}

	@Override
	public String toString() {
		return (sb.toString());
	}

}
