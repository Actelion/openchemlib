package com.actelion.research.util;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;

import java.awt.*;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;


/**
 *
 * <p>Title: Actelion Java Library </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: Actelion Ltd.</p>
 * @author Christian Rufener
 * @version 1.0
 */
public class Sketch
{
    private static final boolean debug_ = false;

    ///////////////////////////////////////////////////////////////////////////////
    //  A sketch unit is one decipoint. One decipoint is 1/10 of a point; one point is
    //	1/72 of a inch. The origin of the coordinate axes is at the  top left corner of the
    //	screen. The x coordinates increase fro left to right, and the y coordinates
    //	increase from top to bottom of the screen.
    //								- ISIS Sketch File Format
    ///////////////////////////////////////////////////////////////////////////////
    public static final double PI = 3.1415926;
    public static final int MAXMOLS = 12;
    public static final byte $Version = (byte)1;
    public static final byte $Totobjs = (byte)2;
    public static final byte $Obj = (byte)3;
    public static final byte $Locked = (byte)4;
    public static final byte $Pen_width = (byte)5;
    public static final byte $Pen_style = (byte)6;
    public static final byte $Pen_color = (byte)7;
    public static final byte $Transparent = (byte)8;
    public static final byte $Fill_style = (byte)9;
    public static final byte $Fill_color = (byte)10;
    public static final byte $Font = (byte)11;
    public static final byte $Parent = (byte)12;
    public static final byte $Obj_coords = (byte)13;
    public static final byte $Crop_coords = (byte)14;
    public static final byte $Roundrect_curve = (byte)15;
    public static final byte $Arc_endpts = (byte)16;
    public static final byte $Poly_points = (byte)17;
    public static final byte $Poly_smoothed = (byte)18;
    public static final byte $Begsketch = (byte)19;
    public static final byte $Endsketch = (byte)20;
    public static final byte $MDLEditText = (byte)21;
    public static final byte $Atom_coords = (byte)22;
    public static final byte $Atom_type = (byte)23;
    public static final byte $Atom_list = (byte)24;
    public static final byte $Atom_alias = (byte)25;
    public static final byte $Atom_number = (byte)26;
    public static final byte $Atom_chg = (byte)27;
    public static final byte $Atom_rad = (byte)28;
    public static final byte $Atom_msdif = (byte)29;
    public static final byte $Atom_valence = (byte)30;
    public static final byte $Atom_rbcount = (byte)31;
    public static final byte $Atom_substcount = (byte)32;
    public static final byte $Atom_stereo_care = (byte)33;
    public static final byte $Atom_h0bit = (byte)34;
    public static final byte $Atom_unsat = (byte)35;
    public static final byte $Atom_value = (byte)36;
    public static final byte $Atom_dispflags = (byte)37;
    public static final byte $Atom_hpos = (byte)38;
    public static final byte $Bond_atoms = (byte)39;
    public static final byte $Bond_type = (byte)40;
    public static final byte $Bond_stereo_type = (byte)41;
    public static final byte $Bond_topo = (byte)42;
    public static final byte $Bond_qtopo = (byte)43;
    public static final byte $Bond_rxn_center = (byte)44;
    public static final byte $Bond_stereo_care = (byte)45;
    public static final byte $Bond_dbl_side = (byte)46;
    public static final byte $Bond_dbl_width = (byte)47;
    public static final byte $RxnAtch = (byte)48;
    public static final byte $RGroupNo = (byte)49;
    public static final byte $RLogic = (byte)50;
    public static final byte $SGroupAtch = (byte)51;
    public static final byte $SGroupName = (byte)52;
    public static final byte $Atom_rgroupAtch = (byte)53;
    public static final byte $SGroupType = (byte)54;
    public static final byte $SGroupLinkVal = (byte)55;
    public static final byte $ArrowDir = (byte)56;
    public static final byte $ArrowStyle = (byte)57;
    public static final byte $MetaData = (byte)58;
    public static final byte $Mol_type = (byte)59;
    public static final byte $Abbrev_atch = (byte)60;
    public static final byte $SGroupAtchPt = (byte)61;
    public static final byte $Atom_npos = (byte)62;
    public static final byte $Atom_aamapped = (byte)63;
    public static final byte $Name = (byte)64;
    public static final byte $Comment = (byte)65;
    public static final byte $Atom_fixed = (byte)66;
    public static final byte $3D_num_basis_objs = (byte)67;
    public static final byte $3D_basis_objs = (byte)68;
    public static final byte $3D_name = (byte)69;
    public static final byte $3D_minval = (byte)70;
    public static final byte $3D_maxval = (byte)71;
    public static final byte $3D_tolerance = (byte)72;
    public static final byte $3D_point_dist = (byte)73;
    public static final byte $3D_dihed_chiral = (byte)74;
    public static final byte $3D_exclus_radius = (byte)75;
    public static final byte $3D_point_dir = (byte)76;
    public static final byte $3D_atom_query = (byte)77;
    public static final byte $Atom_zcoord = (byte)78;
    public static final byte $Atom_exact_change = (byte)79;
    public static final byte $Atom_rxn_stereo = (byte)80;
    public static final byte $Bond_crossed = (byte)81;
    public static final byte $Bond_alt_stereo = (byte)82;
    public static final byte $3D_exclus_ignore = (byte)83;
    public static final byte $Pen_style_token = (byte)84;
    public static final byte $BigMetaData = (byte)85;
    public static final byte $Pen_RGB2Color = (byte)88;
    public static final byte $Fill_RGB2Color = (byte)89;
    public static final byte $Atom_tplatchpt = (byte)91;
    public static final byte $Is_A_Model = (byte)92;
    public static final byte $Atom_aamap_num = (byte)96;
    public static final byte $Atom_hshow = (byte)98;
    public static final byte $SGroupNewAtch = (byte)164;
    public static final byte $SGroupContext = (byte)202;
    public static final byte $Nostruct_label = (byte)203;
    public static final byte $Nostruct_regno = (byte)204;
    public static final byte $Circ_Arc_Points = (byte)205;
    public static final byte $Bond_hash_spacing = (byte)206;
    public static final byte $Bond_bond_spacing = (byte)207;
    public static final byte $Atom_margin_width = (byte)208;
    public static final byte $ArrowType = (byte)220;
    public static final byte $Atom_Orig_coords = (byte)223;
    public static final byte $Atom_AttachLen = (byte)224;
    public static final byte $Atom_DotWidth = (byte)225;
    public static final byte $SGroupBracketLen = (byte)226;
    public static final byte $3D_marker_rf = (byte)227;
    public static final byte $Model_Rotated = (byte)228;
    public static final byte $ArrowSize = (byte)229;
    public static final byte $ArrowShaftSp = (byte)230;
    public static final byte $Pen_widthUnit = (byte)231;
    public static final byte $Atom_symbol = (byte)232;
    public static final byte $Atom_NumSize = (byte)233;
    public static final byte $3D_exclus_dradius = (byte)234;
    public static final byte $Atom_can_reverse = (byte)235;
    public static final byte $3D_point_ddist = (byte)236;
    public static final byte $Obj_Mol = (byte)12;
    public static final byte $Obj_Chiral = (byte)19;
    public static final byte $Obj_Bond = (byte)32;
    public static final byte $Obj_Atom = (byte)33;
    public static final int MSDIFF_OFFSET = 19;
    private static String DefaultFontName = "Arial";
//	private static final int MAXMFCOMMENT = 5000;
    private static final int BYTESIZE = 1;
    private static final int WORDSIZE = 2;

    private static final int MOLSIZE = 1440;
    private static final int REACTIONSIZE = 4000;
    private static final int PLUSSIZE = 20;
    
    
    //        public static final int TAGNAMELEN (sizeof("MDLSK")/sizeof(char)+1);
    public static boolean createMolFromSketchFile(StereoMolecule mol, String szFileName)
        throws IOException
    {
        File f = new File(szFileName);
        byte[] buffer = new byte[(int)f.length()];
        FileInputStream is = new FileInputStream(f);
        is.read(buffer);
        is.close();

        return createMolFromSketchBuffer(mol, buffer);
    }

    public static boolean createMolFromSketchBuffer(StereoMolecule mol, byte[] buffer)
        throws IOException
    {
        ByteArrayInputStream b = new ByteArrayInputStream(buffer);
        LittleEndianDataInputStream fp = new LittleEndianDataInputStream(b);

        return getMolObjects(mol, fp);
    }

    public static boolean writeMolSketchFile(Molecule mol, String filename)
        throws IOException
    {
        FileOutputStream os = new FileOutputStream(filename);
        boolean ok = writeMolSketchFile(mol,os);
        os.close();
        return ok;
    }

    public static boolean writeMolSketchFile(Molecule mol, OutputStream os)
        throws IOException
    {
        boolean ok = false;
        byte[] buffer = createSketchFromMol(mol);
        if (buffer != null) {
            os.write(buffer);
            ok = true;
        }

        return ok;
    }

    public static boolean createReactionFromSketchFile(Reaction rxn, String szFileName)
        throws IOException
    {
        File f = new File(szFileName);
        byte[] buffer = new byte[(int)f.length()];
        FileInputStream is = new FileInputStream(f);
        is.read(buffer);
        is.close();

        ByteArrayInputStream b = new ByteArrayInputStream(buffer);
        LittleEndianDataInputStream fp = new LittleEndianDataInputStream(b);

        return createReactionFromSketchBuffer(rxn, fp);
    }

    public static boolean createReactionFromSketchBuffer(Reaction rxn, byte[] buffer)
        throws IOException
    {
        ByteArrayInputStream b = new ByteArrayInputStream(buffer);
        LittleEndianDataInputStream fp = new LittleEndianDataInputStream(b);

        return createReactionFromSketchBuffer(rxn, fp);
    }

    public static boolean writeReactionSketchFile(Reaction rxn, String filename)
        throws IOException
    {
        boolean ok = false;
        FileOutputStream os = new FileOutputStream(filename);
        if (writeReactionSketchFile(rxn,os)) {
            os.flush();
            os.close();
            ok = true;
        }

        return ok;
    }

    public static boolean writeReactionSketchFile(Reaction rxn, OutputStream os)
    throws IOException
    {
        boolean ok = false;
        byte[] buffer = createSketchFromReaction(rxn);
    
        if (buffer != null) {
            os.write(buffer);
            ok = true;
        }
    
        return ok;
    }

    
    /*
        private int getMolObjectsOl(StereoMolecule mol,LittleEndianDataInputStream fp,int molAtoms[]) throws IOException
        {
            byte buff[] = new byte[1024];
            byte c;
            byte bval;
            short at1,at2,AtomType=0;
            int i;
            short t,val = 0,num=0;
            DblPoint dblCoords = new DblPoint();
            int ObjectCount=0;
            int MolObject=-1,BondObject=-1,AtomObject=-1;
            String s;

            do {
                c = fp.readByte();

                switch (c) {
                    case $Version:
                        ObjectCount=0;
                        MolObject = -1;
                        t = fp.readShort();
                        bval = fp.readByte();
                        break;

                    case $Totobjs:
                        t = fp.readShort();
                        val = fp.readShort();
                        ObjectCount = -1; AtomObject=-1;BondObject=-1;
                        break;

                    case $Obj:
                        t = fp.readShort();
                        bval = fp.readByte();
                        switch(bval) {
                            case $Obj_Mol:
                                MolObject++;AtomObject=-1;BondObject=-1;
                                AtomType = 0;
                                ObjectCount = -1;
                                break;
                            case $Obj_Chiral:
                                //                        mol.SetChiral(true);
                                mol.setChiral(true);
                                ObjectCount++;
                                break;
                            case $Obj_Bond:
                                BondObject++;
                                //System.out.println("mol.vectBonds.insert(mol.vectBonds.end(),bond);");
                                //                        mol.vectBonds.insert(mol.vectBonds.end(),bond);
                                ObjectCount++;
                                AtomType = 0;
                                break;
                            case $Obj_Atom:
                                AtomObject++;
                                //CXR ChemDraw Problem C atoms are not defined with symbol and number!
                                //                        strcpy(atom.Symbol,"C");
                                //                        atom.Number = 6;
                                //                        mol.vectAtoms.insert(mol.vectAtoms.end(),atom);
                                mol.addAtom("C");
                                ObjectCount++;
                                AtomType = 0;
                                break;

                            default:

                                ObjectCount++;
                            break;
                        }
                        if (MolObject > 1) {
                            System.err.println("Sorry Only One Molecule supported\n");
                            return 1;
                        }
                        break;

                    case $Begsketch: // $BeginSketch
                        t = fp.readShort();
                        fp.readFully(buff,0,(int)t-2);
                        //                fp.BufRead(buff,t-2);
                        break;

                    case $Endsketch: //$EndSketch
                        t = fp.readShort();
    //                    return 0;
                        break;
                    case $Atom_coords: //$Atom_coords
                        t = fp.readShort();
                        dblCoords.x = fp.readFloat();
                        dblCoords.y = fp.readFloat();
                        //                fp.BufRead(&dblCoords,sizeof(dblCoords));
                        mol.setAtomX(AtomObject,dblCoords.x);
                        //                mol.vectAtoms[AtomObject].x = dblCoords.x;
                        mol.setAtomY(AtomObject,dblCoords.y);
                        //                mol.vectAtoms[AtomObject].y = dblCoords.y;
                        mol.setAtomZ(AtomObject,0);
                        //                mol.vectAtoms[AtomObject].z = 0;
                        molAtoms[AtomObject] = ObjectCount;

                        //                mol.vectAtoms[AtomObject].Number = ObjectCount;
                        //                mol.vectAtoms[AtomObject].Charge = 0;
                        break;

                    case $Atom_type: //$Atom Type
                        t = fp.readShort();
                        AtomType = fp.readShort();
                        //                fp.BufRead(&AtomType,2);
                        break;

                    case $Atom_symbol:
                        t = fp.readShort();
                        fp.readFully(buff,0,t-2);
                        //                fp.BufRead(buff,t-2);
                        s = P2CStr(buff);
                        System.out.println("Atom Symbol Length " + t + " sym = " + P2CStr(buff));
                        //                                   %d Sym=[%s]\n",t,(char *)P2CStr(buff));
                        switch (AtomType) {
                            case 270:
                                System.out.println("Atom list");
                                //                        strcpy(mol.vectAtoms[AtomObject].Symbol,"L");
                                break;
                            case 271:
                                System.out.println("Atom not list");
                                //                        mol.vectAtoms[AtomObject].bNotList = TRUE;
                                //                        strcpy(mol.vectAtoms[AtomObject].Symbol,"L");
                                break;
                            default:
                                System.out.println("Atom Symbol ");
                            mol.setAtomicNo(AtomObject,Molecule.getAtomicNoFromLabel(s));
                            //                        strcpy(mol.vectAtoms[AtomObject].Symbol,P2CStr(buff));
                            break;
                        }
                        break;

                    case $Atom_list: //$Atom_List
                        t = fp.readShort();
                        num = fp.readShort();
                        //                fp.BufRead(&num,1); // read number of atoms in list
                        for (i = 0; i < num ; i++) {
                            val = fp.readShort(); // read atoms in list
                            //                    mol.vectAtoms[AtomObject].List[i] =  (int)val;
                            //TRACE("Atom %d in List\n",(int)val);
                        }
                        //                mol.vectAtoms[AtomObject].List[i] = -1;
                        break;

                    case $Atom_alias: //$Atom Alias
                        t = fp.readShort();
                        fp.readFully(buff,0,t-2);
                        //                fp.BufRead(buff,t-2);
                        break;

                    case $Atom_number: //Atom_number
                        t = fp.readShort();
                        val = fp.readShort();
                        break;

                    case $Atom_chg: //$Atom_Charge
                        t = fp.readShort();
                        bval = fp.readByte();
                        mol.setAtomCharge(AtomObject,((int)bval)-16);
    //                    System.out.println(" mol.vectAtoms[AtomObject].Charge = ((int)bval)-16;");
                        //                mol.vectAtoms[AtomObject].Charge = ((int)bval)-16;
                        break;

                    case $Atom_rad: //Atom_rad
                        t = fp.readShort();
                        bval = fp.readByte();
                        //TRACE("Atom Radical %d \n",bval);
                        break;

                    case $Atom_msdif: //$Atom_msdif
                        t = fp.readShort();
                        bval = fp.readByte();
                        //TRACE("Atom Mass Diff %d \n",bval);
                        break;

                    case $Atom_valence: //$Atom_valence
                        t = fp.readShort();
                        bval = fp.readByte();
                        //TRACE("Atom Valence %d \n",bval);
                        break;


                    case $Bond_atoms: //$Bond_atoms
                        t = fp.readShort();
                        at1 = fp.readShort();
                        at2 = fp.readShort();
                        mol.addBond(at1,at2,0);
                        //                fp.BufRead(&at1,2);
                        //                fp.BufRead(&at2,2);
    //                    System.out.println(" mol.vectBonds[BondObject].Atoms[0] = at1;");
                        //                mol.vectBonds[BondObject].Atoms[0] = at1;
                        //                mol.vectBonds[BondObject].Atoms[1] = at2;
                        //                mol.vectBonds[BondObject].Type = 1;                                // Default
                        break;

                    case $Bond_type: //$Bond_type
                        t = fp.readShort();
                        bval = fp.readByte();
                        setBondFeatures(mol,BondObject,(int)bval,0,0);
    //                    mol.setBondType(BondObject,getBondType((int)bval));
    //                    System.out.println("mol.vectBonds[BondObject].Type = (int)bval;");
                        //                mol.vectBonds[BondObject].Type = (int)bval;
                        break;

                    case $Bond_stereo_type: //Bond_stereo_type
                        t = fp.readShort();
                        bval = fp.readByte();
                        setBondFeatures(mol,BondObject,0,(int)bval,0);
    //                    System.out.println(" mol.vectBonds[BondObject].Stereo = (int)bval;");
                        //                mol.vectBonds[BondObject].Stereo = (int)bval;
                        break;


                    case $Atom_aamap_num: // Atom_aamap_num
                        t = fp.readShort();
                        val = fp.readShort();
    //                    System.out.println("mol.vectAtoms[AtomObject].AtomMap = (int)val;");
                        mol.setAtomMapNo(AtomObject,(int)val,false);
                        //                mol.vectAtoms[AtomObject].AtomMap = (int)val;
                        break;

                    case -1:
                        break;

                    default:
                        t = fp.readShort();
                    fp.readFully(buff,0,t-2);
                    //                fp.BufRead(buff,t-2);
                    break;

                }
            } while (c != -1);
            CorrectBondAtomNumbering(mol,AtomObject+1, BondObject+1,molAtoms);
            return 1;
        }
    */
    private static boolean getMolObjects(StereoMolecule molMain, LittleEndianDataInputStream fp)
        throws IOException
    {
        byte c;
        byte bval;
        short t;
        int ObjectCount = 0;
        boolean ok = false;
        byte[] buff = new byte[1024];
        boolean chiralFound = false;

        do {
            ;
            c = fp.readByte();

            switch (c) {
                case $Version:
                    ObjectCount = 0;
                    t = fp.readShort();
                    bval = fp.readByte();

                    break;

                case $Totobjs:
                    t = fp.readShort();
                    fp.readShort();    // val
                    ObjectCount = -1;

                    break;

                case $Font:
                    t = fp.readShort();
                    bval = fp.readByte(); // String length
                    fp.readFully(buff, 0, (int)bval); // Font name
                    fp.readShort(); // font size in decipoints
                    fp.readByte(); // emphasis

                    break;

                case $Endsketch: //$EndSketch

                    try {
                        t = fp.readShort();
                    } catch (IOException e) {
                        System.err.println("Error found $Endsketch, but no data...");
                    }

                    break;

                case $Obj:
                    t = fp.readShort();
                    bval = fp.readByte();

                    switch (bval) {
                        case $Obj_Mol:

                            StereoMolecule m = new StereoMolecule();

                            if (!getMolObject(m, fp)) {
                                return false;
                            }

                            molMain.addMolecule(m);

                            if (m.isFragment()) {
                                molMain.setFragment(true);
                            }

                            ok = true;

                            break;

                        case $Obj_Chiral:
                            // Maybe needs a FIX by doing also setChirality(cChiralityPure) ???
                            chiralFound = true;
                            ObjectCount++;

                            break;

                        default:
                            System.out.println(
                                "Warning: Sketch.getMolObjects() Object not supported: " + bval);
                            ObjectCount++;

                            break;
                    }

                    break;

                case 0:

                    //                    System.out.println("Read a byte of value '0'");
                    break;

                default:
                    t = fp.readShort();

                    if (t > 2) {
                        fp.readFully(buff, 0, (int)t - 2);
                    }

                    break;
            }
        } while ((c != -1) && (fp.in.available() > 0));

        if (!chiralFound)
            molMain.setToRacemate();

        return ok;
    }

    private static boolean createReactionFromSketchBuffer(
        Reaction rxn, LittleEndianDataInputStream fp) throws IOException
    {
        int numObjects;
        String strText;
        StereoMolecule mol = null;
        Rect rect = new Rect();
        Arrow arrow = new Arrow();
        ArrayList<StereoMolecule> v = new ArrayList<StereoMolecule>();

        // Read Header
        numObjects = readSketchHeader(fp);

        if (numObjects > 1) {
            for (int i = 0; i < numObjects; i++) {
                switch (getNextObject(fp)) {
                    case 9: // Text
                        strText = getTextObject(rect, fp);

                        if (strText == "+") {
                        }

                        break;

                    case $Obj_Mol: // Mol
                        mol = new StereoMolecule();
                        getMolObject(mol, fp);
                        mol.ensureHelperArrays(Molecule.cHelperNeighbours);
                        v.add(mol);

                        break;

                    case 14: //Arrow
                        getArrowObject(arrow, fp);

                        break;

                    case 4: //LineArrow
                        getLineArrowObject(arrow, fp);

                        break;

                    default:
                        System.err.println(
                            "Warning: Sketch.createReactionFromSketchBuffer(): UnKnown Object");

                        break;
                }
            }

            buildReaction(rxn, v, arrow);

            return true;
        } else {
            return false;
        }
    }

    private static void buildReaction(Reaction r, ArrayList<StereoMolecule> v, Arrow a)
    {
        int size = v.size();
        Point pt;
        StereoMolecule m = null;

        for (int i = 0; i < size; i++) {
            m = v.get(i);
            pt = getMoleculeCenter(m);

            if (a.left > pt.x) {
                r.addReactant(m);
            } else {
                r.addProduct(m);
            }
        }
    }

    private static Rect getBoundingRect(Molecule m)
    {
        Rect r = new Rect(Short.MAX_VALUE, Short.MAX_VALUE, Short.MIN_VALUE, Short.MIN_VALUE);
        int atoms = m.getAllAtoms();
        float x;
        float y;

        for (int i = 0; i < atoms; i++) {
            x = (short)m.getAtomX(i);
            y = (short)m.getAtomY(i);
            r.left = (short)Math.min(x, r.left);
            r.top = (short)Math.min(y, r.top);
            r.right = (short)Math.max(x, r.right);
            r.bottom = (short)Math.max(y, r.bottom);
        }
        debug("Bounding Rect for Molecule " + m + " is " + r.left + "," 
                           + r.top + "," + r.right + "," + r.bottom);
        return r;
    }

    private static Rect2D getBoundingRectangle(Molecule m)
    {
        Rect2D r = new Rect2D(
                Double.MAX_VALUE, Double.MAX_VALUE, Double.MIN_VALUE, Double.MIN_VALUE);
        int atoms = m.getAllAtoms();
        double x;
        double y;

        for (int i = 0; i < atoms; i++) {
            x = m.getAtomX(i);
            y = m.getAtomY(i);
            r.left = Math.min(x, r.left);
            r.top = Math.min(y, r.top);
            r.right = Math.max(x, r.right);
            r.bottom = Math.max(y, r.bottom);
        }

        debug("Bounding Rect2D for Molecule " + m + " is " + r.left + "," 
                           + r.top + "," + r.right + "," + r.bottom);
        return r;
    }

    private static Rect2D getBoundingRectangle(Reaction rxn)
    {
        Rect2D r = new Rect2D(
                Double.MAX_VALUE, Double.MAX_VALUE, Double.MIN_VALUE, Double.MIN_VALUE);
        Rect2D rm = null;
        int mols = rxn.getMolecules();

        for (int i = 0; i < mols; i++) {
            rm = getBoundingRectangle(rxn.getMolecule(i));
            r.left = Math.min(r.left, rm.left);
            r.top = Math.min(r.top, rm.top);
            r.right = Math.max(r.right, rm.right);
            r.bottom = Math.max(r.bottom, rm.bottom);
        }

        debug("Bounding Rect2D for Reaction " + rxn + " is " + r.left + "," 
                           + r.top + "," + r.right + "," + r.bottom);
        return r;
    }

    private static Point getMoleculeCenter(Molecule m)
    {
        Rect r = getBoundingRect(m);

        return getCenter(r);
    }

    private static Point getCenter(Rect r)
    {
        int dx = r.right - r.left;
        int dy = r.bottom - r.top;

        return new Point(r.left + (dx / 2), r.top + (dy / 2));
    }

    private static boolean getMolObject(StereoMolecule mol, LittleEndianDataInputStream fp)
        throws IOException
    {
        byte[] buff = new byte[1024];
        byte c;
        byte bval;
        short at1;
        short at2;
        short AtomType = 0;
        int i;
        short t;
        short val = 0;
        short num = 0;
        DblPoint dblCoords = new DblPoint();
        int ObjectCount = 0;
        int MolObject = -1;
        int BondObject = -1;
        int AtomObject = -1;
        String s;
        int[] molAtoms = new int[mol.getMaxAtoms()];
        boolean chiralFound = false;

        do {
            c = fp.readByte();
//			System.out.println("GetMolObject " + (int)c);
            switch (c) {
                case $Totobjs:
                    t = fp.readShort();
                    val = fp.readShort();
                    ObjectCount = -1;
                    AtomObject = -1;
                    BondObject = -1;

                    break;

                case $Obj:
                    t = fp.readShort();
                    bval = fp.readByte();

                    switch (bval) {
                        /*
                                                case $Obj_Mol:
                                                    MolObject++;AtomObject=-1;BondObject=-1;
                                                    AtomType = 0;
                                                    ObjectCount = -1;
                                                    break;
                        */
                        case $Obj_Bond:
                            BondObject++;
                            ObjectCount++;

                            break;

                        case $Obj_Atom:
                            AtomObject++;
                            mol.addAtom("C");
                            ObjectCount++;
                            AtomType = 0;

                            break;

                        case $Obj_Chiral:
                            ObjectCount++;
                            chiralFound = true;
                            AtomType = 0;

                            break;

                        default:
                            ObjectCount++;

                            break;
                    }

                    if (MolObject > 1) {
                        System.err.println("Sorry Only One Molecule supported\n");
                        throw new RuntimeException("Only One Molecule supported");
                    }

                    break;

                case $Begsketch: // $BeginSketch
                    t = fp.readShort();
                    fp.readFully(buff, 0, (int)t - 2);

                    //                fp.BufRead(buff,t-2);
                    break;

                case $Endsketch: //$EndSketch
                    t = fp.readShort();
                    correctBondAtomNumbering(mol, AtomObject + 1, BondObject + 1, molAtoms);

                    if (!chiralFound)
                        mol.setToRacemate();
                        
                    return true;

                case $Atom_coords: //$Atom_coords
                    t = fp.readShort();
                    dblCoords.x = fp.readFloat();
                    dblCoords.y = fp.readFloat();
                    mol.setAtomX(AtomObject, dblCoords.x);
                    mol.setAtomY(AtomObject, dblCoords.y);
                    mol.setAtomZ(AtomObject, 0);
                    molAtoms[AtomObject] = ObjectCount;

//    System.out.println("R #" + AtomObject + " " + dblCoords.x + " " + dblCoords.y);
                    break;

                case $Atom_type: //$Atom Type
                    t = fp.readShort();
                    AtomType = fp.readShort();

                    break;

                case $Atom_symbol:
                    t = fp.readShort();
                    fp.readFully(buff, 0, t - 2);
                    s = convertPascalString(buff);
                    mol.setAtomicNo(AtomObject, Molecule.getAtomicNoFromLabel(s));

                    break;

                case $Atom_list: //$Atom_List

                    boolean bNotList = false;

                    if (AtomType == 271) {
                        bNotList = true;
                    }

                    // Number of list entries
                    t = fp.readShort();
                    num = fp.readByte();

                    int[] v = new int[num];

                    for (i = 0; i < num; i++) {
                        val = fp.readShort(); // read atoms in list
                        v[i] = (int)val;
                    }

                    mol.setAtomList(AtomObject, v, bNotList);
                    mol.setFragment(true);

                    break;

                case $Atom_alias: //$Atom Alias
                    t = fp.readShort();
                    fp.readFully(buff, 0, t - 2);

                    break;

                case $Atom_number: //Atom_number
                    t = fp.readShort();
                    val = fp.readShort();

                    break;

                case $Atom_chg: //$Atom_Charge
                    t = fp.readShort();
                    bval = fp.readByte();
                    mol.setAtomCharge(AtomObject, ((int)bval) - 16);

                    break;

                case $Atom_rad: //Atom_rad
                    t = fp.readShort();
                    bval = fp.readByte();
                    mol.setAtomRadical(AtomObject, getRadicalType(bval));

                    //System.out.println("Atom Radical %d \n",bval);
                    break;

                case $Atom_msdif: //$Atom_msdif
                    t = fp.readShort();
                    bval = fp.readByte();

                    // ignore is if msdiff is 0
                    if (bval != 0) {
                        int iso = bval - MSDIFF_OFFSET;
                        int atomicNo = mol.getAtomicNo(AtomObject);
                        mol.setAtomMass(AtomObject, Molecule.cRoundedMass[atomicNo] + iso);
                    }

                    break;

                case $Atom_valence: //$Atom_valence
                    t = fp.readShort();
                    bval = fp.readByte();

                    //System.out.println("Atom Valence %d \n",bval);
                    break;

                case $Bond_atoms: //$Bond_atoms
                    t = fp.readShort();
                    at1 = fp.readShort();
                    at2 = fp.readShort();
                    mol.addBond(at1, at2, Molecule.cBondTypeSingle);

                    break;

                case $Bond_type: //$Bond_type
                    t = fp.readShort();
                    bval = fp.readByte();
                    setBondFeatures(mol, BondObject, (int)bval, 0, 0);

                    break;

                case $Bond_stereo_type: //Bond_stereo_type
                    t = fp.readShort();
                    bval = fp.readByte();
                    setBondFeatures(mol, BondObject, 0, (int)bval, 0);

                    break;

                case $Atom_aamap_num: // Atom_aamap_num
                    t = fp.readShort();
                    val = fp.readShort();
                    mol.setAtomMapNo(AtomObject, (int)val, false);

                    break;

                case -1:
                    break;

                default:
                    t = fp.readShort();

                    if (t > 2) {
                        fp.readFully(buff, 0, t - 2);
                    }

                    break;
            }
        } while (c != -1);

        return false;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Function:	FindAtomIndex
    // Purpose:		This function is used by the CorrectBondAtomNumbering Function
    //				and finds the 0 based index of the Atoms in the Molecule based
    //				on the sketch internal numbering
    // Parameters:	Molecule & to operate on
    //				int Sketch Atom number
    //
    // Return Val:	0 based Atom index
    ///////////////////////////////////////////////////////////////////////////////
    private static int findAtomIndex(StereoMolecule mol, int number, int[] molAtoms)
    {
        int i;
        int cnt = molAtoms.length;

        for (i = 0; i < cnt; i++) {
            if (molAtoms[i] == number) {
                return i;
            }
        }

        return -1;
    }

    private static boolean correctBondAtomNumbering(
            StereoMolecule mol, int atoms, int bonds, int[] molAtoms)
    {
        int i;
        int cnt = bonds;
        int a1;
        int a2;

        for (i = 0; i < cnt; i++) {
            a1 = findAtomIndex(mol, mol.getBondAtom(0, i), molAtoms);

            //            a1 = FindAtomIndex(mol,mol.vectBonds[i].Atoms[0],int molAtoms[]);
            a2 = findAtomIndex(mol, mol.getBondAtom(1, i), molAtoms);

            //            a2 = FindAtomIndex(mol,mol.vectBonds[i].Atoms[1],int molAtoms[]);
            if ((a1 == -1) || (a2 == -1)) {
                return false;
            }

            mol.setBondAtom(0, i, a1);

            //            mol.vectBonds[i].Atoms[0] = a1;
            mol.setBondAtom(1, i, a2);

            //            mol.vectBonds[i].Atoms[1] = a2;
        }

        cnt = atoms;

        for (i = 0; i < cnt; i++) {
            //            mol.vectAtoms[i].Number = mol.GetAtomNumberFromSymbol(i);
        }

        return true;
    }

    private static String convertPascalString(byte[] pstr)
    {
        if ((pstr == null) || (pstr[0] > 255)) {
            return null;
        }

        return new String(pstr, 1, pstr[0], StandardCharsets.UTF_8);
    }

    private static void setBondFeatures(
            StereoMolecule mMol, int bond, int bondType, int stereo, int topology)
    {
        int realBondType = 0;
        boolean isAtomESRAnd = false;

        switch (stereo) {
            case 1:
                realBondType = Molecule.cBondTypeUp;

                break;

            case 3:
                realBondType = Molecule.cBondTypeCross;

                break;

            case 4:
                realBondType = Molecule.cBondTypeUp;
                isAtomESRAnd = true;
                break;

            case 6:
                realBondType = Molecule.cBondTypeDown;

                break;

            default:

                switch (bondType) {
                    case 1:
                        realBondType = Molecule.cBondTypeSingle;

                        break;

                    case 2:
                        realBondType = Molecule.cBondTypeDouble;

                        break;

                    case 3:
                        realBondType = Molecule.cBondTypeTriple;

                        break;

                    case 4:
                        realBondType = Molecule.cBondTypeDelocalized;

                        break;
                }

                break;
        }

        if (realBondType > 0) {
            mMol.setBondType(bond, realBondType);
        }

        if (isAtomESRAnd)
            mMol.setAtomESR(mMol.getBondAtom(0, bond), Molecule.cESRTypeAnd, -1);

        int queryFeatures = 0;

        if (bondType > 4) {
            switch (bondType) {
                case 5:
                    queryFeatures |= (Molecule.cBondTypeSingle | Molecule.cBondTypeDouble);

                    break;

                case 6:
                    queryFeatures |= (Molecule.cBondTypeSingle | Molecule.cBondTypeDelocalized);

                    break;

                case 7:
                    queryFeatures |= (Molecule.cBondTypeDouble | Molecule.cBondTypeDelocalized);

                    break;

                case 8:
                    queryFeatures |= Molecule.cBondQFBondTypes;

                    break;
            }
        }

        if (topology == 1) {
            queryFeatures |= Molecule.cBondQFRing;
        }

        if (topology == 2) {
            queryFeatures |= Molecule.cBondQFNotRing;
        }

        if (queryFeatures != 0) {
            mMol.setBondQueryFeature(bond, queryFeatures, true);
        }
    }

    private static int getRadicalType(int rad)
    {
        int radical = 0;

        switch (rad) {
            case 1:
                radical = Molecule.cAtomRadicalStateS;

                break;

            case 2:
                radical = Molecule.cAtomRadicalStateD;

                break;

            case 3:
                radical = Molecule.cAtomRadicalStateT;

                break;

            default:
                radical = 0;

                break;
        }

        return radical;
    }

    private static int readSketchHeader(LittleEndianDataInputStream fp)
        throws IOException
    {
        byte[] buff = new byte[1024];
        byte c;
        short t;
        short val;

        while ((c = fp.readByte()) != -1) {
            switch (c) {
                case $Version: //Version
                    t = fp.readShort();
                    fp.readByte(); // bval

                    break;

                case $Totobjs: // TotObjects
                    t = fp.readShort();
                    val = fp.readShort();

                    //System.out.println("Total Objects %d\n",(int)val);
                    return (int)val;

                case $Comment:
                case $Name:
                    t = fp.readShort();
                    fp.readFully(buff, 0, t - 2);

                    break;

                default:
                    t = fp.readShort();
                    fp.readFully(buff, 0, t - 2);

                    return 0;
            }
        }

        return 0;
    }

    private static int getNextObject(LittleEndianDataInputStream fp)
        throws IOException
    {
        byte[] buff = new byte[1024];
        byte c;
        byte bval;
        short t;

        while ((c = fp.readByte()) != -1) {
            switch (c) {
                case $Obj:
                    t = fp.readShort();
                    bval = fp.readByte();

                    return bval;

                default:
                    t = fp.readShort();
                    fp.readFully(buff, 0, t - 2);

                    break;
            }
        }

        return 0;
    }

    private static String getTextObject(Rect rect, LittleEndianDataInputStream fp)
        throws IOException
    {
        String s = null;
        byte[] buff = new byte[1024];
        byte c;
        short t;

        fp.in.mark(1);

        while ((c = fp.readByte()) != -1) {
            switch (c) {
                case $Obj_coords: //$Obj_coords
                    t = fp.readShort();
                    rect.left = fp.readShort();
                    rect.top = fp.readShort();
                    rect.right = fp.readShort();
                    rect.bottom = fp.readShort();

                    break;

                case $MDLEditText: //$MDLEditText

                    //TRACE("CODE %u ",(unsigned int)c);
                    t = fp.readShort();

                    //fp.readFully(buff,0,t-2);
                    s = getMDLText(fp);

                    break;

                case $Transparent:
                case $Crop_coords:
                case $RxnAtch:
                case $Pen_RGB2Color:

                    //TRACE("CODE %u ",(unsigned int)c);
                    t = fp.readShort();
                    fp.readFully(buff, 0, t - 2);

                    break;

                default:

                    //TRACE("End of Text...\n");
                    // fp.PutBack(1);
                    // return nFound;
                    fp.in.reset();

                    return s;
            }
        }

        return s;
    }

    private static String getMDLText(LittleEndianDataInputStream fp)
        throws IOException
    {
        byte[] buff = new byte[1024];
        byte[] text = null;
        short nNumChar;

        fp.readShort();   // len
        fp.readFully(buff, 0, 12); // ignore 12 bytes
        nNumChar = fp.readShort();
        fp.readShort(); //ignore                //2
        fp.readShort(); //ignore num styles     //2
        fp.readShort(); //ignore                //2
        fp.readInt(); //ignore outline        //4
        fp.readShort(); //ignore                //2
        fp.readFully(buff, 0, 32); // ignore 32 bytes Font //32
        fp.readByte(); //ignore font styles        //1
        fp.readShort(); //ignore Font size;         //2

        text = new byte[nNumChar];

        for (int i = 0; i < nNumChar; i++) {
            text[i] = fp.readByte();
            fp.readByte(); // Jump over the attribute!!
        }

        return new String(text, StandardCharsets.UTF_8);
    }

    private static boolean getArrowObject(Arrow arrow, LittleEndianDataInputStream fp)
        throws IOException
    {
        int maxcount = 12;
        int count = 0;
        boolean bFound = false;
        byte c;
        short t;
        byte[] buff = new byte[1024];

        fp.in.mark(1);

        while ((c = fp.readByte()) != -1) {
            switch (c) {
                case $Obj_coords: //$Obj_coords
                    t = fp.readShort();

                    //fp.BufRead(&rc,sizeof(Rect));
                    arrow.left = fp.readShort();
                    arrow.top = fp.readShort();
                    arrow.right = fp.readShort();
                    arrow.bottom = fp.readShort();
                    bFound = true;

                    break;

                case $Pen_width:
                case $Crop_coords:
                case $ArrowDir:
                case $ArrowStyle:
                case $Pen_RGB2Color:
                case $Pen_style_token:
                case $ArrowType:
                case $ArrowSize:
                case $ArrowShaftSp:
                case $Pen_widthUnit:
                    t = fp.readShort();
                    fp.readFully(buff, 0, t - 2);

                    break;

                case $RxnAtch:
                    t = fp.readShort();
                    arrow.attach[0] = fp.readShort();
                    arrow.attach[1] = fp.readShort();

                    break;

                default:

                    //TRACE("End of Arrow...\n");
                    fp.in.reset();

                    return bFound;
            }

            if (++count >= maxcount) {
                return bFound;
            }
        }

        return bFound;
    }

    private static boolean getLineArrowObject(Arrow arrow, LittleEndianDataInputStream fp)
        throws IOException
    {
        int maxcount = 13;
        int count = 0;
        boolean bFound = false;
        byte c;
        short t;
        byte[] buff = new byte[1024];

        fp.in.mark(1);

        while ((c = fp.readByte()) != -1) {
            switch (c) {
                case $Obj_coords: //$Obj_coords
                    t = fp.readShort();

                    //fp.BufRead(&rc,sizeof(Rect));
                    arrow.left = fp.readShort();
                    arrow.top = fp.readShort();
                    arrow.right = fp.readShort();
                    arrow.bottom = fp.readShort();
                    bFound = true;

                    break;

                case $Pen_width:
                case $Fill_style:
                case $ArrowDir:
                case $ArrowStyle:
                case $Pen_style_token:
                case $Pen_RGB2Color:
                case $Fill_RGB2Color:
                case $ArrowType:
                case $ArrowSize:
                case $ArrowShaftSp:
                case $Pen_widthUnit:
                    t = fp.readShort();
                    fp.readFully(buff, 0, t - 2);

                    break;

                case $RxnAtch:
                    t = fp.readShort();
                    arrow.attach[0] = fp.readShort();
                    arrow.attach[1] = fp.readShort();

                    break;

                default:

                    //TRACE("End of Arrow...\n");
                    fp.in.reset();

                    return bFound;
            }

            if (++count >= maxcount) {
                return bFound;
            }
        }

        return bFound;
    }

    // Write interface
    private static int writeSketchHeader(LittleEndianDataOutputStream fp, int nObjects)
        throws IOException
    {
        byte c;

        writeProperty(fp, $Version, (byte)4); // $Version
        writeProperty(fp, $Totobjs, (short)nObjects); // $TotObjects
        Write97(fp);

        c = $Pen_RGB2Color; // $Pen_RGB2Color
        fp.writeByte(c);
        fp.writeShort((short)(3 + WORDSIZE));
        c = 0;
        fp.writeByte(c);
        fp.writeByte(c);
        fp.writeByte(c);

        c = $Fill_RGB2Color; // $Fill_RGB2Color
        fp.writeByte(c);
        fp.writeShort((short)(3 + WORDSIZE));
        c = 0;
        fp.writeByte(c);
        fp.writeByte(c);
        fp.writeByte(c);

        //System.out.println("@TODO: The font handling in the sketch needs a fix");
        int nLen = DefaultFontName.length();
        fp.writeByte($Font);
        fp.writeShort((short)(nLen + 4 + WORDSIZE)); // len + sizeof(byte) + sizeof(short) + sizeof(byte)
                                                     // string size + length byte + font size + emphasis
        fp.writeByte((byte)nLen);
        fp.writeBytes(DefaultFontName);
        fp.writeShort((short)120);
        fp.writeByte((byte)0);

        return 1;
    }

    private static int writeProperty(LittleEndianDataOutputStream fp, byte c, byte d)
        throws IOException
    {
        fp.writeByte(c);
        fp.writeShort((short)(1 + WORDSIZE)); // sizeof(byte)
        fp.writeByte(d);

        return 4;
    }

    private static int writeProperty(LittleEndianDataOutputStream fp, byte c, String sz)
        throws IOException
    {
        int nLen = sz.length();
        fp.writeByte(c);
        fp.writeShort((short)(nLen + 1 + WORDSIZE)); // sizeof(byte) + strlen()
        fp.writeByte((byte)nLen);
        fp.writeBytes(sz);
        return nLen + 5; // 3 byte + 1 short
    }

    private static int writeProperty(LittleEndianDataOutputStream fp, byte c, short d)
        throws IOException
    {
        fp.writeByte(c);
        fp.writeShort((short)(WORDSIZE + WORDSIZE)); // sizeof(short)
        fp.writeShort(d);

        return 5; // 2 short + byte
    }

    private static int writeProperty(LittleEndianDataOutputStream fp, byte c, DblPoint d)
        throws IOException
    {
        fp.writeByte(c);
        fp.writeShort((short)(8 + WORDSIZE)); // 4 byte per float
        fp.writeFloat(d.x);
        fp.writeFloat(d.y);

        return 8 + 3;
    }

    private static int writeProperty(LittleEndianDataOutputStream fp, byte c, Rect r)
        throws IOException
    {
        fp.writeByte(c);
        fp.writeShort((short)((4 * WORDSIZE) + WORDSIZE)); // sizeof(Rect) = 4 * short
        fp.writeShort(r.left);
        fp.writeShort(r.top);
        fp.writeShort(r.right);
        fp.writeShort(r.bottom);

        return 8 + 3;
    }

    private static int writeProperty(LittleEndianDataOutputStream fp, byte c)
        throws IOException
    {
        fp.writeByte(c);
        fp.writeShort((short)WORDSIZE);

        return 3;
    }

    private static int Write97(LittleEndianDataOutputStream fp)
        throws IOException
    {
        char c = 97;
        fp.writeByte((byte)c);
        fp.writeShort((short)3);
        c = 0;
        fp.writeByte((byte)c);

        return 4;
    }

    private static void translateMolecule(Molecule m)
    {
        Rect2D rc = getBoundingRectangle(m);
        double dx = 0;
        double dy = 0;

        if (rc.left < 0) {
            dx = 0 - rc.left;
        }

        if (rc.top < 0) {
            dy = 0 - rc.top;
        }

        translateMolecule(m, dx, dy);
    }

    private static void translateMolecule(Molecule m, double dx, double dy)
    {
        int atoms = m.getAllAtoms();
        debug("Translate molecule " + dx + " " + dy);
        for (int i = 0; i < atoms; i++) {
            m.setAtomX(i, m.getAtomX(i) + dx);
            m.setAtomY(i, m.getAtomY(i) + dy);
        }
    }

    private static void translateReaction(Reaction m)
    {
        Rect2D rc = getBoundingRectangle(m);
        double dx = 0;
        double dy = 0;

        if (rc.left < 0) {
            dx = 0 - rc.left;
        }

        if (rc.top < 0) {
            dy = 0 - rc.top;
        }

        int mols = m.getMolecules();
        for (int i = 0; i < mols; i++) {
            translateMolecule(m.getMolecule(i),dx,dy);
        }
    }

    private static void scaleMolecule(Molecule m, int maxsize)
    {
        Rect2D rc = getBoundingRectangle(m);
        double dx = Math.abs(rc.right - rc.left);
        double dy = Math.abs(rc.bottom - rc.top);
        double scale = (double)maxsize / Math.max(dx, dy);
        scaleMolecule(m, scale);
    }

    private static void scaleMolecule(Molecule m, double scale)
    {
        int atoms = m.getAllAtoms();

        for (int i = 0; i < atoms; i++) {
            m.setAtomX(i, m.getAtomX(i) * scale);
            m.setAtomY(i, m.getAtomY(i) * scale);
            m.setAtomZ(i, m.getAtomZ(i) * scale);
        }
    }

    private static double getScaleFactor(Reaction m)
    {
        Rect2D rc = getBoundingRectangle(m);
        double scale = REACTIONSIZE / (rc.right - rc.left);
        return scale;
    }
    private static void scaleReaction(Reaction m, double scale)
    {
        debug("Scaling Reaciton at " + scale);
        int mols = m.getMolecules();
        for (int i = 0; i < mols; i++) {
            scaleMolecule(m.getMolecule(i),scale);
        }
    }

    public static byte[] createSketchFromMol(Molecule m)
    {
        byte[] buffer = null;
        StereoMolecule mol = new StereoMolecule(m);
        try {
            if (mol.getAllAtoms() > 0) {
                ByteArrayOutputStream b = new ByteArrayOutputStream();
                LittleEndianDataOutputStream fp = new LittleEndianDataOutputStream(b);

                // Rect rc = new Rect((short)0,(short)0,(short)2000,(short)2000);
                writeSketchHeader(fp, 1);
                translateMolecule(mol);
                scaleMolecule(mol, MOLSIZE);

                writeMoleculeToSketch(fp, mol, 1);
                writeProperty(fp, $Endsketch); // EndSketch

                buffer = b.toByteArray();
            }
        } catch (IOException e) {
            System.err.println("Sketch.SketchFromMol() " + e);
        }

        return buffer;
    }



    public static byte[] createSketchFromReaction(Reaction r)
        throws IOException
    {
        Reaction rxn = new Reaction(r);
        int reactants = rxn.getReactants();
        int prods = rxn.getProducts();
        int cnt = prods + reactants;

        ByteArrayOutputStream b = new ByteArrayOutputStream();
        LittleEndianDataOutputStream fp = new LittleEndianDataOutputStream(b);

        translateReaction(rxn);
        double scale = getScaleFactor(rxn);
        scaleReaction(rxn,scale);
        writeSketchHeader(fp, ((cnt + reactants) - 1 + prods) - 1 + 1);

        for (int i = 0; i < reactants; i++) {
            Molecule mol = rxn.getMolecule(i);
            Rect rc = getBoundingRect(mol);
            debug("Drawing reactant " + i + " @ " + rc.left + "," + rc.top + "," + rc.right + "," + rc.bottom);
            writeMoleculeToSketch(fp, mol, (short)-1, (short)-1, 1);

            if (i < (reactants - 1)) {
                Rect rn = getBoundingRect(rxn.getMolecule(i + 1));
                Rect rp = calcMiddleRect(rc, rn, PLUSSIZE);
                debug("Writing plus at " + 
                                   rp.left + "," + rp.top + "," + rp.right + "," + rp.bottom);
                writePlus(fp, rp, -1, -1);
            }
        }

        if ((reactants > 0) && (prods > 0)) {
            Rect rr = getBoundingRect(rxn.getMolecule(reactants - 1));
            Rect rp = getBoundingRect(rxn.getMolecule(reactants));
            Rect ra = calcMiddleRect(rr, rp, 1);
            int dx = (rp.left - rr.right) / 4;
            ra.left -= dx;
            ra.right += dx;
            writeArrowObject(fp, ra, -1, -1);
        }

        for (int i = reactants; i < cnt; i++) {
            Molecule mol = rxn.getMolecule(i);
            Rect rc = getBoundingRect(mol);
            debug("Drawing prod " + i + " @ " + rc.left + "," + rc.top + "," + rc.right + "," + rc.bottom);
            writeMoleculeToSketch(fp, mol, (short)-1, (short)-1, 1);

            if (i < (cnt - 1)) {
                Rect rn = getBoundingRect(rxn.getMolecule(i + 1));
                Rect rp = calcMiddleRect(rc, rn, PLUSSIZE);
                writePlus(fp, rp, -1, -1);
            }
        }

        writeProperty(fp, $Endsketch); // EndSketch

        return b.toByteArray();
    }


    public static byte[] createSketchFromReactionOLd(Reaction rxn)
        throws IOException
    {
        int reactants = rxn.getReactants();
        int prods = rxn.getProducts();
        int cnt = prods + reactants;

        float factor;
        ByteArrayOutputStream b = new ByteArrayOutputStream();
        LittleEndianDataOutputStream fp = new LittleEndianDataOutputStream(b);

        //        for (int i = 0; i < cnt ; i++) {
        //            Molecule mol = rxn.getMolecule(i);
        //            scaleMolecule(mol,MOLSIZE/cnt);
        //        }
        Rect2D rxnRect = getBoundingRectangle(rxn);

        factor = (float)(rxnRect.right - rxnRect.left);
        factor = 5000 / factor;

        writeSketchHeader(fp, ((cnt + reactants) - 1 + prods) - 1 + 1);

        for (int i = 0; i < reactants; i++) {
            Molecule mol = rxn.getMolecule(i);
            Rect rc = getBoundingRect(mol);
            writeMoleculeToSketch(fp, mol, (short)-1, (short)-1, factor);

            if (i < (reactants - 1)) {
                Rect rn = getBoundingRect(rxn.getMolecule(i + 1));
                Rect rp = calcMiddleRect(rc, rn, 5);
                scaleRect(rp, factor);
                writePlus(fp, rp, -1, -1);
            }
        }

        if ((reactants > 0) && (prods > 0)) {
            Rect rr = getBoundingRect(rxn.getMolecule(reactants - 1));
            Rect rp = getBoundingRect(rxn.getMolecule(reactants));
            Rect ra = calcMiddleRect(rr, rp, 1);
            int dx = (rp.left - rr.right) / 4;
            ra.left -= dx;
            ra.right += dx;
            scaleRect(ra, factor);
            writeArrowObject(fp, ra, -1, -1);
        }

        for (int i = reactants; i < cnt; i++) {
            Molecule mol = rxn.getMolecule(i);
            Rect rc = getBoundingRect(mol);
            writeMoleculeToSketch(fp, mol, (short)-1, (short)-1, factor);

            if (i < (cnt - 1)) {
                Rect rn = getBoundingRect(rxn.getMolecule(i + 1));
                Rect rp = calcMiddleRect(rc, rn, 5);
                scaleRect(rp, factor);
                writePlus(fp, rp, -1, -1);
            }
        }

        writeProperty(fp, $Endsketch); // EndSketch

        return b.toByteArray();
    }

    private static void scaleRect(Rect rc, double scale)
    {
        rc.left *= scale;
        rc.top *= scale;
        rc.right *= scale;
        rc.bottom *= scale;
    }

    private static int writeMoleculeToSketch(
        LittleEndianDataOutputStream fp, Molecule mol, double scale)
        throws IOException
    {
        writeProperty(fp, $Obj, (byte)$Obj_Mol); // Mol Object
        writeSubSketch(fp, mol, scale);

        return 1;
    }

    private static int writeMoleculeToSketch(
        LittleEndianDataOutputStream fp, Molecule mol, short wFrom, short wTo, double scale)
        throws IOException
    {
        byte c;

        writeProperty(fp, $Obj, (byte)$Obj_Mol); // Mol Object

        c = $RxnAtch; // RxnAtch
        fp.writeByte(c);
        fp.writeShort((short)((2 * 2) + 2));
        fp.writeShort(wFrom);
        fp.writeShort(wTo);

        writeSubSketch(fp, mol, scale);

        return 1;
    }

    private static void writeSubSketch(LittleEndianDataOutputStream fp, Molecule mol, double scale)
        throws IOException
    {
        int nAtoms = mol.getAllAtoms();
        int nBonds = mol.getAllBonds();
        byte c;

        writeProperty(fp, $Begsketch); // BeginSketch
        writeProperty(fp, $Totobjs, (short)(nAtoms + nBonds)); // $TotObjects
        Write97(fp);

        for (int i = 0; i < nAtoms; i++) {
            writeAtom(fp, mol, i, scale);

            if (mol.getAtomCharge(i) != 0) {
                writeAtomCharge(fp, mol, i);
            }

            if (!mol.isNaturalAbundance(i)) {
                writeAtomIsotope(fp, mol, i);
            }

            if (mol.getAtomRadical(i) != 0) {
                writeAtomRadical(fp, mol, i);
            }

            if (mol.getAtomList(i) != null) {
                writeAtomList(fp, mol, i);
            }
        }

        for (int j = 0; j < nBonds; j++) {
            writeProperty(fp, $Obj, (byte)$Obj_Bond); // Bond Object
            writeProperty(fp, $Pen_width, (byte)7); // Pen Width

            c = $Bond_atoms; // Bond_Atoms
            fp.writeByte(c);
            fp.writeShort((short)((2 * 2) + 2));
            fp.writeShort((short)mol.getBondAtom(0, j));
            fp.writeShort((short)mol.getBondAtom(1, j));
            writeProperty(fp, $Bond_type, (byte)getMDLBondType(mol.getBondType(j), false)); // Bond Type
            writeProperty(fp, $Bond_stereo_type, (byte)getMDLBondType(mol.getBondType(j), true)); // Bond Stereo
        }

        if (mol.getChirality() != Molecule.cChiralityNotChiral
         && mol.getChirality() != Molecule.cChiralityMeso
         && mol.getChirality() != Molecule.cChiralityRacemic
         && mol.getChirality() != Molecule.cChiralityUnknown) {
            writeChiral(fp, mol, scale);
        }

        writeProperty(fp, $Endsketch); // EndSketch
    }

    private static void writeAtom(
        LittleEndianDataOutputStream fp, Molecule mol, int atom, double scale)
        throws IOException
    {
        DblPoint pt = new DblPoint();
        pt.x = (float)(mol.getAtomX(atom) * scale);
        pt.y = (float)(mol.getAtomY(atom) * scale);

//		debug("W # " + atom + " " + mol.getAtomX(atom) + " " + mol.getAtomY(atom));
        writeProperty(fp, $Obj, (byte)$Obj_Atom); // Atom Object
        writeProperty(fp, $Atom_coords, pt); // Coords
        writeProperty(fp, $Atom_symbol, mol.getAtomLabel(atom)); // Symbol
        writeProperty(fp, $Atom_number, (short)(atom + 1)); // Atom Number within Molecule
        writeProperty(fp, $Atom_dispflags, (short)0x0004); // Disp Flags
        writeProperty(fp, $Atom_margin_width, (short)15); // Atom Space
        if (mol.getAtomMapNo(atom) != 0)
            writeProperty(fp, $Atom_aamap_num, (short)mol.getAtomMapNo(atom)); // Atom-Atom Map
    }

    private static void writeChiral(LittleEndianDataOutputStream fp, Molecule mol, double scale)
        throws IOException
    {
        Rect r = getBoundingRect(mol);
        scaleRect(r, scale);
        writeProperty(fp, $Obj, (byte)$Obj_Chiral);
        r.left = (short)(r.right + 10);
        r.right = (short)(r.left + 300);
        r.top -= (short)50;
        r.bottom = (short)(r.top + 40);
        writeProperty(fp, $Obj_coords, r);
    }

    private static void writeAtomCharge(LittleEndianDataOutputStream fp, Molecule mol, int atom)
        throws IOException
    {
        int charge = mol.getAtomCharge(atom);
        writeProperty(fp, $Atom_chg, (byte)(charge + 16)); // Atom charge
    }

    private static void writeAtomIsotope(LittleEndianDataOutputStream fp, Molecule mol, int atom)
        throws IOException
    {
        int mass = mol.getAtomMass(atom);
        int atomicNo = mol.getAtomicNo(atom);
        int diff = mass - Molecule.cRoundedMass[atomicNo];
        int iso = diff + MSDIFF_OFFSET;
        writeProperty(fp, $Atom_msdif, (byte)(iso)); // Atom isotope
    }

    private static void writeAtomRadical(LittleEndianDataOutputStream fp, Molecule mol, int atom)
        throws IOException
    {
        int rad = 0;

        switch (mol.getAtomRadical(atom)) {
            case Molecule.cAtomRadicalStateS:
                rad = 1;

                break;

            case Molecule.cAtomRadicalStateD:
                rad = 2;

                break;

            case Molecule.cAtomRadicalStateT:
                rad = 3;

                break;

            default:
                break;
        }

        writeProperty(fp, $Atom_rad, (byte)(rad)); // Atom radical
    }

    private static void writeAtomList(LittleEndianDataOutputStream fp, Molecule mol, int atom)
        throws IOException
    {
        int[] atomList = mol.getAtomList(atom);

        if (atomList != null) {
            boolean notlist = (mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) == Molecule.cAtomQFAny;
            int size = atomList.length;

            if (size > 0) {
                writeProperty(fp, $Atom_type, (short)(notlist ? 271 : 270));
                writeProperty(fp, $Atom_list, (byte)size);

                for (int i = 0; i < size; i++) {
                    int atomno = atomList[i];
                    fp.writeShort(atomno);
                }
            }
        }
    }

    private static int writeArrowObject(
        LittleEndianDataOutputStream fp, Rect coords, int wFrom, int wTo)
        throws IOException
    {
        // Write Object Type
        writeProperty(fp, $Obj, (byte)4); // Arrow Object

        // Write Pen Width
        writeProperty(fp, $Pen_width, (byte)7); // Pen Width

        // Write Arrow Direction
        writeProperty(fp, $ArrowDir, (byte)1);

        // Write Arrow Style
        writeProperty(fp, $ArrowStyle, (byte)1);

        // Write Arrow Size
        fp.writeByte($ArrowSize);
        fp.writeShort(((2 * WORDSIZE) + (2 * BYTESIZE) + WORDSIZE)); // word + word + byte + byte
        fp.writeShort(10);
        fp.writeShort(10);
        fp.writeByte(2);
        fp.writeByte(2);

        // Write RxnAttch
        fp.writeByte($RxnAtch);
        fp.writeShort(((2 * WORDSIZE) + WORDSIZE));
        fp.writeShort(wFrom);
        fp.writeShort(wTo);

        // Write Object Coords
        writeProperty(fp, $Obj_coords, coords);

        return 1;
    }

    static int writeMDLText(LittleEndianDataOutputStream buf, String szSource)
        throws IOException
    {
        int i;
        int nFontLen = DefaultFontName.length();
        int TextLen = szSource.length();
        int TotalLen = (TextLen * 2) + 63;

        buf.writeShort((TotalLen - 2));

        for (i = 0; i < 12; i++)
            buf.writeByte(0);

        buf.writeShort(TextLen);
        buf.writeShort(0);
        buf.writeShort(1);
        buf.writeShort(0);

        buf.writeInt(0xC000); // fixed
        buf.writeShort(0x11); // Left | Top

        buf.writeBytes(DefaultFontName); //

        for (i = 0; i < (32 - nFontLen); i++)
            buf.writeByte(0);

        buf.writeByte(0); // Normal
        buf.writeShort(120); // Size

        for (i = 0; i < TextLen; i++) {
            buf.writeByte((byte)szSource.charAt(i));
            buf.writeByte(0);
        }

        return TotalLen;
    }

    static int writePlus(LittleEndianDataOutputStream fp, Rect rect, int wFrom, int wTo)
        throws IOException
    {
        String s = "+";

        //        LPSTR szMDLText = new char [strlen(s) * 2 + 100];
        //        int iTextLen = CreateMDLText(fp,s);
        int TextLen = s.length();
        int TotalLen = (TextLen * 2) + 63;

        // Write Object Type
        writeProperty(fp, $Obj, (byte)9); // Text

        // Write Transparent Tag
        writeProperty(fp, $Transparent); // $Transparent

        // Write RxnAttch
        fp.writeByte($RxnAtch);
        fp.writeShort(((2 * WORDSIZE) + WORDSIZE));
        fp.writeShort(wFrom);
        fp.writeShort(wTo);

        // Write ObjectCoords
        writeProperty(fp, $Obj_coords, rect);

        // $Mdl Text
        fp.writeByte($MDLEditText);
        fp.writeShort((TotalLen + 2));
        writeMDLText(fp, s);

        //        fp.BufWrite((void *)szMDLText,iTextLen);
        //        delete [] szMDLText;
        return 0;
    }

    private static Rect calcMiddleRect(Rect rc, Rect rn, int size)
    {
        Rect rp = new Rect(rc.right, rc.top, rn.left, rc.bottom);
        Point p = getCenter(rp);

        return new Rect(
            (short)(p.x - (size / 2)), (short)(p.y - (size / 2)), (short)(p.x + (size / 2)),
            (short)(p.y + (size / 2)));
    }

    static int getMDLBondType(int btype, boolean bstereo)
    {
        int order = 1;
        int stereo = 0;

        switch (btype) {
            case Molecule.cBondTypeSingle:
                order = 1;
                stereo = 0;

                break;

            case Molecule.cBondTypeDouble:
                order = 2;
                stereo = 0;

                break;

            case Molecule.cBondTypeTriple:
                order = 3;
                stereo = 0;

                break;

            case Molecule.cBondTypeDown:
                order = 1;
                stereo = 6;

                break;

            case Molecule.cBondTypeUp:
                order = 1;
                stereo = 1;

                break;

/*			case Molecule.cBondTypeMix:    // discontinued bond type
                order = 1;
                stereo = 4;

                break;  */

            case Molecule.cBondTypeCross:
                order = 2;
                stereo = 3;

                break;

            case Molecule.cBondTypeDelocalized:
                order = 4;
                stereo = 0;

                break;

            default:
                order = 1;
                stereo = 0;

                break;
        }

        if (bstereo) {
            return stereo;
        } else {
            return order;
        }
    }
    
    
    private static void debug(String s)
    {
        if (debug_)
            System.out.println(s);
    }
}


class DblPoint
{
    float x;
    float y;
}


class Rect
{
    public Rect(short l, short t, short r, short b)
    {
        left = l;
        top = t;
        right = r;
        bottom = b;
    }

    public Rect()
    {
    }

    short left;
    short top;
    short right;
    short bottom;
}
;


class Rect2D
{
    public Rect2D(double l, double t, double r, double b)
    {
        left = l;
        top = t;
        right = r;
        bottom = b;
    }

    public Rect2D()
    {
    }

    double left;
    double top;
    double right;
    double bottom;
}
;


class Arrow extends Rect
{
    int[] attach = new int[2];
}
