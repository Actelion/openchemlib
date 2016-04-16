/*
* Copyright (c) 1997 - 2016
* Actelion Pharmaceuticals Ltd.
* Gewerbestrasse 16
* CH-4123 Allschwil, Switzerland
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 3. Neither the name of the the copyright holder nor the
*    names of its contributors may be used to endorse or promote products
*    derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/

package com.actelion.research.gui.clipboard.external;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.gui.clipboard.NativeClipboardAccessor;
import com.actelion.research.util.LittleEndianDataOutputStream;

import java.io.*;

public class ChemDrawCDX
{
	   private static byte COLORTABLE[] = {
	        (byte)0x00,(byte)0x03,(byte)0x32,(byte)0x00,(byte)0x08,(byte)0x00,(byte)0xFF,(byte)0xFF,
	        (byte)0xFF,(byte)0xFF,(byte)0xFF,(byte)0xFF,(byte)0x00,(byte)0x00,(byte)0x00,(byte)0x00,
	        (byte)0x00,(byte)0x00,(byte)0xFF,(byte)0xFF,(byte)0x00,(byte)0x00,(byte)0x00,(byte)0x00,
	        (byte)0xFF,(byte)0xFF,(byte)0xFF,(byte)0xFF,(byte)0x00,(byte)0x00,(byte)0x00,(byte)0x00,
	        (byte)0xFF,(byte)0xFF,(byte)0x00,(byte)0x00,(byte)0x00,(byte)0x00,(byte)0xFF,(byte)0xFF,
	        (byte)0xFF,(byte)0xFF,(byte)0x00,(byte)0x00,(byte)0x00,(byte)0x00,(byte)0xFF,(byte)0xFF,
	        (byte)0xFF,(byte)0xFF,(byte)0x00,(byte)0x00,(byte)0xFF,(byte)0xFF,
	    };

	    private static final byte[] HEADER = {
	        (byte)0x56,(byte)0x6A,(byte)0x43,(byte)0x44,(byte)0x30,(byte)0x31,(byte)0x30,(byte)0x30,
	        (byte)0x04,(byte)0x03,(byte)0x02,(byte)0x01,(byte)0x00,(byte)0x00,(byte)0x00,(byte)0x00,
	        (byte)0x00,(byte)0x00,(byte)0x00,(byte)0x00,(byte)0x00,(byte)0x00,(byte)0x00,(byte)0x00,
	        (byte)0x00,(byte)0x00,(byte)0x00,(byte)0x00,
	    };

	private static final boolean debug = false;

    public static final int kCDXObj_Page = 0x8001;
    public static final int kCDXObj_Group = 0x8002;
    public static final int kCDXObj_Fragment = 0x8003;
    public static final int kCDXObj_Node = 0x8004;
    public static final int kCDXObj_Bond = 0x8005;
    
    private static void debug(String s)
    {
        if (debug)
            System.out.println(s);
    }

    public boolean parse(BufferedInputStream is)
    {
        return false;
    }

    public String build()
    {
        StringBuilder sb = new StringBuilder();

        return sb.toString();
    }

    public byte[] getChemDrawBuffer(StereoMolecule mol)
    {

        CDXDocument doc = new CDXDocument();
        doc.add(new CDPShowEnhAtomStereo(true));
        doc.add(new CDPShowAtomStereo(true));
        doc.add(new CDPShowBondStereo(true));
        CDXNode page = new CDXPage();
        CDXNode frag = new CDXFragment();

        doc.add(page);
        page.add(frag);

        MolGeom g = getTransform(mol);
        writeMolecule(mol,frag,g);

        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        LittleEndianDataOutputStream l = new LittleEndianDataOutputStream(bos);
        write(l,doc);
        try{
            bos.close();
        } catch(IOException ex){
            System.err.println("ChemDrawCDX::getChemDrawBuffer(): Error closing stream : " + ex);
        }

//        dump(doc);

        return bos.toByteArray();
    }

    public byte[] getChemDrawBuffer(Reaction rxn)
    {
        CDXDocument doc = new CDXDocument();
        int r = rxn.getReactants();
        int p = rxn.getProducts();
        int mn = rxn.getMolecules();

        doc.add(new CDPShowEnhAtomStereo(true));
        doc.add(new CDPShowAtomStereo(true));
        doc.add(new CDPShowBondStereo(true));
        CDXNode page = new CDXPage();

        doc.add(page);
        
        MolGeom g = getTransform(rxn);
        int rids[] = new int[r];
        for (int i = 0; i < r; i++) {
            StereoMolecule m = rxn.getReactant(i);
            CDXNode frag = new CDXFragment();
            writeMolecule(m,frag,g);
            rids[i] =  frag.getID();
            page.add(frag);
        }

        int pids[] = new int[p];
        for (int i = 0; i < p; i++) {
            StereoMolecule m = rxn.getProduct(i);
            CDXNode frag = new CDXFragment();
            writeMolecule(m,frag,g);
            pids[i] =  frag.getID();
            page.add(frag);
        }
        
        CDXReactionStep step = new CDXReactionStep(rids,pids);
        page.add(step);
        
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        LittleEndianDataOutputStream l = new LittleEndianDataOutputStream(bos);
        write(l,doc);
        try{
            bos.close();
        } catch(IOException ex){
            System.err.println("ChemDrawCDX::getChemDrawBuffer(): Error closing stream : " + ex);
        }

//        dump(doc);

        return bos.toByteArray();
    }

    class MolGeom
    {
        MolGeom(double x,double y,double s)
        {
            xmin = x;
            ymin = y;
            scale = s;
        }
        double xmin;
        double ymin;
        double scale;
    }

    private MolGeom getTransform(Reaction rxn)
    {
        double xmin = Double.MAX_VALUE;
        double xmax = Double.MIN_VALUE;
        double ymin = Double.MAX_VALUE;
        double ymax = Double.MIN_VALUE;
        double scale = Double.MAX_VALUE;
        int cnt = rxn.getMolecules();
        
        if (cnt < 1)
            return new MolGeom(0,0,1);
        
        for (int mi = 0; mi < cnt; mi++) {
            StereoMolecule m = rxn.getMolecule(mi);
            m.ensureHelperArrays(Molecule.cHelperCIP);
            int na = m.getAllAtoms();
             // if more than one atom....
             for(int i = 0;i < na;i++){
                 double x = m.getAtomX(i);
                 double y = m.getAtomY(i);
                 xmin = Math.min(xmin,x);
                 xmax = Math.max(xmax,x);
                 ymin = Math.min(ymin,y);
                 ymax = Math.max(ymax,y);
             }
             double dx = xmax - xmin;
             double dy = xmax - ymin;
             debug(String.format("Mofile dx = %f dy = %f",dx,dy));

             double xscale = (double)0x01000000 / dx;
             double yscale = (double)0x01000000 / dy;
             debug(String.format("Mofile scale x = %.16f y = %.16f",xscale,yscale));
             scale = Math.min(xscale,yscale);
        }
        return new MolGeom(xmin,ymin,scale);
    }
    
    private MolGeom getTransform(StereoMolecule mol)
    {
        StereoMolecule m = new StereoMolecule(mol);

        m.ensureHelperArrays(Molecule.cHelperCIP);
        int na = m.getAllAtoms();
        double xmin = m.getAtomX(0);
        double xmax = m.getAtomX(0);
        double ymin = m.getAtomY(0);
        double ymax = m.getAtomY(0);

        // if more than one atom....
        for(int i = 1;i < na;i++){
            double x = m.getAtomX(i);
            double y = m.getAtomY(i);
            xmin = Math.min(xmin,x);
            xmax = Math.max(xmax,x);
            ymin = Math.min(ymin,y);
            ymax = Math.max(ymax,y);
        }
        double dx = xmax - xmin;
        double dy = xmax - ymin;
        debug(String.format("Mofile dx = %f dy = %f",dx,dy));

        double xscale = (double)0x00A00000 / dx;
        double yscale = (double)0x00A00000 / dy;
        debug(String.format("Mofile scale x = %.16f y = %.16f",xscale,yscale));

        return new MolGeom(xmin,ymin,Math.min(xscale,yscale));
    }
    
    private void writeMolecule(StereoMolecule mol,CDXNode frag,MolGeom g )
    {
        StereoMolecule m = new StereoMolecule(mol);
        m.ensureHelperArrays(Molecule.cHelperCIP);
        int na = m.getAllAtoms();
        int nb = m.getAllBonds();
        CDXAtom atoms[] = new CDXAtom[na];
if (false) {
        double xmin = m.getAtomX(0);
        double xmax = m.getAtomX(0);
        double ymin = m.getAtomY(0);
        double ymax = m.getAtomY(0);

        // if more than one atom....
        for(int i = 1;i < na;i++){
            double x = m.getAtomX(i);
            double y = m.getAtomY(i);
            xmin = Math.min(xmin,x);
            xmax = Math.max(xmax,x);
            ymin = Math.min(ymin,y);
            ymax = Math.max(ymax,y);
        }
        double dx = xmax - xmin;
        double dy = xmax - ymin;
        debug(String.format("Mofile dx = %f dy = %f",dx,dy));

        double xscale = (double)0x00800000 / dx;
        double yscale = (double)0x00800000 / dy;
        debug(String.format("Mofile scale x = %.16f y = %.16f",xscale,yscale));

        debug("Before translate");
        for(int i = 0;i < na;i++){
            double x = m.getAtomX(i);
            double y = m.getAtomY(i);
            debug(String.format("O Atom x = %.16f y = %.16f",x,y));
        }
        double scale = Math.min(xscale,yscale);
        m.translateCoords((float)xmin, (float)ymin);
        m.scaleCoords((float)scale);
} else {
//        MolGeom g = getTransform(mol);
        m.translateCoords((float)-g.xmin, (float)-g.ymin);
        m.scaleCoords((float)g.scale);
}
//        debug("After translate");
//        for(int i = 0;i < na;i++){
//            double x = m.getAtomX(i);
//            double y = m.getAtomY(i);
//            debug(String.format("Atom x = %.16f y = %.16f",x,y));
//        }

        boolean isMethane = na == 1 && m.getAtomicNo(0) == 6;
        
        for(int i = 0;i < na;i++){
            double x = m.getAtomX(i);
            double y = m.getAtomY(i);
            debug(String.format("C Atom x = %.0f y = %.0f",x,y));
            CDXAtom atom = null;
            if (isMethane) {
                atom = new CDXMethane(0x00100000 + (int)x,0x00100000 + (int)y);
            } else {
                atom = new CDXAtom((short)m.getAtomicNo(i),0x00100000 + (int)x,0x00100000 + (int)y);
            }
            if(m.isAtomStereoCenter(i)){
                int grp = m.getAtomESRGroup(i);
                int esr = m.getAtomESRType(i);
                if (!m.isAtomConfigurationUnknown(i) && !m.getStereoProblem(i)) {
                    atom.setAtomESRStereo(esr);
                }
                debug("Group is " + grp);
                if(grp != -1){
                    atom.setAtomESRGroup(grp + 1);
                }
                if(m.getAtomCIPParity(i) != 0){
                    atom.setAtomCIP(m.getAtomCIPParity(i),m.getAtomPi(i) == 2,m.isAtomParityPseudo(i));
                }

                if(m.getAtomMass(i) != 0){
                    atom.setAtomMass((short)m.getAtomMass(i));
                }
                if(m.getAtomCharge(i) != 0){
                    atom.setAtomCharge((byte)m.getAtomCharge(i));
                }
                if(m.getAtomRadical(i) != 0){
                    atom.setAtomRadical((byte)m.getAtomRadical(i));
                }

            }
            //			CDXAtom at2 = new CDXAtom((short)8, 0x02000000, 0x00AC03BA);
            //            CDXAtom atom = new CDXAtom((short)7,(int)x,(int)y);
            atoms[i] = atom;
            frag.add(atom);
        }
        // For all bonds
        for(int i = 0;i < nb;i++){
            int a1 = m.getBondAtom(0,i);
            int a2 = m.getBondAtom(1,i);
            int o = m.getBondType(i);
            CDXBond b = new CDXBond(atoms[a1].getID(),atoms[a2].getID());
            b.setBondType(o);
            frag.add(b);
        }
    }

    private static int NUM = 0;

    private abstract class CDXNode
    {
        protected int id;
        protected int type;
        private java.util.List<CDXNode> children;
        protected java.util.Set<CDPProperty> properties;

        public CDXNode(int type)
        {
            id = NUM++;
            this.type = type;
        }

        public int getID()
        {
            return id;
        }

        public abstract String start();

        public abstract String end();

        public void startWrite(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeInt(id);
            writeProperties(o);
        }

        public final void endWrite(DataOutput o) throws IOException
        {
            o.writeShort(0);
        }

        protected void writeProperties(DataOutput o)
        {
            if(properties != null){
                for(CDPProperty p : properties){
                    try{
                        p.write(o);
                    } catch(IOException e){
                        System.err.println("Error writin property: " + e);
                    }
                }
            }
        }

        public void walk(Functor functor)
        {
            functor.start(this);
            if(children != null){
                for(CDXNode n : children){
                    n.walk(functor);
                }
            }
            functor.end(this);
        }

        public void add(CDXNode node)
        {
            if(children == null){
                children = new java.util.ArrayList<CDXNode>();
            }
            children.add(node);
        }

        public void add(CDPProperty prop)
        {
            if(properties == null){
                properties = new java.util.HashSet<CDPProperty>();
            }
            if(properties.contains(prop)){
                properties.remove(prop);
            }
            properties.add(prop);
        }

        // Output
        protected String startString()
        {
            try{
                ByteArrayOutputStream bos = new ByteArrayOutputStream();
                LittleEndianDataOutputStream l = new LittleEndianDataOutputStream(bos);
                startWrite(l);
                bos.close();
                byte buffer[] = bos.toByteArray();
                StringBuilder sb = new StringBuilder();
                for(int i = 0;i < buffer.length;i++){
                    byte b = buffer[i];
                    if(i != 0){
                        sb.append(" ");
                    }
                    if(i > 0 && (i % 8) == 0){
                        sb.append("\n");
                    }
                    sb.append(String.format("%02X",(int)(b & 0xFF)));
                }
                return sb.toString();
            } catch(IOException e){
                System.err.println("Error: " + e);
                return "#";
            }
        }
    }

     public class CDXDocument extends CDXNode
    {
        public CDXDocument()
        {
            super(0x44436A56);
        }

        public void startWrite(DataOutput o) throws IOException
        {
            //		Header
            //		* 8 bytes with the value "VjCD0100" (56 6A 43 44 30 31 30 30)
            //		* 4 bytes reserved (04 03 02 01)
            //		* 16 bytes reserved, set to zero (00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00).
            o.writeInt(type); // 4 byte !		

            o.writeInt(0x30303130); // 4 byte
            o.writeInt(0x01020304); // 4 byte
            o.writeLong(0); // 8 byte
            o.writeLong(0); // 8 byte

            // Color Table
            o.write(COLORTABLE);
//            o.writeShort(0x0300);
//            o.writeInt(0x00080032);
//            o.writeLong(0x0000FFFFFFFFFFFFL);
//            o.writeLong(0x0000FFFF00000000L);
//            o.writeLong(0x0000FFFFFFFF0000L);
//            o.writeLong(0x00000000FFFF0000L);
//            o.writeLong(0x00000000FFFFFFFFL);
//            o.writeLong(0xFFFF0000FFFFFFFFL);

            //	ShowBondStereo
//            o.writeShort(0x060D);
//            o.writeShort(1);
//            o.writeByte(1);

        }

        /*
                public void endWrite(DataOutput o) throws IOException
                {
                    o.writeShort(0);
                }
         */

        public String start()
        {
            return startString();
        }

        public String end()
        {
            return String.format("00 00");
        }
    }

    public class CDXPage extends CDXNode
    {
        public CDXPage()
        {
            super(kCDXObj_Page);
        }

        public String start()
        {
            return startString();
        }

        public String end()
        {
            return String.format("00 00");
        }
    }

    private class CDXFragment extends CDXNode
    {
        public CDXFragment()
        {
            super(kCDXObj_Fragment);
        }

        /*
                public void startWrite(DataOutput o) throws IOException
                {
                    o.writeShort(type);
                    o.writeInt(id);
                }
                public void endWrite(DataOutput o) throws IOException
                {
                    o.writeShort(0);
                }
         */

        public String start()
        {
            return startString();
        }

        public String end()
        {
            return String.format("00 00");
        }

    }

    private class CDXText extends CDXNode
    {
        String s ;
        CDXText(String s, int x, int y)
        {
            super(0x8006);
            this.s = s;
            this.add(new CDPPoint2D(x,y));
            this.add(new CDPText(s));
        }
        public String start()
        {
            return startString();
        }
        public String end()
        {
            return String.format("00 00");
        }
    }

    private class CDXMethane extends CDXAtom
    {
        public CDXMethane(int x,int y)
        {
            super((short)6,x,y);	
            this.add(new CDXText("CH4",x,y));
        }

        public String start()
        {
            return startString();
        }
        public String end()
        {
            return String.format("00 00");
        }

    }

    private class CDXAtom extends CDXNode
    {

        public CDXAtom(short ord,int x,int y)
        {
            super(kCDXObj_Node);
            this.add(new CDPPoint2D(x,y));
            this.add(new CDPElement(ord));
        }

        public void setAtomESRStereo(int stereo)
        {
            debug("setAtomESRStereo " + stereo);
            switch(stereo){
                case Molecule.cESRTypeAbs: //0;
                    this.add(new CDPAtomEnhStereoType((byte)0x02));
                    break;
                case Molecule.cESRTypeAnd: //1;
                    this.add(new CDPAtomEnhStereoType((byte)0x04));
                    break;
                case Molecule.cESRTypeOr: //2;
                    this.add(new CDPAtomEnhStereoType((byte)0x05));
                    break;
                default:
                    break;
            }
        }

        public void setAtomESRGroup(int grp)
        {
            this.add(new CDPAtomEnhGroup((short)grp));
        }

        public void setAtomCIP(int parity,boolean pi,boolean pseudo)
        {
            /*			
                        R = 2
                        S = 3;
                        r = 4
                        s = 5;
                        U = 0;
                        u = 6
                        N = 1;
             */
            if(!pi){
                switch(parity){
                    case Molecule.cAtomCIPParityRorM:
                        this.add(new CDPAtomCIP(pseudo ? (byte)4 : (byte)2));
                        break;
                    case Molecule.cAtomCIPParitySorP:
                        this.add(new CDPAtomCIP(pseudo ? (byte)5 : (byte)3));
                        break;
                    default:
                        this.add(new CDPAtomCIP((byte)6));
                }
            } else{
                // not yet supported	
            }

        }

        public void setAtomMass(short mass)
        {
            add(new CDPAtomIsotope(mass));
        }

        public void setAtomCharge(byte charge)
        {
            add(new CDPAtomCharge(charge));
        }

        public void setAtomRadical(byte rad)
        {
            add(new CDPAtomRadical(rad));
        }

        /*
                public void startWrite(DataOutput o) throws IOException
                {
                    o.writeShort(type);
                    o.writeInt(id);
                    writeProperties(o);
                }
                public void endWrite(DataOutput o) throws IOException
                {
                    o.writeShort(0);
                }
         */

        public String start()
        {
            return startString();
//			return String.format("%02X %02X\n",(type & 0x00FF),(type >> 8)) +
//			String.format("%02X %02X %02X %02X\n",(id & 0x00FF),(id >> 8),(id >> 16),(id >> 24)) +
//			pt.encode();

        }

        public String end()
        {
            return String.format("00 00");
        }

    }

    private class CDXBond extends CDXNode
    {
//        CDPBondBegin b;
//        CDPBondEnd e;
        public CDXBond(int id1,int id2)
        {
            super(kCDXObj_Bond);
            add(new CDPBondBegin(id1));
            add(new CDPBondEnd(id2));
//            b = new CDPBondBegin(id1);
//            e = new CDPBondEnd(id2);
        }

        /*
                public void startWrite(DataOutput o) throws IOException
                {
                    o.writeShort(type);
                    o.writeInt(id);
                    b.write(o);
                    e.write(o);
                }
                public void endWrite(DataOutput o) throws IOException
                {
                    o.writeShort(0);
                }
         */

        public String start()
        {
            return startString();
//			return String.format("%02X %02X\n",(type & 0x00FF),(type >> 8)) +
//			String.format("%02X %02X %02X %02X\n",(id & 0x00FF),(id >> 8),(id >> 16),(id >> 24)) +
//			b.encode() + "\n" + e.encode();

        }

        public void setBondType(int type)
        {
            CDPBondType t = new CDPBondType(type);
            this.add(t);
            switch(type){
                case Molecule.cBondTypeDown: //			= 0x0009;
                    this.add(new CDPBondDisplay((short)3)); // WedgeBegin
                    break;
                case Molecule.cBondTypeUp: //				= 0x0011;
                    this.add(new CDPBondDisplay((short)6)); //WedgedHashBegin
                    break;

            }
        }

        public String end()
        {
            return String.format("00 00");
        }

    }

    public class CDXReactionStep extends CDXNode
    {
        public CDXReactionStep(int rids[], int pids[])
        {
            super(	0x800E);
            CDPReactants rts = new CDPReactants(rids);
            CDPProducts pts = new CDPProducts(pids);
            add(rts);
            add(pts);

        }

        
        public String start()
        {
            return startString();
        }

        public String end()
        {
            return String.format("00 00");
        }
    }

// Properties
    private abstract class CDPProperty
    {
        protected int type;
        public CDPProperty(int type)
        {
            this.type = type;
        }

        public abstract void write(DataOutput o) throws IOException;
//		{
//			o.writeShort(type);
//		}

    }

    /*
        private class CDXNodeTpe extends CDPProperty
        {
            int nodetype;
            public CDXNodeTpe(int nodetype)
            {
                super(0x0400);
                this.nodetype = nodetype;
            }
        }
     */
    public class CDPAtomCharge extends CDPProperty
    {
        byte charge;
        CDPAtomCharge(byte charge)
        {
            super(0x0421);
            this.charge = charge;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(1);
            o.writeByte(charge);
        }
    }

    public class CDPAtomRadical extends CDPProperty
    {
        byte rad;
        CDPAtomRadical(byte rad)
        {
            super(0x0422);
            this.rad = rad;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(1);
            o.writeByte(rad);
        }
    }

    public class CDPAtomIsotope extends CDPProperty
    {
        short iso;
        CDPAtomIsotope(short iso)
        {
            super(0x0420);
            this.iso = iso;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(2);
            o.writeShort(iso);
        }
    }

    public class CDPAtomCIP extends CDPProperty
    {
        byte val;
        CDPAtomCIP(byte v)
        {
            super(0x0437);
            val = v;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(1);
            o.writeByte(val);
        }
    }

    public class CDPAtomEnhStereoType extends CDPProperty
    {
        byte val;
        CDPAtomEnhStereoType(byte v)
        {
            super(0x0446);
            val = v;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(1);
            o.writeByte(val);
        }
    }

    public class CDPShowEnhAtomStereo extends CDPProperty
    {
        byte val;
        CDPShowEnhAtomStereo(boolean v)
        {
            super(0x0445);
            val = v ? (byte)1 : (byte)0;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(1);
            o.writeByte(val);
        }
    }

    public class CDPShowAtomStereo extends CDPProperty
    {
        byte val;
        CDPShowAtomStereo(boolean v)
        {
            super(0x043B);
            val = v ? (byte)1 : (byte)0;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(1);
            o.writeByte(val);
        }
    }

    public class CDPShowBondStereo extends CDPProperty
    {
        byte val;
        CDPShowBondStereo(boolean v)
        {
            super(0x060D);
            val = v ? (byte)1 : (byte)0;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(1);
            o.writeByte(val);
        }
    }

    public class CDPAtomEnhGroup extends CDPProperty
    {
        short val;
        CDPAtomEnhGroup(short g)
        {
            super(0x0447);
            val = g;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(2);
            o.writeShort(val);
        }
    }

    public class CDPPoint2D extends CDPProperty
    {
        int x,y;
        public CDPPoint2D(int x,int y)
        {
            super(0x0200);
            this.x = x;
            this.y = y;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(8);
            o.writeInt(y);
            o.writeInt(x);
        }
        /*
                public String encode()
                {
                    return String.format("%02X %02X 08 00", (type & 0x00FF), (type >> 8)) +   // type and length
                    String.format(" %02X %02X %02X %02X",(y & 0x00FF),(y >> 8),(y >> 16),(y >> 24)) +
                    String.format(" %02X %02X %02X %02X",(x & 0x00FF),(x >> 8),(x >> 16),(x >> 24));
                }
         */
    }

    public class CDPBondType extends CDPProperty
    {
        short btype;
        CDPBondType(int bond)
        {
            super(0x0600);
            this.btype = (short)0xFFFF;
            switch(bond){
                case Molecule.cBondTypeSingle: //			= 0x0001;
                    btype = 1;
                    break;
                case Molecule.cBondTypeDouble: //			= 0x0002;
                    btype = 2;
                    break;
                case Molecule.cBondTypeTriple: //			= 0x0004;
                    btype = 4;
                    break;
                case Molecule.cBondTypeDown: //			= 0x0009;
                    btype = 1;
                    break;
                case Molecule.cBondTypeUp: //				= 0x0011;
                    btype = 1;
                    break;
                case Molecule.cBondTypeCross: //			= 0x001A;
                    btype = 2;
                    break;
                case Molecule.cBondTypeMetalLigand: //    = 0x0020;
                    btype = 1;
                    break;
                case Molecule.cBondTypeDelocalized: //	= 0x0040;
                    btype = 0x0080;
                    break;
                case Molecule.cBondTypeDeleted: //		= 0x0080;
                case Molecule.cBondTypeIncreaseOrder: //  = 0x007F;
                    break;
            }

        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(2);
            o.writeShort(btype);
        }

        // There's only one of these Property in a set!
        public boolean equals(Object o)
        {
            if(o instanceof CDPBondType){
                return true;
            } else{
                return false;
            }
        }

        public int hascode()
        {
            return type;
        }
    }

    public class CDPElement extends CDPProperty
    {
        short ord;
        public CDPElement(short ord)
        {
            super(0x0402);
            this.ord = ord;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(2);
            o.writeShort(ord);
        }
        /*
                public String encode()
                {
                    return String.format("%02X %02X 08 00", (type & 0x00FF), (type >> 8)) +   // type and length
                    String.format(" %02X %02X %02X %02X",(y & 0x00FF),(y >> 8),(y >> 16),(y >> 24)) +
                    String.format(" %02X %02X %02X %02X",(x & 0x00FF),(x >> 8),(x >> 16),(x >> 24));
                }
         */
    }

    public abstract class CDPBondAttach extends CDPProperty
    {
        int ref;
        public CDPBondAttach(int t,int ref)
        {
            super(t);
            this.ref = ref;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(4);
            o.writeInt(ref);
        }

        /*		
                public String encode()
                {
                    return String.format("%02X %02X 04 00", (type & 0x00FF), (type >> 8)) +  // type and length
                    String.format(" %02X %02X %02X %02X", (ref & 0x00FF), (ref >> 8), (ref >> 16), (ref >> 24));
                }
         */
    }

    public class CDPBondBegin extends CDPBondAttach
    {
        public CDPBondBegin(int ref)
        {
            super(0x0604,ref);
        }
    }

    public class CDPBondEnd extends CDPBondAttach
    {
        public CDPBondEnd(int ref)
        {
            super(0x0605,ref);
        }

    }

    public class CDPBondDisplay extends CDPProperty
    {
        short val;
        CDPBondDisplay(short val)
        {
            super(0x0601);
            this.val = val;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(2);
            o.writeShort(val);
        }
    }


    public class CDPReactants extends CDPProperty
    {
        int elems[];
        CDPReactants(int elems[])
        {
            super(0x0C01);
            this.elems = elems;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort(elems.length*4);
            for (int e : elems) {
                o.writeInt(e);
            }
        }
    }

    public class CDPProducts extends CDPProperty
     {
         int elems[];
         CDPProducts(int elems[])
         {
             super(0x0C02);
             this.elems = elems;
         }
 
         public void write(DataOutput o) throws IOException
         {
             o.writeShort(type);
             o.writeShort(elems.length*4);
             for (int e : elems) {
                 o.writeInt(e);
             }
         }
    }

    public class CDPText extends CDPProperty
    {
        String s;
        CDPText(String s)
        {
            super(0x0700);
            this.s = s;
        }

        public void write(DataOutput o) throws IOException
        {
            o.writeShort(type);
            o.writeShort((short)(2 + s.length()));
            o.writeShort((short)0);
            o.write(s.getBytes());
        }
    }

    /*	
        private class CDXObject
        {
            short type;
            java.util.List<CDPProperty> props = null;

            public CDXObject(short t)
            {
                type = t;
            }

            public void addProperty(CDPProperty p)
            {
                if (props == null)
                    props = new java.util.ArrayList<CDPProperty>();
                props.add(p);
            }
        }
     */
    private interface Functor
    {
        void start(CDXNode o);

        void end(CDXNode o);
    }

    /*
        private class Tree
        {
            CDXNode data;
            java.util.List<Tree> children;

            public Tree(CDXNode n)
            {
                data = n;
            }

            public void walk(Functor functor)
            {
                functor.start(data);
                if (children != null) {
                    for (Tree n : children)
                        n.walk(functor);
                }
                functor.end(data);
            }
            public void add(Tree e)
            {
                if (children == null)
                    children = new java.util.ArrayList<Tree>();
                children.add(e);
            }

            public void add(CDXNode node)
            {
                if (children == null)
                    children = new java.util.ArrayList<Tree>();
                children.add(new Tree(node));
            }

        }
     */
    public void dump(CDXNode tn)
    {

        tn.walk(new Functor()
        {
            public void start(CDXNode n)
            {
                System.out.println(n.start());
            }

            public void end(CDXNode n)
            {
                System.out.println(n.end());
            }
        });
    }

    public void write(final DataOutput d,CDXNode tn)
    {
        tn.walk(new Functor()
        {
            public void start(CDXNode n)
            {
                try{
                    n.startWrite(d);
                } catch(IOException e){
                    System.err.println("Error " + e);
                }
            }

            public void end(CDXNode n)
            {
                try{
                    n.endWrite(d);
                } catch(IOException e){
                    System.err.println("Error " + e);
                }
            }
        });

    }

    public void test()
    {
        CDXDocument doc = new CDXDocument();
        CDXNode page = new CDXPage();
        CDXNode frag = new CDXFragment();
//		CDXAtom at1 = new CDXAtom((short)7, 0x00948000, 0x00A44000);
//		CDXAtom at2 = new CDXAtom((short)8, 0x00B17A4E, 0x00AC03BA);
        CDXAtom at1 = new CDXAtom((short)7,0x00,0x00);
        CDXAtom at2 = new CDXAtom((short)8,0x02000000,0x00AC03BA);
        CDXBond b1 = new CDXBond(at1.getID(),at2.getID());
        doc.add(page);
        page.add(frag);
        frag.add(at1);
        frag.add(at2);
        frag.add(b1);
        dump(doc);

        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        LittleEndianDataOutputStream l = new LittleEndianDataOutputStream(bos);
        write(l,doc);
        try{
            bos.close();
        } catch(IOException e){}
        byte buffer[] = bos.toByteArray();
        StringBuilder sb = new StringBuilder();
        for(int i = 0;i < buffer.length;i++){
            byte b = buffer[i];
            if(i != 0){
                sb.append(" ");
            }
            if(i > 0 && (i % 16) == 0){
                sb.append("\n");
            }
            sb.append(String.format("%02X",(int)(b & 0xFF)));
        }
        System.out.println("\nByte Dump:\n" + sb);

        NativeClipboardAccessor.setClipBoardData("ChemDraw Interchange Format",buffer);

        try{
            BufferedOutputStream os = new BufferedOutputStream(new FileOutputStream("d:\\dev\\console\\test.cdx"));
            l = new LittleEndianDataOutputStream(os);
            write(l,doc);
            os.close();
        } catch(IOException e){}

    }
    /*
    public static void main(String args[])
    {
        ChemDrawCDX cdx = new ChemDrawCDX();
//        cdx.test();

        String FILES[] = {
            "out1.mol",
            "1.mol",
            "2.mol",
            "3.mol",
            "4.mol",
            "5.mol",
            "6.mol",
            "7.mol",
        };

        String PATH = "D:\\Dev\\Java\\ClipboardDLL\\tests\\";
        String filename = PATH + FILES[7];
        MolfileParser p = new MolfileParser();
        StereoMolecule m = new StereoMolecule();
        if(p.parse(m,new java.io.File(filename))){
            byte array[] = cdx.getChemDrawBuffer(m);
            NativeClipboardAccessor.setClipBoardData("ChemDraw Interchange Format",array);
        } else{
            System.err.println("Error parsing molfile ");
        }

    }
*/
}
