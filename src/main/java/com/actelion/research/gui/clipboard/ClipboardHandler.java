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
package com.actelion.research.gui.clipboard;

import com.actelion.research.chem.*;
import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.chem.io.RXNFileCreator;
import com.actelion.research.chem.io.RXNFileParser;
import com.actelion.research.chem.name.StructureNameResolver;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.gui.clipboard.external.ChemDrawCDX;
import com.actelion.research.gui.wmf.WMF;
import com.actelion.research.gui.wmf.WMFGraphics2D;
import com.actelion.research.util.Sketch;

import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.geom.Rectangle2D;
import java.io.*;

/**
 * <p>Title: Actelion Library</p>
 * <p>Description: Actelion Java Library</p>
 * <p>Copyright: Copyright (c) 2002-2003</p>
 * <p>Company: Actelion Ltd</p>
 *
 * @author Thomas Sander, Christian Rufener
 * @version 1.0
 */
public class ClipboardHandler implements IClipboardHandler
{
    private static final byte MDLSK[] = {(byte) 'M', (byte) 'D', (byte) 'L', (byte) 'S', (byte) 'K', 0, 0};

    /**
     * Get a Molecule from the Clipboard. The supported formats are: MDLSK,MDLCT,MDL_MOL,CF_ENHMETAFILE with embedded sketch
     * If the clipboard molecule has 3D coordinates, then new 2D-coords are invented and used instead.
     *
     * @return Molecule found or null if no molecule present on the clipboard
     */
    public StereoMolecule pasteMolecule()
    {
        return pasteMolecule(true);
    }

        /**
		 * Get a Molecule from the Clipboard. The supported formats are: MDLSK,MDLCT,MDL_MOL,CF_ENHMETAFILE with embedded sketch
		 *
         * @param prefer2D if true and if the clipboard molecule has 3D coordinates, then new 2D-coords are invented
		 * @return Molecule found or null if no molecule present on the clipboard
		 */
    public StereoMolecule pasteMolecule(boolean prefer2D)
    {
        byte[] buffer = null;
        StereoMolecule mol = null;

        if ((buffer = NativeClipboardHandler.getClipboardData(NativeClipboardHandler.NC_SERIALIZEMOLECULE)) != null) {
            try {
                ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(buffer));
                Object o = is.readObject();
                System.out.println("Object read from Bytearray input " + o);
                if (o instanceof StereoMolecule) {
                    mol = (StereoMolecule) o;
                }
                is.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.out.println("NativeClipboardAccessor.pasteMolecule(): Exception " + e);
            }
        }
        System.out.println("Mol is " + mol);
        if (mol == null) {
            if ((buffer = NativeClipboardHandler.getClipboardData(NativeClipboardHandler.NC_CTAB)) != null || (buffer = NativeClipboardHandler.getClipboardData(NativeClipboardHandler.NC_MOLFILE)) != null) {
                MolfileParser p = new MolfileParser();
                mol = new StereoMolecule();
                if (!p.parse(mol, new String(buffer))) {
                    mol = null;
                    System.err.println("Error Parsing CTAB during clipboard paste");
                }
            }
            if (mol == null) {
                if ((buffer = NativeClipboardHandler.getClipboardData(NativeClipboardHandler.NC_SKETCH)) != null || (buffer = NativeClipboardHandler.getClipboardData(NativeClipboardHandler.NC_EMBEDDEDSKETCH)) != null) {
                    try {
                        mol = new StereoMolecule();
                        if (!Sketch.createMolFromSketchBuffer(mol, buffer)) {
                            mol = null;
                        }
                    } catch (IOException e) {
                        mol = null;
                        e.printStackTrace();
                        System.out.println("NativeClipboardAccessor.pasteMolecule(): Exception " + e);
                    }
                }
            }
        }
        String clipboardText = null;
        if (mol == null) {
            if ((buffer = NativeClipboardHandler.getClipboardData(NativeClipboardHandler.NC_IDCODE)) != null) {
                clipboardText = new String(buffer);
                try {
                    mol = new StereoMolecule();
                    IDCodeParser parser = new IDCodeParser(prefer2D);
                    System.out.printf("Pasted string '%s'\n",clipboardText);
                    parser.parse(mol,buffer);
                    if (mol.getAllAtoms() == 0)
                    	mol = null;
                } catch (Exception e) {
                    mol = null;
                    System.out.println("NativeClipboardAccessor.pasteMolecule(): Exception " + e);
                }
            }
        }
        if (mol == null) {
            // get StringFlavor from clipboard and try parsing it as molfile, smiles, or (if NameResolver exists) as name
            if (clipboardText == null) {
                Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
                try {
                    Object o = t.getTransferData(DataFlavor.stringFlavor);
                    if (o != null)
                        clipboardText = o.toString();
                }
                catch (Exception ioe) {}
            }
            if (clipboardText != null) {
                mol = new MolfileParser().getCompactMolecule(clipboardText);
                if (mol == null) {
                    mol = new StereoMolecule();
                    try {
                        new SmilesParser().parse(mol, clipboardText);
                    }
                    catch (Exception e) {
                        mol = null;
                    }
                }

                if (mol == null)
                    mol = StructureNameResolver.resolve(clipboardText);
            }
        }
        if (prefer2D && mol != null && is3DMolecule(mol)) {
            mol.ensureHelperArrays(Molecule.cHelperParities);    // to ensure stereo parities
            new CoordinateInventor().invent(mol);
//			mol.setStereoBondsFromParity(); not needed anymore
        }

        System.out.println("returned Mol is " + mol);
        return mol;
    }


    /**
     * Get a Reaction from the Clipboard
     *
     * @return Reaction or null if no reaction present
     */
    public Reaction pasteReaction()
    {
        byte[] buffer = null;
        Reaction rxn = null;


        if ((buffer = NativeClipboardHandler.getClipboardData(NativeClipboardHandler.NC_SERIALIZEREACTION)) != null) {
            try {
                ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(buffer));
                Object o = is.readObject();
                if (o instanceof Reaction) {
                    rxn = (Reaction) o;
                }
                is.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.out.println("NativeClipboardAccessor.pasteMolecule(): Exception " + e);
            }
        } else if ((buffer = NativeClipboardHandler.getClipboardData(NativeClipboardHandler.NC_CTAB)) != null
                || (buffer = NativeClipboardHandler.getClipboardData(NativeClipboardHandler.NC_MOLFILE)) != null) {
            RXNFileParser p = new RXNFileParser();
            rxn = new Reaction();
            try {
                if (!p.parse(rxn, new String(buffer)))
                    rxn = null;
            } catch (Exception e) {
                System.err.println("Error parsing Reaction Buffer " + e);
                rxn = null;
            }
        } else if ((buffer = NativeClipboardHandler.getClipboardData(NativeClipboardHandler.NC_SKETCH)) != null
                || (buffer = NativeClipboardHandler.getClipboardData(NativeClipboardHandler.NC_EMBEDDEDSKETCH)) != null) {
            try {
                rxn = new Reaction();
                if (!Sketch.createReactionFromSketchBuffer(rxn, buffer)) {
                    rxn = null;
                }
            } catch (IOException e) {
                rxn = null;
            }
        }
        return rxn;
    }

    public boolean copyMolecule(String molfile)
    {
        StereoMolecule m = new StereoMolecule();
        MolfileParser p = new MolfileParser();
        p.parse(m, molfile);
        return copyMolecule(m);
    }


    /**
     * Copies a molecule to the clipboard in various formats:
     * ENHMETAFILE with an embedded sketch
     * MDLSK Sketch
     * MDLCT MDL molfile
     */
    public boolean copyMolecule(StereoMolecule mol)
    {
        boolean ok = false;
        try {
            StereoMolecule m = mol.getCompactCopy();
            for (int atom=0; atom<m.getAllAtoms(); atom++)
            	m.setAtomMapNo(atom, 0, false);

            byte buffer[] = Sketch.createSketchFromMol(m);

            File temp = File.createTempFile("actnca", ".wmf");
            temp.deleteOnExit();

            if (writeMol2Metafile(temp, m, buffer)) {
                // Serialize to a byte array
                System.out.println("CopyMolecule");
                com.actelion.research.gui.clipboard.external.ChemDrawCDX cdx = new com.actelion.research.gui.clipboard.external.ChemDrawCDX();
                byte[] cdbuffer = cdx.getChemDrawBuffer(m);

                ByteArrayOutputStream bos = new ByteArrayOutputStream();
                ObjectOutputStream out = new ObjectOutputStream(bos);

	            // Changed from m to mol, because writeMol2Metafile() may have scaled xy-coords of m,
	            // which is unacceptable for 3D molecules.
	            // If an application needs coordinate scaling, then this should be done after pasting. TLS 07Feb2016
                out.writeObject(mol);

                out.close();
                bos.close();
//				ok = NativeClipboardHandler.copyMoleculeToClipboard(temp.getAbsolutePath(),buffer,bos.toByteArray());
                ok = NativeClipboardHandler.copyMoleculeToClipboard(temp.getAbsolutePath(), cdbuffer, bos.toByteArray());
                temp.delete();
            }
        } catch (IOException e) {
            System.err.println("Error copying Molecule " + e);
        }
        return ok;
    }


    /**
     * Copies a molecule to the clipboard in various formats:
     * ENHMETAFILE with an embedded sketch
     * MDLSK Sketch
     * MDLCT MDL molfile
     */
    public boolean copyReaction(Reaction r)
    {
        boolean ok = false;
        try {
            Reaction rxn = makeRXNCopy(r);
            RXNFileCreator mc = new RXNFileCreator(rxn);
            String ctab = mc.getRXNfile();
            ok = copyReactionToClipboard(rxn, ctab);
        } catch (IOException e) {
            System.err.println("Error Copying Reaction " + e);
        }
        return ok;
    }

    private Reaction makeRXNCopy(Reaction r)
    {
        Reaction rxn = new Reaction(r);
        int mols = rxn.getMolecules();
        for (int i = 0; i < mols; i++) {
            rxn.getMolecule(i).ensureHelperArrays(Molecule.cHelperCIP);
        }
        return rxn;
    }
    /**
     * Copies a molecule to the clipboard in various formats:
     * ENHMETAFILE with an embedded sketch
     * MDLSK Sketch
     * MDLCT MDL molfile
     */
    public boolean copyReaction(String ctab)
    {
        boolean ok = false;
        try {
            Reaction rxn = new Reaction();
            RXNFileParser p = new RXNFileParser();
            p.parse(rxn, ctab);
            ok = copyReactionToClipboard(rxn, ctab);
        } catch (IOException e) {
            System.err.println("Error copy reaction " + e);
        } catch (Exception e) {
            System.err.println("Error copy reaction " + e);
        }
        return ok;
    }

    private boolean copyReactionToClipboard(Reaction m, String molfile) throws IOException
    {
        boolean ok = false;
        byte buffer[] = Sketch.createSketchFromReaction(m);
        File temp = File.createTempFile("actnca", ".wmf");
        temp.deleteOnExit();
        // This is currently not used
        //if (writeRXN2Metafile(temp, buffer, m))
        {
            ChemDrawCDX cdx = new com.actelion.research.gui.clipboard.external.ChemDrawCDX();
            byte[] cdbuffer = cdx.getChemDrawBuffer(m);
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream out = new ObjectOutputStream(bos);
            out.writeObject(m);
            out.close();
            bos.close();
            byte t[] = bos.toByteArray();
            System.out.println("Reaction copy with serialized object " + (t != null) + " " + (t != null ? t.length : 0));

//            ok = NativeClipboardHandler.copyReactionToClipboard(temp.getAbsolutePath(),buffer,bos.toByteArray());
            ok = NativeClipboardHandler.copyReactionToClipboard(temp.getAbsolutePath(), cdbuffer, bos.toByteArray());
            temp.delete();
        }
        return ok;
    }

    private boolean writeMol2Metafile(File temp, StereoMolecule m, byte[] sketch)
    {
        boolean ok = false;
        try {
            ok = writeMol2Metafile(new FileOutputStream(temp), m, sketch);
        } catch (Exception e) {
            System.err.println("error writing molfile " + e);
            e.printStackTrace();
        }
        return ok;
    }


    private boolean writeMol2Metafile(OutputStream out, StereoMolecule m, byte[] sketch) throws IOException
    {
        int w = 300;
        int h = 200;
        WMF wmf = new WMF();
        WMFGraphics2D g = new WMFGraphics2D(wmf, w, h, Color.black, Color.white);

        Depictor d = new Depictor(m);
        d.updateCoords(g, new Rectangle2D.Double(0, 0, w, h), AbstractDepictor.cModeInflateToMaxAVBL);
        d.paint(g);

        if (sketch != null) {
            byte temp[] = new byte[MDLSK.length + sketch.length];
            System.arraycopy(MDLSK, 0, temp, 0, MDLSK.length);
            System.arraycopy(sketch, 0, temp, MDLSK.length, sketch.length);
            wmf.escape(WMF.MFCOMMENT, temp);
        }
        wmf.writeWMF(out);
        out.close();
//		g.dispose();
        return true;
    }

/*

    private boolean writeRXN2Metafile(File temp, byte sketch[], Reaction m)
    {
        try {
            return writeRXN2Metafile(new FileOutputStream(temp), sketch, m);
        } catch (Exception e) {
            System.err.println("Error writeRXN2Metafile " + e);
            e.printStackTrace();
            return false;
        }
    }
*/

/*
    private boolean writeRXN2Metafile(OutputStream out, byte sketch[], Reaction m) throws IOException
    {
        int w = 400;
        int h = 300;
        WMF wmf = new WMF();
        WMFGraphics2D g = new WMFGraphics2D(wmf, w, h, Color.black, Color.white);
        AbstractReactionDepictor d = new ReactionDepictor(m);
        d.updateCoords(new GraphicsContext(g), 8, 8, w - 16, h - 16, AbstractDepictor.cModeInflateToMaxAVBL);
        d.paint(new GraphicsContext(g));
        if (sketch != null) {
//			byte MDLSK[] = {(byte)'M',(byte)'D',(byte)'L',(byte)'S',(byte)'K',0,0};
            byte temp[] = new byte[MDLSK.length + sketch.length];
            System.arraycopy(MDLSK, 0, temp, 0, MDLSK.length);
            System.arraycopy(sketch, 0, temp, MDLSK.length, sketch.length);
            wmf.escape(WMF.MFCOMMENT, temp);
        }
        wmf.writeWMF(out);
        out.close();
        return true;
    }
*/


    /**
     * Copies an Image to the clipboard
     *
     * @param img Image to be copied
     * @return true on success
     */
    public boolean copyImage(java.awt.Image img)
    {
        return ImageClipboardHandler.copyImage(img);
    }

    public java.awt.Image pasteImage()
    {
        return ImageClipboardHandler.pasteImage();
    }

    private boolean is3DMolecule(ExtendedMolecule mol)
    {
        for (int atom = 0; atom < mol.getAllAtoms(); atom++)
            if (mol.getAtomZ(atom) != 0.0)
                return true;

        return false;
    }

    /**
     * @deprecated Use ImageClipboardHandler.pasteImage for consistency reasons
     */
    public static Image getImage()
    {
        return ImageClipboardHandler.pasteImage();
    }

    /**
     * @deprecated You may use ImageClipboardHandler.copyImage for consistency reasons
     */
    public static void putImage(Image image)
    {
        ImageClipboardHandler.copyImage(image);
    }
}
