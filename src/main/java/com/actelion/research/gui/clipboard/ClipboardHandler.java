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
import com.actelion.research.chem.dnd.ChemistryFlavors;
import com.actelion.research.chem.io.RXNFileCreator;
import com.actelion.research.chem.io.RXNFileParser;
import com.actelion.research.chem.name.StructureNameResolver;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.gui.dnd.MoleculeTransferable;
import com.actelion.research.gui.dnd.ReactionTransferable;
import com.actelion.research.gui.wmf.WMF;
import com.actelion.research.gui.wmf.WMFGraphics2D;
import com.actelion.research.util.Platform;
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
	public StereoMolecule pasteMolecule() {
		return pasteMolecule(true);
	}

		/**
		 * Get a Molecule from the Clipboard. The supported formats are: MDLSK,MDLCT,MDL_MOL,CF_ENHMETAFILE with embedded sketch
		 *
		 * @param prefer2D if true and if the clipboard molecule has 3D coordinates, then new 2D-coords are invented
		 * @return Molecule found or null if no molecule present on the clipboard
		 */
	public StereoMolecule pasteMolecule(boolean prefer2D) {
		StereoMolecule mol;

		if (Platform.isWindows())
			mol = pasteMoleculeWindowsNative(prefer2D);
		else
			mol = pasteMoleculeLinux();

		if (mol == null) {
			// get StringFlavor from clipboard and try parsing it as idcode, molfile, smiles, or (if NameResolver exists) as name
			Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
			String text = null;
			try {
                text = (String)t.getTransferData(DataFlavor.stringFlavor);
				}
			catch (Exception ioe) {}
			if (text != null) {
				try { mol = new MolfileParser().getCompactMolecule(text); } catch (Exception e) {}
				if (mol == null)
					try { mol = new IDCodeParser(prefer2D).getCompactMolecule(text); } catch (Exception e) {}
				if (mol == null) {
					mol = new StereoMolecule();
					try {
						new SmilesParser().parse(mol, text);
					}
					catch (Exception e) {
						mol = null;
					}
				}

				if (mol == null)
					mol = StructureNameResolver.resolve(text);
			}
		}
		if (prefer2D && mol != null && mol.is3D()) {
			mol.ensureHelperArrays(Molecule.cHelperParities);    // to ensure stereo parities
			new CoordinateInventor().invent(mol);
		}

		System.out.println("returned Mol is " + mol);
		return mol;
	}

	private StereoMolecule pasteMoleculeWindowsNative(boolean prefer2D) {
		byte[] buffer;
		StereoMolecule mol = null;

		if ((buffer = WindowsClipboardAccessor.getClipboardData(WindowsClipboardAccessor.NC_SERIALIZEMOLECULE)) != null) {
			try {
				ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(buffer));
				Object o = is.readObject();
				if (o instanceof StereoMolecule) {
					mol = (StereoMolecule) o;
				}
				is.close();
			} catch (Exception e) {
				e.printStackTrace();
				System.out.println("Parsing NC_SERIALIZEMOLECULE during clipboard paste: Exception " + e);
			}
		}
		System.out.println("Mol is " + mol);
		if (mol == null) {
			if ((buffer = WindowsClipboardAccessor.getClipboardData(WindowsClipboardAccessor.NC_CTAB)) != null
					|| (buffer = WindowsClipboardAccessor.getClipboardData(WindowsClipboardAccessor.NC_MOLFILE)) != null) {
				MolfileParser p = new MolfileParser();
				mol = new StereoMolecule();
				if (!p.parse(mol, new String(buffer))) {
					mol = null;
					System.out.println("Error Parsing CTAB during clipboard paste");
				}
			}
			if (mol == null) {
				if ((buffer = WindowsClipboardAccessor.getClipboardData(WindowsClipboardAccessor.NC_SKETCH)) != null
						|| (buffer = WindowsClipboardAccessor.getClipboardData(WindowsClipboardAccessor.NC_EMBEDDEDSKETCH)) != null) {
					try {
						mol = new StereoMolecule();
						if (!Sketch.createMolFromSketchBuffer(mol, buffer)) {
							mol = null;
						}
					} catch (IOException e) {
						mol = null;
						e.printStackTrace();
						System.out.println("Parsing NC_SKETCH during clipboard paste: Exception " + e);
					}
				}
			}
		}
		String clipboardText = null;
		if (mol == null) {
			if ((buffer = WindowsClipboardAccessor.getClipboardData(WindowsClipboardAccessor.NC_IDCODE)) != null) {
				clipboardText = new String(buffer);
				try {
					mol = new StereoMolecule();
					IDCodeParser parser = new IDCodeParser(prefer2D);
					parser.parse(mol,buffer);
					if (mol.getAllAtoms() == 0)
						mol = null;
					else
						System.out.printf("NC_IDCODE '%s' successfully interpreted as idcode\n",clipboardText);
				} catch (Exception e) {
					e.printStackTrace();
					System.out.printf("NC_IDCODE '%s' could not be parsed: "+e+"\n", clipboardText);
					mol = null;
				}
			}
		}

		return mol;
	}

	private StereoMolecule pasteMoleculeLinux() {
		try {
			Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
			return (StereoMolecule)t.getTransferData(ChemistryFlavors.DF_SERIALIZED_MOLECULE);
		} catch (Exception e) {
			System.err.println("error getting clipboard data "+ e);
		}
	return null;
	}

	/**
	 * Get a Reaction from the Clipboard
	 *
	 * @return Reaction or null if no reaction present
	 */
	public Reaction pasteReaction() {
		Reaction rxn;

		return Platform.isWindows() ? pasteReactionWindowsNative() : pasteReactionLinux();
	}

	public Reaction pasteReactionWindowsNative() {
		byte[] buffer;
		Reaction rxn = null;

		if ((buffer = WindowsClipboardAccessor.getClipboardData(WindowsClipboardAccessor.NC_SERIALIZEREACTION)) != null) {
			try {
				ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(buffer));
				Object o = is.readObject();
				if (o instanceof Reaction) {
					rxn = (Reaction) o;
				}
				is.close();
			} catch (Exception e) {
				e.printStackTrace();
				System.out.println("ClipboardHandler.pasteReaction(): Exception " + e);
			}
		} else if ((buffer = WindowsClipboardAccessor.getClipboardData(WindowsClipboardAccessor.NC_CTAB)) != null) {
			RXNFileParser p = new RXNFileParser();
			rxn = new Reaction();
			try {
				if (!p.parse(rxn, new String(buffer)))
					rxn = null;
			} catch (Exception e) {
				System.err.println("Error parsing Reaction Buffer " + e);
				rxn = null;
			}
		} else if ((buffer = WindowsClipboardAccessor.getClipboardData(WindowsClipboardAccessor.NC_SKETCH)) != null) {
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

	private Reaction pasteReactionLinux() {
		try {
			Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
			return (Reaction)t.getTransferData(ChemistryFlavors.DF_SERIALIZED_REACTION);
		} catch (Exception e) {
			System.err.println("error getting clipboard data "+ e);
		}
		return null;
	}

	public boolean copyMolecule(String molfile) {
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
	public boolean copyMolecule(StereoMolecule mol) {
		if (!Platform.isWindows()) {
			MoleculeTransferable transferable = new MoleculeTransferable(mol);
			Toolkit.getDefaultToolkit().getSystemClipboard().setContents(transferable, transferable);
			return true;
			}

		// For now we keep the old handling for Windows...
		boolean ok = false;
		try {
			StereoMolecule m = mol.getCompactCopy();
			for (int atom=0; atom<m.getAllAtoms(); atom++)
				m.setAtomMapNo(atom, 0, false);

			byte buffer[] = Sketch.createSketchFromMol(m);

			File temp = File.createTempFile("actnca", ".wmf");
			temp.deleteOnExit();

			String path = null;
			if (writeMol2Metafile(temp, m, buffer))
				path = temp.getAbsolutePath();

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

			ok = WindowsClipboardAccessor.copyMoleculeToClipboard(path, cdbuffer, bos.toByteArray());

			temp.delete();
		} catch (IOException e) {
			System.err.println("ClipboardHandler: Exception copying Molecule " + e);
		}
		return ok;
	}


	/**
	 * Copies a reaction to the clipboard in various formats:
	 * MDLSK Sketch
	 * MDLCT MDL molfile
	 */
	public boolean copyReaction(Reaction r) {
		if (!Platform.isWindows()) {
			ReactionTransferable transferable = new ReactionTransferable(r);
			Toolkit.getDefaultToolkit().getSystemClipboard().setContents(transferable, transferable);
			return true;
		}

		// For now we keep the old handling for Windows...
		boolean ok = false;
		try {
			ok = copyReactionToClipboard(null, r);
		} catch (IOException e) {
			System.err.println("ClipboardHandler: Exception Copying Reaction " + e);
		}
		return ok;
	}

	private Reaction makeRXNCopy(Reaction r)  {
		Reaction rxn = new Reaction(r);
		int mols = rxn.getMolecules();
		for (int i = 0; i < mols; i++) {
			rxn.getMolecule(i).ensureHelperArrays(Molecule.cHelperCIP);
		}
		return rxn;
	}

	/**
	 * Copies a reaction to the clipboard in various formats:
	 * CTAB with an embedded sketch
	 * MDLSK Sketch
	 * serialized
	 */
	public boolean copyReaction(String ctab) {
		boolean ok = false;
		try {
			Reaction rxn = new Reaction();
			RXNFileParser p = new RXNFileParser();
			p.parse(rxn, ctab);
			ok = copyReactionToClipboard(ctab, rxn);
		} catch (IOException e) {
			System.err.println("ClipboardHandler: Exception copying reaction " + e);
		} catch (Exception e) {
			System.err.println("ClipboardHandler: Exception copying reaction " + e);
		}
		return ok;
	}

	private boolean copyReactionToClipboard(String ctab, Reaction rxn) throws IOException {
		if (ctab == null) {
			RXNFileCreator mc = new RXNFileCreator(rxn);
			ctab = mc.getRXNfile();
		}

		byte sketch[] = Sketch.createSketchFromReaction(makeRXNCopy(rxn));

		// Serialize to a byte array
		System.out.println("copyReactionToClipboard");

		ByteArrayOutputStream bos = new ByteArrayOutputStream();
		ObjectOutputStream out = new ObjectOutputStream(bos);
		out.writeObject(rxn);
		out.close();
		bos.close();

		return WindowsClipboardAccessor.copyReactionToClipboard(ctab.getBytes(), sketch, bos.toByteArray());
	}

	private boolean writeMol2Metafile(File temp, StereoMolecule m, byte[] sketch) {
		boolean ok = false;
		try {
			ok = writeMol2Metafile(new FileOutputStream(temp), m, sketch);
		} catch (Exception e) {
			System.err.println("ClipboardHandler: Exception writing molfile " + e);
			e.printStackTrace();
		}
		return ok;
	}


	private boolean writeMol2Metafile(OutputStream out, StereoMolecule m, byte[] sketch) throws IOException {
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

	public static boolean setClipBoardData(String format, byte[] buffer) {
		if (Platform.isWindows()) {
			return WindowsClipboardAccessor.setClipBoardData(format,buffer);
		} else
			return false;
	}

	/**
	 * Copy a windows enhance metafile to the Windows clipboard
	 * @param data byte[]
	 * @return boolean
	 */
	public static boolean copyMetaFile(byte []data) {
		return Platform.isWindows() ? setClipBoardData(WindowsClipboardAccessor.NC_METAFILE,data) : false;
		}

/*	private boolean writeRXN2Metafile(File temp, byte sketch[], Reaction m)
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

/*	private boolean writeRXN2Metafile(OutputStream out, byte sketch[], Reaction m) throws IOException
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
	public boolean copyImage(java.awt.Image img) {
		return ImageClipboardHandler.copyImage(img);
	}

	public java.awt.Image pasteImage() {
		return ImageClipboardHandler.pasteImage();
	}

	/**
	 * @deprecated Use ImageClipboardHandler.pasteImage for consistency reasons
	 */
	public static Image getImage() {
		return ImageClipboardHandler.pasteImage();
	}

	/**
	 * @deprecated You may use ImageClipboardHandler.copyImage for consistency reasons
	 */
	public static void putImage(Image image) {
		ImageClipboardHandler.copyImage(image);
	}
}
