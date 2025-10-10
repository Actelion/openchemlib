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
import com.actelion.research.chem.reaction.ReactionEncoder;
import com.actelion.research.gui.clipboard.external.ChemDrawCDX;
import com.actelion.research.gui.dnd.MoleculeTransferable;
import com.actelion.research.gui.dnd.ReactionTransferable;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.gui.wmf.WMF;
import com.actelion.research.gui.wmf.WMFGraphics2D;
import com.actelion.research.util.Platform;
import com.actelion.research.util.Sketch;

import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.io.*;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;

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
	public static final String NC_SERIALIZEMOLECULE = "ACT_MOLECULE";
	public static final String NC_SERIALIZEREACTION = "ACT_REACTION";
	public static final String NC_CTAB = "MDLCT";
	public static final String NC_MOLFILE = "MDL_MOL";
	public static final String NC_SKETCH = "MDLSK";
	public static final String NC_EMBEDDEDSKETCH = "MDLSK_EMBEDDED";
	public static final String NC_IDCODE = "IDCODE";
	public static final String NC_CHEMDRAWINTERCHANGE = "ChemDraw Interchange Format";
	public static final String NC_METAFILE = "CF_METAFILEPICT";
	public static final String NC_DIB = "CF_DIB";
	private static final byte MDLSK[] = {(byte) 'M', (byte) 'D', (byte) 'L', (byte) 'S', (byte) 'K', 0, 0};
	public static final java.util.List<String> readableReactionFormats = Arrays.asList(NC_SERIALIZEREACTION, NC_CTAB, NC_IDCODE);
	public static final java.util.List<String> readableMoleculeFormats = Arrays.asList(NC_SERIALIZEMOLECULE, NC_CTAB, NC_MOLFILE, NC_SKETCH, NC_EMBEDDEDSKETCH, NC_IDCODE);

	public static java.util.List<String> writableReactionFormats = Arrays.asList(NC_SERIALIZEREACTION, NC_CTAB, NC_CHEMDRAWINTERCHANGE);
	public static java.util.List<String> writableMoleculeFormats;
	private static java.util.List<String> nativeCliphandler;

	static {
		if (Platform.isWindows()) {
			nativeCliphandler = Arrays.asList(
					"com.actelion.research.gui.clipboard.JNAWinClipboardHandler",
					"com.actelion.research.gui.clipboard.NativeClipboardAccessor"
			);
			writableMoleculeFormats = Arrays.asList(NC_METAFILE, NC_SERIALIZEMOLECULE, NC_CHEMDRAWINTERCHANGE);
		} else {
			nativeCliphandler = new ArrayList<>();
			writableMoleculeFormats = Arrays.asList(NC_SERIALIZEMOLECULE, NC_IDCODE);
		}
	}



	/**
	 * Get one or more Molecule(s) from the Clipboard. On all platforms the first choice is a serialized StereoMolecule.
	 * Further supported formats on Windows are: MDLSK,MDLCT,MDL_MOL,CF_ENHMETAFILE with embedded sketch.
	 * If no supported format is found and the clipboard contains text, which can be interpreted as molfile, then the
	 * corresponding molecule is returned. If the clipboard contains one or multiple SMILES, IUPAC name(s) or idcode(s),
	 * then the corresponding molecule(s) is/are returned. These can be multiple if allowMultiple is true.
	 * Otherwise, if a StructureNameResolver is present, then it tries to interpret the name(s) and returns the
	 * corresponding molecules, which may be limited to a certain number.
	 * @return Molecule found or null if no molecule present on the clipboard
	 */
	@Override
	public ArrayList<StereoMolecule> pasteMolecules() {
		return pasteMolecules(true, true, SmilesParser.SMARTS_MODE_GUESS);
	}

	/**
	 * Get one Molecule from the Clipboard. On all platforms the first choice is a serialized StereoMolecule.
	 * Further supported formats on Windows are: MDLSK,MDLCT,MDL_MOL,CF_ENHMETAFILE with embedded sketch.
	 * If no supported format is found and the clipboard contains text, which can be interpreted as molfile,
	 * SMILES (with SMARTS features), IUPAC name or idcode, then the corresponding molecule is returned.
	 * Otherwise, if a StructureNameResolver is present, then it tries to interpret the name.
	 * If the clipboard molecule has 3D coordinates, then new 2D-coords are invented and used instead.
	 * @return Molecule found or null if no molecule present on the clipboard
	 */
	@Override
	public StereoMolecule pasteMolecule() {
		return pasteMolecule(true, SmilesParser.SMARTS_MODE_GUESS);
	}

	@Override
	public StereoMolecule pasteMolecule(boolean prefer2D, int smartsMode) {
		ArrayList<StereoMolecule> molList = pasteMolecules(prefer2D, false, smartsMode);
		return molList.isEmpty() ? null : molList.get(0);
	}

	/**
	 * Get one or more Molecule(s) from the Clipboard. On all platforms the first choice is a serialized StereoMolecule.
	 * Further supported formats on Windows are: MDLSK,MDLCT,MDL_MOL,CF_ENHMETAFILE with embedded sketch.
	 * If no supported format is found and the clipboard contains text, which can be interpreted as molfile, then the
	 * corresponding molecule is returned. If the clipboard contains one or multiple SMILES, IUPAC name(s) or idcode(s),
	 * then the corresponding molecule(s) is/are returned. These can be multiple if allowMultiple is true.
	 * Otherwise, if a StructureNameResolver is present, then it tries to interpret the name(s) and returns the
	 * corresponding molecules, which may be limited to a certain number.
	 * @param prefer2D if true and if the clipboard molecule has 3D coordinates, then new 2D-coords are invented
	 * @param allowMultiple whether multiple molecules may be generated from clipboard text, if no serialized mol or special molecule format present
	 * @param smartsMode SmilesParser.SMARTS_MODE for parsing strings as SMILES, i.e. how SMARTS features are considered
	 * @return list of molecules found or generated from SMILES, names, etc; empty list if no molecule present on the clipboard
	 */
	private ArrayList<StereoMolecule> pasteMolecules(boolean prefer2D, boolean allowMultiple, int smartsMode) {
		ArrayList<StereoMolecule> molList = new ArrayList<>();

		StereoMolecule mol = /*Platform.isWindows()*/ isNativeClipHandler() ? pasteMoleculeNative(prefer2D) : pasteMoleculeLinux();
		if (mol != null)
			molList.add(mol);

		if (molList.isEmpty()) {
			// get StringFlavor from clipboard and try parsing it as idcode, molfile, smiles, or (if NameResolver exists) as name
			Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
			String text = null;
			try {
				text = (String)t.getTransferData(DataFlavor.stringFlavor);
			}
			catch (Exception ioe) {}
			if (text != null) {
				try {
					mol = new MolfileParser().getCompactMolecule(text);
					if (mol != null && mol.getAllAtoms() != 0)
						molList.add(mol);
				}
				catch (Exception e) {}

				ArrayList<String> unresolvedNameList = null;

				if (molList.isEmpty()) {
					boolean isFirstLine = true;
					int column = -1;    // in case we have a TAB-delimited table with '[idcode]' tags, use first tagged column
					BufferedReader reader = new BufferedReader(new StringReader(text));
					try {
						String line = reader.readLine();
						while (line != null) {
							line = line.trim();
							if (isFirstLine) {
								String[] header = line.split("\\t");
								for (int i=0; i<header.length; i++) {
									if (header[i].endsWith("[idcode]")) {
										column = i;
										break;
									}
								}
								isFirstLine = false;
								if (column != -1)
									continue;
							}
							if (line.isEmpty())
								continue;
							try {
								String idcode = line;
								if (column != -1) {
									String[] entry = line.split("\\t");
									if (column < entry.length)
										idcode = entry[column];
								}
								mol = new IDCodeParser(prefer2D).getCompactMolecule(idcode);
							} catch (Exception e) {
								mol = null;
							}

							if (mol == null || mol.getAllAtoms() == 0) {
								mol = new StereoMolecule();
								try {
									new SmilesParser(smartsMode).parse(mol, line);
								}
								catch (Exception e) {
									mol = null;
								}
							}

							if (mol == null || mol.getAllAtoms() == 0)
								mol = StructureNameResolver.resolveLocal(line);

							if ((mol == null || mol.getAllAtoms() == 0) && !allowMultiple)
								mol = StructureNameResolver.resolveRemote(line);

							if (mol != null && mol.getAllAtoms() != 0) {
								molList.add(mol);

								if (!allowMultiple)
									break;
							}
							else if (allowMultiple) {
								if (unresolvedNameList == null)
									unresolvedNameList = new ArrayList<>();
								unresolvedNameList.add(line);
							}

							line = reader.readLine();
						}
					}
					catch (IOException ioe) {}
				}

				if (unresolvedNameList != null && !unresolvedNameList.isEmpty()) {
					String[] idcodes = StructureNameResolver.resolveRemote(unresolvedNameList.toArray(new String[0]));
					for (String idcode:idcodes) {
						try {
							mol = new IDCodeParser(prefer2D).getCompactMolecule(idcode);
							if (mol != null && mol.getAllAtoms() != 0)
								molList.add(mol);
						} catch (Exception e) {}
					}
				}
			}
		}

		if (prefer2D) {
			for (StereoMolecule m:molList) {
				if (m.is3D()) {
					m.ensureHelperArrays(Molecule.cHelperParities);    // to ensure stereo parities
					new CoordinateInventor().invent(m);
				}
			}
		}

//		System.out.println("returned mol(s): " + molList.size());
		return molList;
	}

	private StereoMolecule pasteMoleculeNative(boolean prefer2D) {
		StereoMolecule mol = null;
		Class clipHandlerClz = getNativeClipHandler();
		if (clipHandlerClz != null) {
			for (String f : readableMoleculeFormats) {
				try {
					mol = rawToMol((byte[]) clipHandlerClz.getMethod("getClipboardData", String.class).invoke(clipHandlerClz, f), f, prefer2D);
				} catch (NoSuchMethodException | InvocationTargetException | IllegalAccessException e) {
					e.printStackTrace();
					mol = pasteMoleculeLinux();
					break;
				}
				if (mol != null) break;
			}
		} else {
			mol = pasteMoleculeLinux();
		}
		return mol;
	}

	private static Class getNativeClipHandler() {

		Class clipHandlerClz = null;

		for (String clz : nativeCliphandler) {
			try {
				clipHandlerClz = Class.forName(clz);
			} catch (ClassNotFoundException | NoClassDefFoundError | UnsatisfiedLinkError e1) {
				e1.printStackTrace();
			}
			try {
				if (clipHandlerClz != null && Boolean.TRUE.equals(clipHandlerClz.getMethod("isInitOK").invoke(clipHandlerClz))) {
					System.out.println("WinClipHandler class: " + clz);
					break;
				}
			} catch (NoSuchMethodException | SecurityException | IllegalAccessException | InvocationTargetException  e2) {
				e2.printStackTrace();
			}
		}
		return clipHandlerClz;
	}

	public static boolean isNativeClipHandler() {
		return nativeCliphandler != null && nativeCliphandler.size() > 0;
	}

	private StereoMolecule pasteMoleculeLinux() {
		try {
			Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
			return (StereoMolecule)t.getTransferData(ChemistryFlavors.DF_SERIALIZED_MOLECULE);
		} catch (Exception e) {}
		return null;
	}

	/**
	 * Get a Reaction from the Clipboard
	 *
	 * @return Reaction or null if no reaction present
	 */
	public Reaction pasteReaction() {
		Reaction rxn = isNativeClipHandler() ? pasteReactionWindowsNative() : pasteReactionLinux();

		if (rxn == null) {
			// get StringFlavor from clipboard and try parsing it as rxn-idcode, rxn-file, or reaction-smiles
			Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
			String text = null;
			try {
				text = (String)t.getTransferData(DataFlavor.stringFlavor);
			}
			catch (Exception ioe) {}
			if (text != null) {
				try {
					rxn = ReactionEncoder.decode(text, true);
					if (rxn != null && rxn.isEmpty())
						rxn = null;
				}
				catch (Exception e) {}
				if (rxn == null) {
					try {
						rxn = new RXNFileParser().getReaction(text);
						if (rxn != null && rxn.isEmpty())
							rxn = null;
					}
					catch (Exception e) {}
				}
				if (rxn == null) {
					try {
						rxn = new SmilesParser().parseReaction(text);
						if (rxn != null && rxn.isEmpty())
							rxn = null;
					}
					catch (Exception e) {}
				}
			}
		}

		return rxn;
	}
	public Reaction pasteReactionWindowsNative() {
		Reaction rxn = null;
		Class clipHandlerClz = getNativeClipHandler();
		if (clipHandlerClz != null) {
			for (String f : readableReactionFormats) {
				try {
					rxn = rawToRxn((byte[]) clipHandlerClz.getMethod("getClipboardData", String.class).invoke(clipHandlerClz, f), f);
				} catch (NoSuchMethodException | InvocationTargetException | IllegalAccessException e) {
					e.printStackTrace();
					rxn = pasteReactionLinux();
					break;
				}
				if (rxn != null) break;
			}
		} else {
			rxn = pasteReactionLinux();
		}
		return rxn;
	}
	private Reaction pasteReactionLinux() {
		try {
			Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
			return (Reaction)t.getTransferData(ChemistryFlavors.DF_SERIALIZED_REACTION);
		} catch (Exception e) {}
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
		if (isNativeClipHandler()) {
			//return copyMoleculeNative(mol);
			Class clipHandlerClz;
			if ((clipHandlerClz = getNativeClipHandler())!= null) {
				try {
					Method method = clipHandlerClz.getMethod("copyMoleculeToClipboard", String.class, byte[].class, byte[].class);

					StereoMolecule m = mol.getCompactCopy();
					for (int atom = 0; atom < m.getAllAtoms(); atom++)
						m.setAtomMapNo(atom, 0, false);

					byte buffer[] = Sketch.createSketchFromMol(m);

					File temp = File.createTempFile("actnca", ".wmf");
					temp.deleteOnExit();

					String path = null;
					if (writeMol2Metafile(temp, m, buffer))
						path = temp.getAbsolutePath();

					// Serialize to a byte array
					/*System.out.println("CopyMolecule");
					com.actelion.research.gui.clipboard.external.ChemDrawCDX cdx = new com.actelion.research.gui.clipboard.external.ChemDrawCDX();
					byte[] cdbuffer = cdx.getChemDrawBuffer(m);

					ByteArrayOutputStream bos = new ByteArrayOutputStream();
					ObjectOutputStream out = new ObjectOutputStream(bos);
					// Changed from m to mol, because writeMol2Metafile() may have scaled xy-coords of m,
					// which is unacceptable for 3D molecules.
					// If an application needs coordinate scaling, then this should be done after pasting. TLS 07Feb2016
					out.writeObject(mol);
					out.close();
					bos.close();*/

					boolean ok = (boolean) method.invoke(clipHandlerClz, path, molToRaw(mol, NC_CHEMDRAWINTERCHANGE), molToRaw(mol,NC_SERIALIZEMOLECULE));

					temp.delete();

					return ok;
				} catch (NoSuchMethodException | IOException | IllegalAccessException | InvocationTargetException e) {
					e.printStackTrace();
				}
			}
		}
		MoleculeTransferable transferable = new MoleculeTransferable(mol);
		Toolkit.getDefaultToolkit().getSystemClipboard().setContents(transferable, transferable);
		return true;
	}

	public static void emptyClipboard() {
		// This Transferable clears the clipboard
		Toolkit.getDefaultToolkit().getSystemClipboard().setContents(new Transferable() {
			public DataFlavor[] getTransferDataFlavors() {
				return new DataFlavor[0];
			}

			public boolean isDataFlavorSupported(DataFlavor flavor) {
				return false;
			}

			public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException {
				throw new UnsupportedFlavorException(flavor);
			}
		}, null);
	}
	public boolean copyMoleculeNative(StereoMolecule mol) {
		Class winClipHandler;
		boolean ok = true;
		if ((winClipHandler = getNativeClipHandler()) != null) {

			emptyClipboard();
			for (String f : writableMoleculeFormats) {
				try{
					byte[] bytes = molToRaw(mol, f);
					if (bytes != null) {
						ok &= (boolean)winClipHandler.getMethod("setClipBoardData", String.class, byte[].class, boolean.class).invoke(winClipHandler, f, bytes, false);
					}
				} catch (InvocationTargetException |IllegalAccessException | NoSuchMethodException e) {
    				e.printStackTrace();
					ok = false;
                }
            }
		} else ok = false;
		return ok;
	}

	public static StereoMolecule rawToMol(byte[] buffer, String nativeFormat, boolean prefer2D) {

		StereoMolecule mol = null;

		if (buffer != null) {
			switch (nativeFormat) {
				case NC_SERIALIZEMOLECULE:
					try (ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(buffer))) {
						Object o = is.readObject();
						is.close();
						if (o instanceof StereoMolecule) {
							mol = (StereoMolecule) o;
						}
					} catch (IOException e) {
						e.printStackTrace();
					} catch (ClassNotFoundException e) {
						e.printStackTrace();
					}
					break;
				case NC_CTAB:
				case NC_MOLFILE:
					MolfileParser p = new MolfileParser();
					mol = new StereoMolecule();
					if (!p.parse(mol, new String(buffer))) {
						mol = null;
						System.out.println("Error Parsing CTAB during clipboard paste");
					}
					break;
				case NC_SKETCH:
				case NC_EMBEDDEDSKETCH:
					try {
						mol = new StereoMolecule();
						if (!Sketch.createMolFromSketchBuffer(mol, buffer)) {
							mol = null;
						}
					} catch (IOException ex) {
						mol = null;
						ex.printStackTrace();
						System.out.println("Parsing NC_SKETCH during clipboard paste: Exception " + ex);
					}
					break;
				case NC_IDCODE:
					String clipboardText = new String(buffer);

					try {
						mol = new StereoMolecule();
						IDCodeParser parser = new IDCodeParser(prefer2D);
						parser.parse(mol, buffer);
						if (mol.getAllAtoms() == 0) {
							mol = null;
						} else {
							System.out.printf("NC_IDCODE '%s' successfully interpreted as idcode\n", clipboardText);
						}
					} catch (Exception ex) {
						ex.printStackTrace();
						System.out.printf("NC_IDCODE '%s' could not be parsed: " + ex + "\n", clipboardText);
						mol = null;
					}
					break;
				default:
					System.out.printf("No known format provided: ");
					break;

			}
		}
		return mol;
	}

	public static Reaction rawToRxn(byte[] buffer, String nativeFormat) {
		Reaction rxn = null;
		if (buffer != null) {
			switch (nativeFormat) {
				case NC_SERIALIZEREACTION:
					try {
						ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(buffer));
						Object o = is.readObject();
						if (o instanceof Reaction) {
							rxn = (Reaction)o;
						}
						is.close();
					} catch (Exception var7) {
						var7.printStackTrace();
						System.out.println("ClipboardHandler.pasteReaction(): Exception " + var7);
					}
					break;
				case NC_CTAB:
					RXNFileParser p = new RXNFileParser();
					rxn = new Reaction();

					try {
						if (!p.parse(rxn, new String(buffer))) {
							rxn = null;
						}
					} catch (Exception var6) {
						System.err.println("Error parsing Reaction Buffer " + var6);
						rxn = null;
					}
					break;
				case NC_SKETCH:
					try {
						rxn = new Reaction();
						if (!Sketch.createReactionFromSketchBuffer(rxn, buffer)) {
							rxn = null;
						}
					} catch (IOException var5) {
						rxn = null;
					}
					break;
				case NC_IDCODE:
					String clipboardText = new String(buffer);
					rxn = ReactionEncoder.decode(clipboardText, true);
					break;
				default:
					System.out.printf("No known format provided: " + nativeFormat);
					break;
			}
		}
		return rxn;
	}
	public byte[] molToRaw(StereoMolecule mol, String nativeFormat) {
		byte[] bytes = null;
		if (mol != null) {
			switch (nativeFormat) {
				case NC_SERIALIZEMOLECULE:
					try {
						ByteArrayOutputStream bos = new ByteArrayOutputStream();
						ObjectOutputStream out = new ObjectOutputStream(bos);
						// Changed from m to mol, because writeMol2Metafile() may have scaled xy-coords of m,
						// which is unacceptable for 3D molecules.
						// If an application needs coordinate scaling, then this should be done after pasting. TLS 07Feb2016
						out.writeObject(mol);
						out.close();
						bos.close();
						bytes = bos.toByteArray();
					} catch (IOException e) {
						e.printStackTrace();
					}
					break;
				case NC_CHEMDRAWINTERCHANGE:
					StereoMolecule m = mol.getCompactCopy();
					for (int atom = 0; atom < m.getAllAtoms(); atom++)
						m.setAtomMapNo(atom, 0, false);
					com.actelion.research.gui.clipboard.external.ChemDrawCDX cdx = new com.actelion.research.gui.clipboard.external.ChemDrawCDX();
					bytes = cdx.getChemDrawBuffer(m);
					break;
				case NC_METAFILE:
					try {
						m = mol.getCompactCopy();
						for (int atom = 0; atom < m.getAllAtoms(); atom++)
							m.setAtomMapNo(atom, 0, false);
						byte buffer[] = Sketch.createSketchFromMol(m);

						File temp = File.createTempFile("actnca", ".wmf");
						temp.deleteOnExit();

						String path = null;
						if (writeMol2Metafile(temp, m, buffer))
							bytes = Files.readAllBytes(temp.toPath());
						temp.delete();
					} catch (IOException e) {
						e.printStackTrace();
					}
					break;
				case NC_IDCODE:
					bytes = mol.getIDCode().getBytes();
				default:
					System.err.println("Unkown Format " + nativeFormat);
			}
		}
		return bytes;
	}

	public byte[] rxnToRaw(Reaction rxn, String ctab, String nativeFormat) {
		byte[] bytes = null;
		if (rxn != null) {
			switch (nativeFormat) {
				case NC_SERIALIZEREACTION:
					try {
						ByteArrayOutputStream bos = new ByteArrayOutputStream();
						ObjectOutputStream out = new ObjectOutputStream(bos);
						out.writeObject(rxn);
						out.close();
						bos.close();
						bytes = bos.toByteArray();
					} catch (IOException e) {
						e.printStackTrace();
					}
					break;
				case NC_CHEMDRAWINTERCHANGE:
					ChemDrawCDX cdx = new ChemDrawCDX();
					bytes = cdx.getChemDrawBuffer(rxn);
					break;
				case NC_CTAB:
					if (ctab == null) {
						RXNFileCreator mc = new RXNFileCreator(rxn);
						ctab = mc.getRXNfile();
					}
					bytes = ctab.getBytes();
					break;
				case NC_IDCODE:
					bytes = ReactionEncoder.encode(rxn, true, ReactionEncoder.INCLUDE_DEFAULT).getBytes();
				default:
					System.err.println("Unkown Format " + nativeFormat);
			}
		}
		return bytes;
	}

	public boolean copyReaction(Reaction r) {
		return copyReactionToClipboard(null, r);
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
		} catch (Exception e) {
			System.err.println("ClipboardHandler: Exception copying reaction " + e);
		}
		return ok;
	}

	/**
	 * Copies a reaction to the clipboard in various formats:
	 * MDLSK Sketch
	 * MDLCT MDL molfile
	 */
	public boolean copyReactionToClipboard(String ctab, Reaction rxn) {
		if (isNativeClipHandler()) {
			//return copyReactionNative(ctab, rxn);
			Class clipHandlerClz;
			if ((clipHandlerClz = getNativeClipHandler())!= null) {
				try {
					Method method = clipHandlerClz.getMethod("copyReactionToClipboard", byte[].class, byte[].class, byte[].class);
					/*if (ctab == null) {
						RXNFileCreator mc = new RXNFileCreator(rxn);
						ctab = mc.getRXNfile();
					}
					ChemDrawCDX cdx = new ChemDrawCDX();
					byte[] cdxBuffer = cdx.getChemDrawBuffer(rxn);
					System.out.println("copyReactionToClipboard");
					ByteArrayOutputStream bos = new ByteArrayOutputStream();
					ObjectOutputStream out = new ObjectOutputStream(bos);
					out.writeObject(rxn);
					out.close();
					bos.close();*/

					return (boolean) method.invoke(clipHandlerClz, rxnToRaw(rxn, ctab, NC_CTAB), rxnToRaw(rxn, ctab, NC_CHEMDRAWINTERCHANGE), rxnToRaw(rxn, ctab, NC_SERIALIZEREACTION));

				} catch (NoSuchMethodException | IllegalAccessException | InvocationTargetException e) {
					e.printStackTrace();
				}
			}
		}

		ReactionTransferable transferable = new ReactionTransferable(rxn);
		Toolkit.getDefaultToolkit().getSystemClipboard().setContents(transferable, transferable);
		return true;

	}
	public boolean copyReactionNative(String ctab, Reaction rxn) {
		Class clipHandlerClz = null;
		boolean ok = true;
		if ((clipHandlerClz = getNativeClipHandler()) != null) {
			emptyClipboard();
			for (String f : writableReactionFormats) {
				try{
					byte[] bytes = rxnToRaw(rxn, ctab, f);
					if (bytes != null) {
						ok &= (boolean)clipHandlerClz.getMethod("setClipBoardData", String.class, byte[].class, boolean.class).invoke(clipHandlerClz, f, bytes, false);
					}
				} catch (InvocationTargetException |IllegalAccessException | NoSuchMethodException e) {
					e.printStackTrace();
					ok = false;
				}
			}
		} else ok = false;
		return ok;
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

		Depictor2D d = new Depictor2D(m);
		d.updateCoords(g, new GenericRectangle(0, 0, w, h), AbstractDepictor.cModeInflateToMaxAVBL);
		d.paint(g);

		if (sketch != null) {
			byte[] temp = new byte[MDLSK.length + sketch.length];
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
		if (isNativeClipHandler()) {
			Class clipHandlerClz = getNativeClipHandler();
			if (clipHandlerClz != null) {
				try {
					return (boolean) clipHandlerClz.getMethod("setClipBoardData", String.class, byte[].class).invoke(clipHandlerClz, format, buffer);
				} catch (NoSuchMethodException | InvocationTargetException | IllegalAccessException e) {
					e.printStackTrace();
				}
			}
		}
		return false;
	}

	/**
	 * Copy a windows enhance metafile to the Windows clipboard
	 * @param data byte[]
	 * @return boolean
	 */
	public static boolean copyMetaFile(byte []data) {
		return Platform.isWindows() && setClipBoardData(NC_METAFILE, data);
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

	public static void setNativeCliphandler(java.util.List<String> nativeCliphandler) {
		ClipboardHandler.nativeCliphandler = nativeCliphandler;
	}

	public static java.util.List<String> getNativeCliphandler() {
		return ClipboardHandler.nativeCliphandler;
	}
}
