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
import java.util.Iterator;

import static com.actelion.research.gui.clipboard.ClipboardFormat.*;
import static com.actelion.research.gui.clipboard.ClipboardFormat.NC_CHEMDRAWINTERCHANGE;
import static com.actelion.research.gui.clipboard.ClipboardFormat.NC_CHEMDRAWINTERCHANGE;

/**
 * <p>Title: Actelion Library</p>
 * <p>Description: Actelion Java Library</p>
 * <p>Copyright: Copyright (c) 2002-2003</p>
 * <p>Company: Actelion Ltd</p>
 *
 * @author Thomas Sander, Christian Rufener
 * @version 1.0
 */
public class ClipboardHandler implements IClipboardHandler {

    private static final byte MDLSK[] = {(byte) 'M', (byte) 'D', (byte) 'L', (byte) 'S', (byte) 'K', 0, 0};


    private static java.util.List<ClipboardFormat> readableMoleculeFormats =
            Arrays.asList(
                    NC_IDCODE, NC_DATAFLAVOR_IDCODE,NC_SERIALIZEMOLECULE,
                    NC_DATAFLAVOR_SERIALIZEDMOLECULE, NC_MOLFILE_V3, NC_MOLFILE,
                    COM_MDLI_MOLFILE, NC_CTAB,
                    NC_SKETCH, COM_MDLI_SKETCHFILE,
                    //NC_CHEMDRAWINTERCHANGE
                    COM_PE_CDX,
                    NC_EMBEDDEDSKETCH,
                    NC_METAFILE,
                    NC_SMILES);

    private static java.util.List<ClipboardFormat> writableMoleculeFormats =
            Arrays.asList(
                    NC_IDCODE, NC_SERIALIZEMOLECULE, NC_DATAFLAVOR_SERIALIZEDMOLECULE,
                    NC_MOLFILE_V3, NC_MOLFILE, NC_DATAFLAVOR_IDCODE,
                    //NC_CHEMDRAWINTERCHANGE,
                    NC_CTAB,NC_SKETCH, NC_METAFILE
                    ,COM_MDLI_MOLFILE, COM_MDLI_SKETCHFILE, COM_PE_CDX, NC_SMILES
            );

    private static java.util.List<ClipboardFormat> readableReactionFormats =
            Arrays.asList(NC_SERIALIZEREACTION, NC_DATAFLAVOR_SERIALIZEDREACTION,
                    NC_IDCODE, NC_CTAB, COM_MDLI_MOLFILE,
                    NC_SKETCH, NC_DATAFLAVOR_RXNSMILES);

    private static java.util.List<ClipboardFormat> writableReactionFormats =
            Arrays.asList(NC_SERIALIZEREACTION, NC_DATAFLAVOR_SERIALIZEDREACTION,
                    NC_CHEMDRAWINTERCHANGE,
                    NC_IDCODE, NC_CTAB, COM_MDLI_MOLFILE,
                    NC_DATAFLAVOR_RXNSMILES);


    private static java.util.List<String> nativeCliphandlerList;

    private static Class nativeCliphandlerClass = null;
    private static boolean compatibilityMode = false;
    private static int sketchwidth = 0;
    private static int sketchheight = 0;

    static {
        initDefaults();
    }

    public static void initDefaults() {
        nativeCliphandlerList = new ArrayList<>();

        if (Platform.isWindows()) {
            nativeCliphandlerList.add("com.actelion.research.gui.clipboard.JNAWinClipboardHandler");
            nativeCliphandlerList.add("com.actelion.research.gui.clipboard.NativeClipboardAccessor");
        } else if (Platform.isMacintosh()) {
            nativeCliphandlerList.add("com.actelion.research.gui.clipboard.JNAMacClipboardHandler");
        } else {
            writableMoleculeFormats = Arrays.asList(NC_SERIALIZEMOLECULE, NC_CHEMDRAWINTERCHANGE, NC_IDCODE);
        }
        sketchheight = 0;
        sketchwidth = 0;
    }

    public static void resetNativeCliphandler() {
        nativeCliphandlerClass = null;
        compatibilityMode = false;
    }

    /**
     * Get one or more Molecule(s) from the Clipboard. On all platforms the first choice is a serialized StereoMolecule.
     * Further supported formats on Windows are: MDLSK,MDLCT,MDL_MOL,CF_ENHMETAFILE with embedded sketch.
     * If no supported format is found and the clipboard contains text, which can be interpreted as molfile, then the
     * corresponding molecule is returned. If the clipboard contains one or multiple SMILES, IUPAC name(s) or idcode(s),
     * then the corresponding molecule(s) is/are returned. These can be multiple if allowMultiple is true.
     * Otherwise, if a StructureNameResolver is present, then it tries to interpret the name(s) and returns the
     * corresponding molecules, which may be limited to a certain number.
     *
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
     *
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
     *
     * @param prefer2D      if true and if the clipboard molecule has 3D coordinates, then new 2D-coords are invented
     * @param allowMultiple whether multiple molecules may be generated from clipboard text, if no serialized mol or special molecule format present
     * @param smartsMode    SmilesParser.SMARTS_MODE for parsing strings as SMILES, i.e. how SMARTS features are considered
     * @return list of molecules found or generated from SMILES, names, etc; empty list if no molecule present on the clipboard
     */
    private ArrayList<StereoMolecule> pasteMolecules(boolean prefer2D, boolean allowMultiple, int smartsMode) {
        ArrayList<StereoMolecule> molList = new ArrayList<>();

        StereoMolecule mol = isNativeClipHandler() ? pasteMoleculeNative(prefer2D) : pasteMoleculeLinux();
        if (mol != null)
            molList.add(mol);

        if (molList.isEmpty()) {
            // get StringFlavor from clipboard and try parsing it as idcode, molfile, smiles, or (if NameResolver exists) as name
            Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
            String text = null;
            try {
                text = (String) t.getTransferData(DataFlavor.stringFlavor);
            } catch (Exception ioe) {
                // ignored!
            }

            if (text != null) {
                try {
                    mol = new MolfileParser().getCompactMolecule(text);
                    if (mol != null && mol.getAllAtoms() != 0)
                        molList.add(mol);
                } catch (Exception e) {
                    // ignored!
                }

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
                                for (int i = 0; i < header.length; i++) {
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
                                } catch (Exception e) {
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
                            } else if (allowMultiple) {
                                if (unresolvedNameList == null)
                                    unresolvedNameList = new ArrayList<>();
                                unresolvedNameList.add(line);
                            }

                            line = reader.readLine();
                        }
                    } catch (IOException ioe) {
                        // ignored!
                    }
                }

                if (unresolvedNameList != null && !unresolvedNameList.isEmpty()) {
                    String[] idcodes = StructureNameResolver.resolveRemote(unresolvedNameList.toArray(new String[0]));
                    for (String idcode : idcodes) {
                        try {
                            mol = new IDCodeParser(prefer2D).getCompactMolecule(idcode);
                            if (mol != null && mol.getAllAtoms() != 0)
                                molList.add(mol);
                        } catch (Exception e) {
                            // ignored!
                        }
                    }
                }
            }
        }

        if (prefer2D) {
            for (StereoMolecule m : molList) {
                if (m.is3D()) {
                    m.ensureHelperArrays(Molecule.cHelperParities);    // to ensure stereo parities
                    new CoordinateInventor().invent(m);
                }
            }
        }
        return molList;
    }

    private StereoMolecule pasteMoleculeNative(boolean prefer2D) {
        StereoMolecule mol = null;
        Class clipHandlerClz = getNativeClipHandler();
        if (clipHandlerClz != null) {
            for (ClipboardFormat format : readableMoleculeFormats) {
                try {
                    mol = rawToMol((byte[]) clipHandlerClz.getMethod("getClipboardData", String.class)
                            .invoke(clipHandlerClz, format.value()), format, prefer2D);
                } catch (NoSuchMethodException | InvocationTargetException | IllegalAccessException e) {
                    e.printStackTrace();
                    break;
                }
                if (mol != null)
                    break;
            }
        }
        if (mol == null) { // Fallback to Java defaults
            mol = pasteMoleculeLinux();
        }
        return mol;
    }

    private static void loadNativeCliphandler() {
        Class clipHandlerClz = null;
        Iterator<String> iterator = nativeCliphandlerList.iterator();

        while (iterator.hasNext()) {
            String clz = iterator.next();
            try {
                clipHandlerClz = Class.forName(clz);
                try {
                    if (clipHandlerClz != null && (boolean) clipHandlerClz.getMethod("isInitOK").invoke(clipHandlerClz)) {
                        nativeCliphandlerClass = clipHandlerClz;
                        break;
                    }
                } catch (NoSuchMethodException | SecurityException | IllegalAccessException |
                         InvocationTargetException e2) {
                    e2.printStackTrace();
                }
            } catch (ClassNotFoundException | NoClassDefFoundError | UnsatisfiedLinkError e) {
                //e1.printStackTrace();
                System.out.printf("Attempt to load %s failed: %s\n", clz, e.getMessage());
            }
            if (nativeCliphandlerClass == null) {
                // avoid to try to use a non-available cliphandler over and over again
                System.err.println("Unable to load " + clz);
                iterator.remove();
            }
        }
    }


    private static Class getNativeClipHandler() {

        if (nativeCliphandlerClass == null) {
            loadNativeCliphandler();
        }
        return nativeCliphandlerClass;
    }

    public static boolean isNativeClipHandler() {
        Class nativeClipHandler = getNativeClipHandler();
        return nativeClipHandler != null;
    }

    private StereoMolecule pasteMoleculeLinux() {
        try {
            Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
            return (StereoMolecule) t.getTransferData(ChemistryFlavors.DF_SERIALIZED_MOLECULE);
        } catch (Exception e) {
            // ignored!
        }
        return null;
    }

    /**
     * Get a Reaction from the Clipboard
     *
     * @return Reaction or null if no reaction present
     */
    public Reaction pasteReaction() {
        Reaction rxn = isNativeClipHandler() ? pasteReactionNative() : pasteReactionLinux();

        if (rxn == null) {
            // get StringFlavor from clipboard and try parsing it as rxn-idcode, rxn-file, or reaction-smiles
            Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
            String text = null;
            try {
                text = (String) t.getTransferData(DataFlavor.stringFlavor);
            } catch (Exception ioe) {
                // ignored!
            }
            if (text != null) {
                try {
                    rxn = ReactionEncoder.decode(text, true);
                    if (rxn != null && rxn.isEmpty())
                        rxn = null;
                } catch (Exception e) {
                    // ignored!
                }
                if (rxn == null) {
                    try {
                        rxn = new RXNFileParser().getReaction(text);
                        if (rxn != null && rxn.isEmpty())
                            rxn = null;
                    } catch (Exception e) {
                        // ignored!
                    }
                }
                if (rxn == null) {
                    try {
                        rxn = new SmilesParser().parseReaction(text);
                        if (rxn != null && rxn.isEmpty())
                            rxn = null;
                    } catch (Exception e) {
                        // ignored!
                    }
                }
            }
        }

        return rxn;
    }

    public Reaction pasteReactionNative() {
        Reaction rxn = null;
        Class clipHandlerClz = getNativeClipHandler();
        if (clipHandlerClz != null) {
            for (ClipboardFormat format : readableReactionFormats) {
                try {
                    rxn = rawToRxn((byte[]) clipHandlerClz.getMethod("getClipboardData", String.class)
                            .invoke(clipHandlerClz, format.value()), format);
                } catch (NoSuchMethodException | InvocationTargetException | IllegalAccessException e) {
                    e.printStackTrace();
                    break;
                }
                if (rxn != null) break;
            }
        }
        if (rxn == null) {
            rxn = pasteReactionLinux();
        }
        return rxn;
    }

    private Reaction pasteReactionLinux() {
        try {
            Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
            return (Reaction) t.getTransferData(ChemistryFlavors.DF_SERIALIZED_REACTION);
        } catch (Exception e) {
            // ignored!
        }
        return null;
    }

    public boolean copyMolecule(String molfile) {
        StereoMolecule m = new StereoMolecule();
        MolfileParser p = new MolfileParser();
        p.parse(m, molfile);
        return copyMolecule(m);
    }

    public boolean copySizedMolecule(StereoMolecule mol, int width, int height) {
        sketchwidth = width;
        sketchheight = height;
        boolean ok = copyMolecule(mol);
        sketchwidth = 0;
        sketchheight = 0;
        return ok;
    }

    /**
     * Copies a molecule to the clipboard in various formats.
     * If available, the copy is executed via a "native" clipboard handler
     * otherwise it falls back to standard Java-based copy.
     */
    public boolean copyMolecule(StereoMolecule mol) {
        boolean ok = false;
        if (isNativeClipHandler()) {
            if (compatibilityMode)
                ok = copyMoleculeNativeClassic(mol);
            else
                ok = copyMoleculeNative(mol);
        }
        if (!ok) { // If not successful or no native handler
            MoleculeTransferable transferable = new MoleculeTransferable(mol);
            Toolkit.getDefaultToolkit().getSystemClipboard().setContents(transferable, transferable);
            ok = true;
        }
        return ok;
    }

    public static void emptyClipboard() {
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

    private boolean copyMoleculeNativeClassic(StereoMolecule mol) {
        Class clipHandlerClz;
        boolean ok = false;

        if ((clipHandlerClz = getNativeClipHandler()) != null) {
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

                ok = (boolean) method.invoke(clipHandlerClz, path,
                        molToRaw(mol, NC_CHEMDRAWINTERCHANGE), molToRaw(mol, NC_SERIALIZEMOLECULE));

                temp.delete();

            } catch (NoSuchMethodException | IOException | IllegalAccessException | InvocationTargetException e) {
                e.printStackTrace();
            }
        }
        return ok;
    }

    private boolean copyMoleculeNative(StereoMolecule mol) {
        Class clipHandlerClz;
        boolean ok = false;
        if ((clipHandlerClz = getNativeClipHandler()) != null) {
            emptyClipboard();
            for (ClipboardFormat format : writableMoleculeFormats) {
                try {
                    byte[] bytes = molToRaw(mol, format);
                    if (bytes != null) {
                        // Don't fail if any of the formats fail!
                        ok |= (boolean) clipHandlerClz.getMethod("setClipBoardData", String.class, byte[].class, boolean.class)
                                .invoke(clipHandlerClz, format.value(), bytes, false);
                    }
                } catch (InvocationTargetException | IllegalAccessException | NoSuchMethodException e) {
                    System.err.println("Can not use method: " + e.getMessage() + " trying fallback");
                    ok = copyMoleculeNativeClassic(mol);
                    if (ok)
                        compatibilityMode = true;
                    else {
                        System.err.println("Loaded cliphandler " + clipHandlerClz.getName() + " but none of the required methods to copy a Molecule was available");
                        // could be useful here (CXR: doubtful!)
                        //useNextnativeCliphandler(true);
                    }
                    break;
                }
            }
        }
        return ok;
    }

    public static StereoMolecule rawToMol(byte[] buffer, ClipboardFormat nativeFormat, boolean prefer2D) {

        StereoMolecule mol = null;
        if (buffer != null) {
            switch (nativeFormat) {
                case NC_SERIALIZEMOLECULE:
                case NC_DATAFLAVOR_SERIALIZEDMOLECULE:
                    mol = parseSerializedMolecule(buffer);
                    break;

                case NC_MOLFILE:
                case NC_MOLFILE_V3:
                    mol = parseMolFile(new String(buffer));
                    break;

                case NC_SKETCH:
                case COM_MDLI_SKETCHFILE:
                case NC_EMBEDDEDSKETCH:
                    mol = parseSketch(buffer);
                    break;

                case NC_DATAFLAVOR_IDCODE:
                    mol = parseSerializedIdCode(buffer, prefer2D);
                    break;

                case NC_IDCODE:
                    mol = parseIdCode(new String(buffer), prefer2D);
                    break;

                case NC_SMILES:
                    mol = parseSmiles(buffer);
                    break;

                case NC_CTAB:
                case COM_MDLI_MOLFILE:
                    mol = parsePascalStringsMolFile(buffer);
                    break;

                case COM_PE_CDX:
                case NC_CHEMDRAWINTERCHANGE:
                case NC_METAFILE:
                    // Unsupported
                    System.out.println("Unsupported format provided: " + nativeFormat);
                    break;

                default:
                    System.out.println("Unknown format provided: " + nativeFormat);
                    break;

            }
        }
        return mol;
    }

    public static Reaction rawToRxn(byte[] buffer, ClipboardFormat nativeFormat) {
        Reaction rxn = null;
        if (buffer != null) {
            switch (nativeFormat) {
                case NC_SERIALIZEREACTION:
                case NC_DATAFLAVOR_SERIALIZEDREACTION:
                    rxn = deSerializeReaction(buffer);
                    break;
                case NC_SKETCH:
                    rxn = parseReactionSketch(buffer);
                    break;
                case NC_IDCODE:
                    rxn = ReactionEncoder.decode(new String(buffer), true);
                    break;
                case NC_CTAB:
                case COM_MDLI_MOLFILE:
                    rxn = parsePascalStringsRxnFile(buffer);
                    break;
                case NC_DATAFLAVOR_RXNSMILES:
                    rxn = parseReactionSmiles(buffer);
                    break;
                default:
                    System.out.println("No known format provided: " + nativeFormat);
                    break;
            }
        }
        return rxn;
    }

    public static byte[] molToRaw(StereoMolecule mol, ClipboardFormat nativeFormat) {
        byte[] bytes = null;
        if (mol != null) {
            switch (nativeFormat) {
                case NC_SERIALIZEMOLECULE:
                case NC_DATAFLAVOR_SERIALIZEDMOLECULE:
                    bytes = createSerializedMolecule(mol);
                    break;
                case COM_PE_CDX:
                case NC_CHEMDRAWINTERCHANGE:
                    bytes = createCDX(mol);
                    break;
                case NC_METAFILE:
                    bytes = createMetafile(mol);
                    break;
                case NC_DATAFLAVOR_IDCODE:
                    bytes = createSerializedString(mol.getIDCode());
                    break;
                case NC_IDCODE:
                    bytes = mol.getIDCode().getBytes();
                    break;
                case NC_MOLFILE_V3:
                    bytes = createV3MolFile(mol);
                    break;
                case NC_MOLFILE:
                    bytes = createV2MolFile(mol);
                    break;
                case NC_CTAB:
                case COM_MDLI_MOLFILE:
                    bytes = createPascalStringsMolFile(mol);
                    break;
                case NC_SMILES:
                    bytes = createSmiles(mol);
                    break;
                case NC_SKETCH:
                case COM_MDLI_SKETCHFILE:
                    bytes = createSketch(mol);
                    break;
                default:
                    System.err.println("Unknown Format " + nativeFormat);
                    break;
            }
        }
        return bytes;
    }

    public static byte[] rxnToRaw(Reaction rxn, String ctab, ClipboardFormat nativeFormat) {
        byte[] bytes = null;
        if (rxn != null) {
            switch (nativeFormat) {
                case NC_SERIALIZEREACTION:
                case NC_DATAFLAVOR_SERIALIZEDREACTION:
                    bytes = serializeReaction(rxn);
                    break;
                case COM_PE_CDX:
                case NC_CHEMDRAWINTERCHANGE:
                    ChemDrawCDX cdx = new ChemDrawCDX();
                    bytes = cdx.getChemDrawBuffer(rxn);
                    break;
                case NC_CTAB:
                case COM_MDLI_MOLFILE:
                    bytes = createPascalStringsRxnFile(rxn);
                    break;
                case NC_IDCODE:
                    bytes = ReactionEncoder.encode(rxn, true, ReactionEncoder.INCLUDE_DEFAULT).getBytes();
                    break;
                case NC_DATAFLAVOR_RXNSMILES:
                    bytes = serializeReactionSmiles(rxn);
                    break;
                default:
                    System.err.println("Unkown Format " + nativeFormat);
                    break;
            }
        }
        return bytes;
    }


    public boolean copyReaction(Reaction r) {
        return copyReactionToClipboard(null, r);
    }

    private Reaction makeRXNCopy(Reaction r) {
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
            if (compatibilityMode)
                return copyReactionNativeClassic(ctab, rxn);
            else
                return copyReactionNative(ctab, rxn);
        }

        ReactionTransferable transferable = new ReactionTransferable(rxn);
        Toolkit.getDefaultToolkit().getSystemClipboard().setContents(transferable, transferable);
        return true;

    }

    public boolean copyReactionNativeClassic(String ctab, Reaction rxn) {
        Class clipHandlerClz;
        if ((clipHandlerClz = getNativeClipHandler()) != null) {
            try {
                Method method = clipHandlerClz.getMethod("copyReactionToClipboard", byte[].class, byte[].class, byte[].class);

                return (boolean) method.invoke(clipHandlerClz, rxnToRaw(rxn, ctab, NC_CTAB), rxnToRaw(rxn, ctab, NC_CHEMDRAWINTERCHANGE), rxnToRaw(rxn, ctab, NC_SERIALIZEREACTION));

            } catch (NoSuchMethodException | IllegalAccessException | InvocationTargetException e) {
                e.printStackTrace();
            }
        }
        return false;
    }


    public boolean copyReactionNative(String ctab, Reaction rxn) {
        Class clipHandlerClz = null;
        boolean ok = false;
        if ((clipHandlerClz = getNativeClipHandler()) != null) {
            emptyClipboard();
            for (ClipboardFormat format : writableReactionFormats) {
                try {
                    byte[] bytes = rxnToRaw(rxn, ctab, format);
                    if (bytes != null) {
                        ok |= (boolean) clipHandlerClz.getMethod("setClipBoardData", String.class, byte[].class, boolean.class)
                                .invoke(clipHandlerClz, format.value(), bytes, false);
                    }
                } catch (InvocationTargetException | IllegalAccessException | NoSuchMethodException e) {
                    System.err.println("Issue with: " + e.getMessage() + ", trying compatibility method");
                    ok = copyReactionNativeClassic(ctab, rxn);
                    if (ok) compatibilityMode = true;
                    else {
                        System.err.println("Loaded cliphandler " + clipHandlerClz.getName() + " but none of the required methods to copy a Reaction was available");
                    }
                    break;
                }
            }
        }
        return ok;
    }

    private static boolean writeMol2Metafile(File temp, StereoMolecule m, byte[] sketch) {
        boolean ok = false;
        try {
            ok = writeMol2Metafile(new FileOutputStream(temp), m, sketch);
        } catch (Exception e) {
            System.err.println("ClipboardHandler: Exception writing molfile " + e);
            e.printStackTrace();
        }
        return ok;
    }

    private static boolean writeMol2Metafile(OutputStream out, StereoMolecule m, byte[] sketch) throws IOException {
        int w = 300;
        int h = 200;
        WMF wmf = new WMF();
        WMFGraphics2D g = new WMFGraphics2D(wmf, w, h, Color.black, Color.white);

        Depictor2D d = new Depictor2D(m);
        d.updateCoords(g, new GenericRectangle(0, 0, w, h), AbstractDepictor.cModeInflateToMaxAVBL);
        if (sketchwidth != 0 && sketchheight != 0)
            g.scale((double) sketchwidth / d.getBoundingRect().getWidth(), (double) sketchheight / d.getBoundingRect().getHeight());
        d.paint(g);

        if (sketch != null) {
            byte[] temp = new byte[MDLSK.length + sketch.length];
            System.arraycopy(MDLSK, 0, temp, 0, MDLSK.length);
            System.arraycopy(sketch, 0, temp, MDLSK.length, sketch.length);
            wmf.escape(WMF.MFCOMMENT, temp);
        }
        wmf.writeWMF(out);
        out.close();
        return true;
    }

    public static boolean setClipBoardData(String format, byte[] buffer) {
        if (isNativeClipHandler()) {
            Class clipHandlerClz = getNativeClipHandler();
            if (clipHandlerClz != null) {
                try {
                    return (boolean) clipHandlerClz.getMethod("setClipBoardData", String.class, byte[].class)
                            .invoke(clipHandlerClz, format, buffer);
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
    public static boolean copyMetaFile(byte[] data) {
        return setClipBoardData(NC_METAFILE.value(), data);
    }

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

    public static void setNativeCliphandlerList(java.util.List<String> nativeCliphandlerList) {
        ClipboardHandler.nativeCliphandlerList = nativeCliphandlerList;
    }

    public static java.util.List<String> getNativeCliphandlerList() {
        return ClipboardHandler.nativeCliphandlerList;
    }

    public static void setCompatibilityMode(boolean compatibilityMode) {
        ClipboardHandler.compatibilityMode = compatibilityMode;
    }

    public static void setCompatibilityModeAuto() {
        Class nativeClipHandlerClz = getNativeClipHandler();
        if (nativeClipHandlerClz != null) {
            try {
                nativeClipHandlerClz.getMethod("setClipBoardData", String.class, byte[].class, boolean.class);
                ClipboardHandler.compatibilityMode = false;
            } catch (NoSuchMethodException e) {
                System.err.println("Compatibility set to true because setClipboardData(String,byte[],boolean) not available: " + e.getLocalizedMessage());
                ClipboardHandler.compatibilityMode = true;
            }
        }
    }

    public static void useNextnativeCliphandler(boolean removeFromList) {
        if (nativeCliphandlerList != null && nativeCliphandlerList.size() > 0) {
            String s = nativeCliphandlerList.remove(0);
            if (!removeFromList) nativeCliphandlerList.add(s);
            ClipboardHandler.resetNativeCliphandler();
        }
    }

    private static byte[] createSketch(StereoMolecule mol) {
        return Sketch.createSketchFromMol(mol);
    }

    private static byte[] createSmiles(StereoMolecule mol) {
        try {
            byte[] bytes = null;
            IsomericSmilesCreator creator = new IsomericSmilesCreator(mol);
            String smiles = creator.getSmiles();
            bytes = smiles.getBytes();
            return bytes;
        } catch (Throwable e) {
            System.out.printf("Unable to create SMILES %s\n",e.getMessage());
            return null;
        }
    }

    private static byte[] createPascalStringsMolFile(StereoMolecule mol) {
        try {
            MolfileCreator mf = new MolfileCreator(mol);
            ByteArrayOutputStream stream = new ByteArrayOutputStream();
            DataOutputStream os = new DataOutputStream(stream);
            String molfile = mf.getMolfile();
            String[] lines = molfile.split("\\r?\\n");
            for (String line : lines) {
                int length = line.length();
                os.write(length);
                if (length > 0)
                    os.write(line.getBytes());
            }
            os.close();
            return stream.toByteArray();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    private static byte[] createPascalStringsRxnFile(Reaction rxn) {
        try {
            RXNFileCreator mc = new RXNFileCreator(rxn);
            ByteArrayOutputStream stream = new ByteArrayOutputStream();
            DataOutputStream os = new DataOutputStream(stream);
            String rxnFile = mc.getRXNfile();
            String[] lines = rxnFile.split("\\r?\\n");
            for (String line : lines) {
                int length = line.length();
                os.write(length);
                if (length > 0)
                    os.write(line.getBytes());
            }
            os.close();
            return stream.toByteArray();
        } catch (IOException e) {
            System.out.printf("Could not create Pascal-string RXN buffer: %s\n", e.getMessage());
        }
        return null;
    }


    private static byte[] createSerializedString(String input) {
        byte[] bytes = null;
        try {
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream out = new ObjectOutputStream(bos);
            out.writeObject(input);
            out.close();
            bos.close();
            bytes = bos.toByteArray();
        } catch (IOException e) {
            e.printStackTrace(); // Report this: This should be fixed!
        }
        return bytes;
    }

    private static byte[] createSerializedMolecule(StereoMolecule mol) {
        byte[] bytes = null;
        try {
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream out = new ObjectOutputStream(bos);
            out.writeObject(mol);
            out.close();
            bos.close();
            bytes = bos.toByteArray();
        } catch (IOException e) {
            e.printStackTrace(); // Fix required, if this happens
        }
        return bytes;
    }

    private static byte[] createCDX(StereoMolecule mol) {
        try {
            byte[] bytes;
            StereoMolecule m = mol.getCompactCopy();
            for (int atom = 0; atom < m.getAllAtoms(); atom++)
                m.setAtomMapNo(atom, 0, false);
            ChemDrawCDX cdx = new ChemDrawCDX();
            bytes = cdx.getChemDrawBuffer(m);
            return bytes;
        } catch (Exception e) {
            System.out.printf("Unable to create ChemDraw Buffer: %s\n",e.getMessage());
            return null;
        }
    }

    private static byte[] createMetafile(StereoMolecule mol) {
        byte[] bytes = null;
        StereoMolecule m;
        try {
            m = mol.getCompactCopy();
            for (int atom = 0; atom < m.getAllAtoms(); atom++)
                m.setAtomMapNo(atom, 0, false);
            byte buffer[] = Sketch.createSketchFromMol(m);

            File temp = File.createTempFile("actnca", ".wmf");
            temp.deleteOnExit();

            if (writeMol2Metafile(temp, m, buffer))
                bytes = Files.readAllBytes(temp.toPath());
            temp.delete();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return bytes;
    }

    private static byte[] createV2MolFile(StereoMolecule mol) {
        byte[] bytes;
        MolfileCreator mf = new MolfileCreator(mol);
        bytes = mf.getMolfile().getBytes();
        return bytes;
    }

    private static byte[] createV3MolFile(StereoMolecule mol) {
        byte[] bytes;
        MolfileV3Creator mf = new MolfileV3Creator(mol);
        bytes = mf.getMolfile().getBytes();
        return bytes;
    }

    private static StringBuffer getMolFromMDLCT(InputStream is) throws IOException {
        StringBuffer sb = new StringBuffer();
        // In MDLCT format lines are prepended with their length.  There are
        // no newline characters (i.e. a line is just like a Pascal string)
        int len;
        while ((len = is.read()) != -1) {
            if (len > 0) {
                byte[] cbuf = new byte[len];
                //  attempt to read all characters in the "line"
                int rd = is.read(cbuf, 0, len);
                // unless we read that number of characters, something is wrong
                if (rd == -1) {
                    throw new IOException("MDLCT line specified length of " + len + " characters for line, found EOF");
                } else if (rd != len) {
                    throw new IOException("MDLCT input specified length of " + len + " characters for line, found " + rd);
                }
                sb.append(new String(cbuf));
            }
            sb.append(System.getProperty("line.separator"));
        }
        return sb;
    }

    private static StereoMolecule parsePascalStringsMolFile(byte[] buffer) {
        try {
            StringBuffer sb = getMolFromMDLCT(new ByteArrayInputStream(buffer));
            return parseMolFile(sb.toString());
        } catch (IOException e) {
            System.out.println("Unable to parse MDLCT: " + e.getMessage());
        }
        return null;
    }

    private static Reaction parsePascalStringsRxnFile(byte[] buffer) {
        try {
            StringBuffer sb = getMolFromMDLCT(new ByteArrayInputStream(buffer));
            return parseRXN(sb.toString());
        } catch (IOException e) {
            System.out.println("Unable to parse MDLCT: " + e.getMessage());
        }
        return null;
    }

    private static Reaction parseRXN(String buffer) {
        RXNFileParser parser = new RXNFileParser();
        Reaction reaction = new Reaction();
        try {
            parser.parse(reaction, buffer);
            return reaction;
        } catch (Exception e) {
            System.out.printf("Unable to parse RXN file: %s\n", e.getMessage());
        }
        return null;
    }

    private static StereoMolecule parseSmiles(byte[] buffer) {
        try {
            StereoMolecule mol = new StereoMolecule();
            SmilesParser parser = new SmilesParser();
            parser.parse(mol, buffer);
            return mol;
        } catch (Exception ex) {
            System.out.println("Unable to parse SMILES: " + ex.getMessage());
        }
        return null;
    }

    private static StereoMolecule parseMolFile(String buffer) {
        MolfileParser p = new MolfileParser();
        StereoMolecule mol = new StereoMolecule();
        if (!p.parse(mol, buffer)) {
            mol = null;
            System.out.println("Warn: Parsing CTAB failed\n");
        }
        return mol; // Can be null
    }

    private static StereoMolecule parseIdCode(String idCode, boolean prefer2D) {
        StereoMolecule mol;
        try {
            mol = new StereoMolecule();
            IDCodeParser parser = new IDCodeParser(prefer2D);
            parser.parse(mol, idCode);
            if (mol.getAllAtoms() == 0) {
                mol = null;
            }
        } catch (Exception ex) {
            System.out.printf("IdCode '%s' could not be parsed: %s\n", idCode, ex.getMessage());
            ex.printStackTrace();
            mol = null;
        }
        return mol;
    }

    private static StereoMolecule parseSerializedMolecule(byte[] buffer) {
        StereoMolecule mol = null;
        try (ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(buffer))) {
            Object o = is.readObject();
            is.close();
            if (o instanceof StereoMolecule) {
                mol = (StereoMolecule) o;
            }
        } catch (Exception e) {
            System.out.printf("Serialized molecule could not be read: %s\n",e.getMessage());
            e.printStackTrace();
        }
        return mol;
    }

    private static StereoMolecule parseSketch(byte[] buffer) {
        StereoMolecule mol;
        try {
            mol = new StereoMolecule();
            if (!Sketch.createMolFromSketchBuffer(mol, buffer)) {
                mol = null;
            }
        } catch (IOException ex) {
            mol = null;
            System.out.printf("Warn: parsing Sketch failed: %s\n" ,ex);
            ex.printStackTrace();
        }
        return mol;
    }

    private static StereoMolecule parseSerializedIdCode(byte[] buffer, boolean prefer2D) {
        StereoMolecule mol = new StereoMolecule();
        try {
            ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(buffer));
            Object o = is.readObject();
            if (o instanceof String) {
                mol = parseIdCode((String) o, prefer2D);
            }
            is.close();
            return mol;
        } catch (Exception ex) {
            System.out.printf("Warn: parsing serialized IDCode failed: %s\n" ,ex);
            ex.printStackTrace();
        }
        return null;
    }

    private static byte[] serializeReaction(Reaction rxn) {
        byte[] bytes = null;
        try {
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream out = new ObjectOutputStream(bos);
            out.writeObject(rxn);
            out.close();
            bos.close();
            bytes = bos.toByteArray();
        } catch (IOException e) {
            System.out.printf("Warn: Serializing Reaction failed: %s\n" ,e);
            e.printStackTrace();
        }
        return bytes;
    }

    private static Reaction parseReactionSketch(byte[] buffer) {
        Reaction rxn;
        try {
            rxn = new Reaction();
            if (!Sketch.createReactionFromSketchBuffer(rxn, buffer)) {
                rxn = null;
            }
        } catch (IOException ex) {
            System.out.printf("Warn: parsing Reaction Sketch failed: %s\n" ,ex);
            rxn = null;
        }
        return rxn;
    }

    private static Reaction deSerializeReaction(byte[] buffer) {
        Reaction rxn = null;
        try {
            ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(buffer));
            Object o = is.readObject();
            if (o instanceof Reaction) {
                rxn = (Reaction) o;
            }
            is.close();
        } catch (Exception ex) {
            System.out.printf("Warn: Deserializing Reaction failed: %s\n" ,ex);
            ex.printStackTrace();
        }
        return rxn;
    }

    private static Reaction parseReactionSmiles(byte[] buffer) {
        Reaction rxn = null;
        try {
            ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(buffer));
            Object o = is.readObject();
            if (o instanceof String) {
                rxn = new SmilesParser().parseReaction((String) o);
            }
            is.close();
        } catch (Exception ex) {
            System.out.printf("Warn: Parsing Reaction SMILES failed: %s\n" ,ex);
            ex.printStackTrace();
        }
        return rxn;

    }

    private static byte[] serializeReactionSmiles(Reaction rxn) {
        String reactionSmiles = IsomericSmilesCreator.createReactionSmiles(rxn);
        byte[] bytes = null;
        try {
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream out = new ObjectOutputStream(bos);
            out.writeObject(reactionSmiles);
            out.close();
            bos.close();
            bytes = bos.toByteArray();
        } catch (IOException e) {
            System.out.printf("Warn: Serializing Reaction SMILES failed: %s\n" ,e);
            e.printStackTrace();
        }
        return bytes;
    }


}
