package com.actelion.research.gui.fx;

import com.actelion.research.chem.*;
import com.actelion.research.chem.name.StructureNameResolver;
import javafx.scene.input.Clipboard;
import javafx.scene.input.ClipboardContent;
import javafx.scene.input.DataFormat;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InvalidClassException;
import java.io.ObjectInputStream;
import java.io.ObjectStreamClass;
import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;

/**
 * Helper to put molecules onto and to read molecules from a JavaFX Clipboard or Dragboard
 * (which extends Clipboard) such that they can be exchanged with Swing components and other applications.
 * The common denominator between JavaFX and Swing/AWT is plain text containing idcode and coordinates,
 * which is also what com.actelion.research.gui.dnd.MoleculeTransferable provides as stringFlavor.
 * Additionally, molecules are provided as idcode, molfile V2, molfile V3 and SMILES using the MIME types
 * of com.actelion.research.chem.dnd.ChemistryFlavors. When reading, these formats are also recognized,
 * if they were put onto the clipboard by AWT, which encodes them as Java serialized Strings
 * using 'JAVA_DATAFLAVOR:' prefixed identifiers.
 */
public class FXMoleculeTransfer {
	private static final String MIME_IDCODE = "chemical/x-openmolecules-idcode";
	private static final String MIME_MOLFILE = "chemical/x-mdl-molfile";
	private static final String MIME_MOLFILE_V3 = "chemical/x-mdl-molfilev3";
	private static final String MIME_SMILES = "chemical/x-daylight-smiles";
	private static final String AWT_PREFIX = "JAVA_DATAFLAVOR:";

	public static final DataFormat DF_IDCODE = getDataFormat(MIME_IDCODE);
	public static final DataFormat DF_MOLFILE = getDataFormat(MIME_MOLFILE);
	public static final DataFormat DF_MOLFILE_V3 = getDataFormat(MIME_MOLFILE_V3);
	public static final DataFormat DF_SMILES = getDataFormat(MIME_SMILES);

	private static final String[] MOLECULE_MIME_TYPES = { MIME_IDCODE, MIME_MOLFILE_V3, MIME_MOLFILE, MIME_SMILES };

	private FXMoleculeTransfer() {}

	/**
	 * DataFormat identifiers must be unique within a JVM. Therefore, we re-use existing ones.
	 */
	private static DataFormat getDataFormat(String mimeType) {
		DataFormat format = DataFormat.lookupMimeType(mimeType);
		return (format != null) ? format : new DataFormat(mimeType);
	}

	/**
	 * Creates clipboard content containing the molecule in all supported formats.
	 * @param mol
	 * @return content for Clipboard.setContent() or Dragboard.setContent()
	 */
	public static ClipboardContent createContent(StereoMolecule mol) {
		ClipboardContent content = new ClipboardContent();
		Canonizer canonizer = new Canonizer(mol);
		String idcode = canonizer.getIDCode()+" "+canonizer.getEncodedCoordinates(true);
		content.putString(idcode);
		content.put(DF_IDCODE, idcode);
		try {
			content.put(DF_MOLFILE_V3, new MolfileV3Creator(mol).getMolfile());
		}
		catch (Exception e) {}
		try {
			content.put(DF_MOLFILE, new MolfileCreator(mol).getMolfile());
		}
		catch (Exception e) {}
		try {
			content.put(DF_SMILES, new IsomericSmilesCreator(mol).getSmiles());
		}
		catch (Exception e) {}
		return content;
	}

	/**
	 * @param clipboard
	 * @return whether the clipboard contains any content, which may be interpreted as molecule
	 */
	public static boolean hasMolecule(Clipboard clipboard) {
		if (clipboard.hasString())
			return true;
		for (DataFormat format : clipboard.getContentTypes())
			for (String mimeType : MOLECULE_MIME_TYPES)
				if (matches(format, mimeType))
					return true;
		return false;
	}

	/**
	 * Creates a molecule from the clipboard content. Chemistry specific formats are preferred.
	 * If includeText is true, plain text is tried to be interpreted as molfile, idcode, SMILES,
	 * or (if a StructureNameResolver exists) as name, in that order.
	 * @param clipboard
	 * @param includeText whether to fall back to interpreting plain text
	 * @return molecule or null
	 */
	public static StereoMolecule readMolecule(Clipboard clipboard, boolean includeText) {
		StereoMolecule mol = parseIDCode(getString(clipboard, MIME_IDCODE));
		if (mol == null)
			mol = parseMolfile(getString(clipboard, MIME_MOLFILE_V3));
		if (mol == null)
			mol = parseMolfile(getString(clipboard, MIME_MOLFILE));
		if (mol == null)
			mol = parseSmiles(getString(clipboard, MIME_SMILES));
		if (mol == null && includeText && clipboard.hasString())
			mol = parseText(clipboard.getString());
		return mol;
	}

	private static StereoMolecule parseText(String text) {
		if (text == null || text.trim().isEmpty())
			return null;

		StereoMolecule mol = parseMolfile(text);	// not trimmed, because the first molfile lines may be empty
		text = text.trim();
		if (mol == null)
			mol = parseIDCode(text);
		if (mol == null)
			mol = parseSmiles(text);
		if (mol == null) {
			mol = StructureNameResolver.resolve(text);
			if (mol != null && mol.getAllAtoms() == 0)
				mol = null;
			}
		return mol;
	}

	private static StereoMolecule parseIDCode(String idcode) {
		if (idcode == null)
			return null;
		try {
			StereoMolecule mol = new IDCodeParser(true).getCompactMolecule(idcode.trim());
			return (mol != null && mol.getAllAtoms() != 0) ? mol : null;
		}
		catch (Exception e) {
			return null;
		}
	}

	private static StereoMolecule parseMolfile(String molfile) {
		if (molfile == null)
			return null;
		try {
			StereoMolecule mol = new MolfileParser().getCompactMolecule(molfile);
			return (mol != null && mol.getAllAtoms() != 0) ? mol : null;
		}
		catch (Exception e) {
			return null;
		}
	}

	private static StereoMolecule parseSmiles(String smiles) {
		if (smiles == null)
			return null;
		try {
			StereoMolecule mol = new StereoMolecule();
			new SmilesParser().parse(mol, smiles.trim());
			return (mol.getAllAtoms() != 0) ? mol : null;
		}
		catch (Exception e) {
			return null;
		}
	}

	/**
	 * @return content of the first data format matching the mime type, either FX style or AWT style; or null
	 */
	private static String getString(Clipboard clipboard, String mimeType) {
		for (DataFormat format : clipboard.getContentTypes()) {
			if (matches(format, mimeType)) {
				try {
					String s = decodeString(clipboard.getContent(format));
					if (s != null && !s.isEmpty())
						return s;
				}
				catch (Exception e) {}
			}
		}
		return null;
	}

	/**
	 * @return true, if one of format's identifiers equals mimeType or is an AWT encoded flavor of mimeType
	 */
	private static boolean matches(DataFormat format, String mimeType) {
		for (String id : format.getIdentifiers()) {
			if (id.startsWith(AWT_PREFIX))
				id = id.substring(AWT_PREFIX.length());
			int index = id.indexOf(';');
			if (index != -1)
				id = id.substring(0, index);
			if (id.trim().equalsIgnoreCase(mimeType))
				return true;
		}
		return false;
	}

	/**
	 * Depending on platform and source, custom formats are delivered as String or as raw bytes.
	 * AWT puts Java serialized String objects onto the clipboard for its custom String based flavors.
	 */
	private static String decodeString(Object content) {
		if (content instanceof String)
			return (String)content;

		byte[] bytes = null;
		if (content instanceof ByteBuffer) {
			ByteBuffer buffer = ((ByteBuffer)content).duplicate();
			bytes = new byte[buffer.remaining()];
			buffer.get(bytes);
		}
		else if (content instanceof byte[]) {
			bytes = (byte[])content;
		}

		if (bytes == null || bytes.length == 0)
			return null;

		if (bytes.length > 2 && (bytes[0] & 0xFF) == 0xAC && (bytes[1] & 0xFF) == 0xED)
			return deserializeString(bytes);

		int length = bytes.length;
		while (length > 0 && bytes[length-1] == 0)	// remove C-string terminators
			length--;
		return new String(bytes, 0, length, StandardCharsets.UTF_8);
	}

	/**
	 * Deserializes a String object. For security reasons no other classes are allowed to be instantiated.
	 */
	private static String deserializeString(byte[] bytes) {
		try (ObjectInputStream is = new StringOnlyObjectInputStream(new ByteArrayInputStream(bytes))) {
			Object o = is.readObject();
			return (o instanceof String) ? (String)o : null;
		}
		catch (Exception e) {
			return null;
		}
	}

	private static class StringOnlyObjectInputStream extends ObjectInputStream {
		public StringOnlyObjectInputStream(InputStream in) throws IOException {
			super(in);
		}

		@Override
		protected Class<?> resolveClass(ObjectStreamClass desc) throws IOException {
			throw new InvalidClassException(desc.getName(), "Only String objects are accepted");
		}

		@Override
		protected Class<?> resolveProxyClass(String[] interfaces) throws IOException {
			throw new InvalidClassException("Proxy classes are not accepted");
		}
	}
}
