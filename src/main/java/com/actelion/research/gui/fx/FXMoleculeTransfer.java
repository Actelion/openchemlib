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
 * Additionally, molecules are provided as serialized StereoMolecule, molfile V2, molfile V3, SMILES and idcode
 * using the MIME types of com.actelion.research.chem.dnd.ChemistryFlavors. The serialized molecule uses the
 * same identifier as AWT, such that Swing and JavaFX can exchange molecule objects directly.
 * When reading, the string formats are also recognized,
 * if they were put onto the clipboard by AWT, which encodes them as Java serialized Strings
 * using 'JAVA_DATAFLAVOR:' prefixed identifiers.
 */
public class FXMoleculeTransfer {
	private static final String MIME_IDCODE = "chemical/x-openmolecules-idcode";
	private static final String MIME_MOLFILE = "chemical/x-mdl-molfile";
	private static final String MIME_MOLFILE_V3 = "chemical/x-mdl-molfilev3";
	private static final String MIME_SMILES = "chemical/x-daylight-smiles";
	private static final String AWT_PREFIX = "JAVA_DATAFLAVOR:";

	// AWT's native name of ChemistryFlavors.DF_SERIALIZED_MOLECULE. JavaFX and AWT both use Java serialization
	// for custom formats. Thus, using the same identifier lets JavaFX and Swing exchange StereoMolecule objects.
	private static final String SERIALIZED_MOLECULE_ID = AWT_PREFIX+"application/x-java-serialized-object; class="+StereoMolecule.class.getName();

	public static final DataFormat DF_SERIALIZED_MOLECULE = getDataFormat(SERIALIZED_MOLECULE_ID);
	public static final DataFormat DF_MOLFILE = getDataFormat(MIME_MOLFILE);
	public static final DataFormat DF_MOLFILE_V3 = getDataFormat(MIME_MOLFILE_V3);
	public static final DataFormat DF_SMILES = getDataFormat(MIME_SMILES);
	public static final DataFormat DF_IDCODE = getDataFormat(MIME_IDCODE);

	// same order as in ChemistryFlavors.MOLECULE_FLAVORS
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
		content.put(DF_SERIALIZED_MOLECULE, new StereoMolecule(mol));	// copy, because local drags may deliver the object itself
		try {
			content.put(DF_MOLFILE, new MolfileCreator(mol).getMolfile());
		}
		catch (Exception e) {}
		try {
			content.put(DF_MOLFILE_V3, new MolfileV3Creator(mol).getMolfile());
		}
		catch (Exception e) {}
		try {
			content.put(DF_SMILES, new IsomericSmilesCreator(mol).getSmiles());
		}
		catch (Exception e) {}
		Canonizer canonizer = new Canonizer(mol);
		String idcode = canonizer.getIDCode()+" "+canonizer.getEncodedCoordinates(true);
		content.put(DF_IDCODE, idcode);
		content.putString(idcode);
		return content;
	}

	/**
	 * @param clipboard
	 * @return whether the clipboard contains any content, which may be interpreted as molecule
	 */
	public static boolean hasMolecule(Clipboard clipboard) {
		if (clipboard.hasString())
			return true;
		for (DataFormat format : clipboard.getContentTypes()) {
			if (isSerializedMolecule(format))
				return true;
			for (String mimeType : MOLECULE_MIME_TYPES)
				if (matches(format, mimeType))
					return true;
			}
		return false;
	}

	/**
	 * Creates a molecule from the clipboard content. Chemistry specific formats are preferred in the same
	 * order as by com.actelion.research.gui.dnd.MoleculeDropAdapter and ClipboardHandler: serialized molecule,
	 * idcode, molfile V3, molfile V2, SMILES. If a format cannot be interpreted, the next one is tried.
	 * If includeText is true, plain text is tried to be interpreted as molfile, idcode, SMILES,
	 * or (if a StructureNameResolver exists) as name, in that order.
	 * @param clipboard
	 * @param includeText whether to fall back to interpreting plain text
	 * @return molecule or null
	 */
	public static StereoMolecule readMolecule(Clipboard clipboard, boolean includeText) {
		StereoMolecule mol = getSerializedMolecule(clipboard);
		if (mol == null)
			mol = parseIDCode(getString(clipboard, MIME_IDCODE));
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

	/**
	 * JavaFX deserializes custom formats itself. Within the same JVM it may even deliver the original object.
	 * @return copy of a serialized molecule put there by JavaFX or AWT; or null
	 */
	private static StereoMolecule getSerializedMolecule(Clipboard clipboard) {
		for (DataFormat format : clipboard.getContentTypes()) {
			if (isSerializedMolecule(format)) {
				try {
					Object o = clipboard.getContent(format);
					if (o instanceof Molecule) {
						StereoMolecule mol = new StereoMolecule((Molecule)o);
						if (mol.getAllAtoms() != 0)
							return mol;
					}
				}
				catch (Exception e) {}
			}
		}
		return null;
	}

	/**
	 * @return true, if one of format's identifiers is the one used by AWT for serialized StereoMolecules
	 */
	private static boolean isSerializedMolecule(DataFormat format) {
		for (String id : format.getIdentifiers())
			if (id.replace(" ", "").equalsIgnoreCase(SERIALIZED_MOLECULE_ID.replace(" ", "")))
				return true;
		return false;
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
