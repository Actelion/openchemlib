package com.actelion.research.chem.dnd;

import javafx.scene.input.DataFormat;

public class ChemistryDataFormats {
	final static DataFormat get(String mimetype) {
		DataFormat df = DataFormat.lookupMimeType(mimetype);
		if (df == null)
			df = new DataFormat(mimetype);
		return df;
	}

	public static final DataFormat DF_SERIALIZED_OBJECT = get("application/x-java-serialized-object");
	public static final DataFormat DF_MDLMOLFILE = get("chemical/x-mdl-molfile");
	public static final DataFormat DF_MDLMOLFILEV3 = get("chemical/x-mdl-molfilev3");
	public static final DataFormat DF_SMILES = get("chemical/x-daylight-smiles");
	public static final DataFormat DF_IDCODE = get("chemical/x-openmolecules-idcode");
}
