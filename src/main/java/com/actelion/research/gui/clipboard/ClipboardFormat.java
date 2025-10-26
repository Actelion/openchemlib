package com.actelion.research.gui.clipboard;

/**
 * Copyright (c) 2009-2025
 * Christian Rufener async.ch
 * All rights reserved
 * User: christian
 * Date: 26.10.25
 * Time: 09:21
 */
public enum ClipboardFormat {
    NC_SERIALIZEMOLECULE("ACT_MOLECULE"),                       // Java serialized Molecule
    NC_SERIALIZEREACTION("ACT_REACTION"),                       // Java serialized Reaction
    NC_CTAB("MDLCT"),                                           // MDL Molfile as Pascal-String encoded Binary block
    NC_MOLFILE("MDL_MOL"),                                      // MDL Molfile as Text
    NC_SKETCH("MDLSK"),                                         // MDL Sketch
    NC_EMBEDDEDSKETCH("MDLSK_EMBEDDED"),                        // Metafile containing MDL Sketch
    NC_IDCODE("IDCODE"),                                        // IdCode as Text
    NC_CHEMDRAWINTERCHANGE("ChemDraw Interchange Format"),      // ChemDraw Clipboard Format
    NC_METAFILE("CF_METAFILEPICT"),                             // Windows Metafile
   // NC_DIB("CF_DIB"),                                         // Window Device Independent Bitmapt
    NC_MOLFILE_V3("MOLFILE_V3"),                                // MDL Version 3 Molfile as Text
    NC_SMILES("SMILES"),                                        // SMILES as Text
    COM_MDLI_SKETCHFILE("com.mdli.sketchfile"),                 // MDL Sketch as used by ChemDraw (Mac only?)
    COM_PE_CDX("com.perkinelmer.chemdraw.cdx-clipboard"),       // ChemDraw Clipboard Format (Mac?)
    COM_MDLI_MOLFILE("com.mdli.molfile"),                       // MDL Molfile as Pascal-String encoded Binary block
    NC_DATAFLAVOR_SERIALIZEDMOLECULE("JAVA_DATAFLAVOR:application/x-java-serialized-object; class=com.actelion.research.chem.StereoMolecule"),   // Java Serialized Stereo Molecule by Java DataTransfer
    NC_DATAFLAVOR_IDCODE("JAVA_DATAFLAVOR:chemical/x-openmolecules-idcode; class=java.lang.String"),                                             // Java Serialized IdCode String by Java DataTransfer
    NC_DATAFLAVOR_SERIALIZEDREACTION("JAVA_DATAFLAVOR:application/x-java-serialized-object; class=com.actelion.research.chem.reaction.Reaction"),// Java Serialized Reaction by Java DataTransfer
    NC_DATAFLAVOR_RXNSMILES("JAVA_DATAFLAVOR:chemical/x-daylight-reactionsmiles; class=java.lang.String")                                        // Java Serialized Reaction SMILES String by Java DataTransfer
    ;

    private final String value;

    ClipboardFormat(String format) {
        this.value = format;
    }

    public String value() {
        return value;
    }

    @Override
    public String toString() {
        return value;
    }
}
