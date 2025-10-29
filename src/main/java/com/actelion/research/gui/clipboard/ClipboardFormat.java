/*
 * Copyright (c) 2024-2025
 * Christian Rufener
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
