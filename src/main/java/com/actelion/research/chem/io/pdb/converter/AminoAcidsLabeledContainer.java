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
 * @author Modest v. Korff
 */

package com.actelion.research.chem.io.pdb.converter;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.IO;

import java.io.InputStream;
import java.util.HashMap;
import java.util.List;

/**
 * AminoAcidsLabeledContainer
 *
 * The txt file was generated with com.actelion.research.chem.parsers.pdb.ProcessLabeledAminoAcids
 */
public enum AminoAcidsLabeledContainer { // enum singleton pattern;
	
	INSTANCE;

    public static final int AMINO_ACIDS = 20;

    private HashMap<String, AminoAcidLabeled> hmAbbreviation_AminoAcidLabeled;

    private AminoAcidsLabeledContainer() {
    	read();
    }

    private void read()  {
    	try {
        InputStream inputStream= this.getClass().getResourceAsStream(ConstantsAminoAcidsLabeled.RESOURCE_AAS_LABELED);

        List<String> li = IO.readLines2List(inputStream);

        inputStream.close();

        // Remove comments.
        for (int i = li.size() - 1; i >= 0; i--) {
            if(li.get(i).startsWith("#")){
                li.remove(i);
            }
        }

        IDCodeParser parser = new IDCodeParser();

        hmAbbreviation_AminoAcidLabeled = new HashMap<>();
        for (String line : li) {

            String [] arr = line.split("\t");

            StereoMolecule mol = parser.getCompactMolecule(arr[0], arr[1]);

            String key = arr[3].toUpperCase();

            AminoAcidLabeled aminoAcidLabeled = new AminoAcidLabeled(mol, arr[2], key);

            hmAbbreviation_AminoAcidLabeled.put(key, aminoAcidLabeled);

        }
    	}
        catch(Exception e) {
        	System.err.println("no Amino-Acid File detected!");
        }
    }

    public AminoAcidLabeled get(String threeLetterCodeUpperCase) {
        return hmAbbreviation_AminoAcidLabeled.get(threeLetterCodeUpperCase);
    }
}
