package com.actelion.research.chem.io.pdb.converter;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.IO;

import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.List;

/**
 * AminoAcidsLabeledContainer
 *
 * The txt file was generated with com.actelion.research.chem.parsers.pdb.ProcessLabeledAminoAcids
 *
 * <p>Copyright: Idorsia Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Idorsia Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * Created by korffmo1 on 13.04.18.
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
