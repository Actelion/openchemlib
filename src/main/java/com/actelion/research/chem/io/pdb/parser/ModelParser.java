package com.actelion.research.chem.io.pdb.parser;

import java.util.ArrayList;
import java.util.List;

/**
 * ModelParser
 * <p>Copyright: Idorsia Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Idorsia Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * Created by korffmo1 on 11.04.18.
 */
public class ModelParser {

    private int indexLine;

    public void parse(List<String> liRaw, int indexLine, List<AtomRecord> protAtomRecords, 
    		List<AtomRecord> hetAtomRecords) {
    	
        String tagAtom = PDBFileParser.TAG_ATOM;
        String tagHeteroAtom = PDBFileParser.TAG_HETATM;

        String l0 = liRaw.get(indexLine);
        if(!l0.startsWith(tagAtom)
                && !l0.startsWith(tagHeteroAtom)
                && !l0.startsWith(PDBFileParser.TAG_MODEL)) {
            throw new RuntimeException("Error in parsing atoms.");
        }

        this.indexLine = indexLine;

        int start = indexLine;
        
        // rather use Optional here! 
        AtomRecord lastRecord = null;
        for (int i = start; i < liRaw.size(); i++) {

            String line = liRaw.get(i);
            if(line.startsWith(tagAtom)) {
                AtomRecord modelAtom = parseAtom(line);
                lastRecord = modelAtom;
                protAtomRecords.add(modelAtom);
            }
             else if(line.startsWith(tagHeteroAtom)) {

                AtomRecord modelAtom = parseAtom(line);
                lastRecord = modelAtom;
                hetAtomRecords.add(modelAtom);

            } else if (line.startsWith(PDBFileParser.TAG_ANISOU)) {

                continue;

            } else if (line.startsWith(PDBFileParser.TAG_TER)) { // End of chain or hetero atom group

                if(lastRecord!=null)
                	lastRecord.setTerminalC(true);

            } else if (line.startsWith(PDBFileParser.TAG_MODEL)){

                continue;

            } else if (line.startsWith(PDBFileParser.TAG_ENDMDL)){
                if(lastRecord!=null)
                	lastRecord.setTerminalC(true);

            } else {
                this.indexLine = i;
                break;
            }
        }


    }

    private AtomRecord parseAtom(String line){
        int serialId = Integer.parseInt(line.substring(6,11).trim());

        String atomName = line.substring(12,16).trim();

        String altLoc = line.substring(16,17).trim();

        String residueName = line.substring(17,20).trim();

        String chainId = line.substring(21,22).trim();

        int resSeq = Integer.parseInt(line.substring(22,26).trim());
        
        String insertionCode = line.substring(26,27).trim();

        double x = Double.parseDouble(line.substring(30,38).trim());
        double y = Double.parseDouble(line.substring(38,46).trim());
        double z = Double.parseDouble(line.substring(46,54).trim());

        double occupancy = Double.parseDouble(line.substring(54,60).trim());

        double tempFactor = Double.parseDouble(line.substring(60,66).trim());

        String element = line.substring(76,78).trim();
        element = element.toLowerCase();
        element = element.substring(0, 1).toUpperCase() + element.substring(1);

        AtomRecord modelAtom = new AtomRecord(serialId,
                atomName,
                altLoc,
                residueName,
                chainId,
                resSeq,
                insertionCode,
                x,
                y,
                z,
                occupancy,
                tempFactor,
                element);

        return modelAtom;
    }

    public int getIndexLine() {
        return indexLine;
    }


}
