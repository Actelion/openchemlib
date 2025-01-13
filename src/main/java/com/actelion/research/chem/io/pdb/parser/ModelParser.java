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

package com.actelion.research.chem.io.pdb.parser;

import java.util.List;
import java.util.TreeSet;

/**
 * ModelParser
 * Created by korffmo1 on 11.04.18.
 */
public class ModelParser {

    private int indexLine;

    public void parse(List<String> liRaw, int indexLine, TreeSet<AtomRecord> protAtomRecords,
                      TreeSet<AtomRecord> hetAtomRecords) {
    	
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

        //invert y and z coordinates for compatibility with Java coordinate system (analogously to Molfileparser & Mol2FileParser)
        double x = Double.parseDouble(line.substring(30,38).trim());
        double y = -Double.parseDouble(line.substring(38,46).trim());
        double z = -Double.parseDouble(line.substring(46,54).trim());

        double occupancy = (line.length() < 60) ? 1.0 : Double.parseDouble(line.substring(54,60).trim());

        double tempFactor = (line.length() < 66) ? 50.0 : Double.parseDouble(line.substring(60,66).trim());

        String element = (line.length() < 78) ? atomName.substring(0, 1) : line.substring(76,78).trim();
        element = element.toLowerCase();
        element = element.substring(0, 1).toUpperCase() + element.substring(1);

        return new AtomRecord(serialId,
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
    }

    public int getIndexLine() {
        return indexLine;
    }


}
