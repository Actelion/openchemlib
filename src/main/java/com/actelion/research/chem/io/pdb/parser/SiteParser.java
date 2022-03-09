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

import java.text.ParseException;
import java.util.HashMap;
import java.util.List;

/**
 * SiteParser
 *
 * Site records specify residues comprising catalytic, co-factor, anti-codon, regulatory or other
 * essential sites or environments surrounding ligands present in the structure.
 *
 * Created by korffmo1 on 11.04.18.
 */
public class SiteParser {



    private int indexLine;

    private HashMap<String, String> hmId_Site;

    void parse(List<String> liRaw, int indexLine) throws ParseException {

        hmId_Site = new HashMap<>();

        String tag = PDBFileParser.TAG_SITE;

        String l0 = liRaw.get(indexLine);
        if(!l0.startsWith(tag)) {
            throw new RuntimeException("Error in parsing " + tag);
        }

        StringBuilder sb = new StringBuilder();

        int start = indexLine;

        String nameId = l0.substring(11, 14);

        for (int i = start; i < liRaw.size(); i++) {

            String l = liRaw.get(i);

            if(!l.startsWith(tag)) {
                break;
            }

            String nameIdLine = l.substring(11, 14);

            if(!nameId.equals(nameIdLine)) {

                String name = sb.toString();

                hmId_Site.put(nameId, name);

                nameId = nameIdLine;

                sb = new StringBuilder();
            }

            if(sb.length()>0){
                sb.append(" ");
            }

            String nameLine = l.substring(18).trim();
            sb.append(nameLine);

            indexLine++;
        }

        if(sb.length()>0){
            String name = sb.toString();

            hmId_Site.put(nameId, name);
        }

        this.indexLine = indexLine;
    }

    public int getIndexLine() {
        return indexLine;
    }

    public HashMap<String, String> getHMId_Site() {
        return hmId_Site;
    }
}