package com.actelion.research.chem.io.pdb.parser;

import java.util.HashMap;
import java.util.List;

/**
 * HetNameParser
 * <p>Copyright: Idorsia Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Idorsia Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * Created by korffmo1 on 10.04.18.
 */
public class HetNameParser {

    private int indexLine;

    private HashMap<String, String> hmId_Name;

    void parse(List<String> liRaw, int indexLine) {

        hmId_Name = new HashMap<>();

        String tag = PDBFileParser.TAG_HETNAM;

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

                hmId_Name.put(nameId, name);

                nameId = nameIdLine;

                sb = new StringBuilder();
            }

            if(sb.length()>0){
                sb.append(" ");
            }

            String nameLine = l.substring(15).trim();
            sb.append(nameLine);

            indexLine++;
        }

        if(sb.length()>0){
            String name = sb.toString();

            hmId_Name.put(nameId, name);
        }

        this.indexLine = indexLine;
    }

    public int getIndexLine() {
        return indexLine;
    }

    public HashMap<String, String> getHMId_Name() {
        return hmId_Name;
    }
}