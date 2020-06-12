package com.actelion.research.chem.io.pdb.parser;


import java.text.ParseException;
import java.util.AbstractMap;
import java.util.HashMap;
import java.util.List;
import java.util.AbstractMap.SimpleEntry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * RemarkParser
 * <p>Copyright: Idorsia Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Idorsia Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * Created by korffmo1 on 09.04.18.
 */
public class RemarkParser {

    private int indexLine;

    private HashMap<Integer, String> hmNo_Remark;

    public RemarkParser() {

    }

    void parse(List<String> liRaw, int indexLine) throws ParseException {

        this.indexLine = indexLine;

        hmNo_Remark = new HashMap<>();

        int ccRemark = 0;

        // Optional Mandatory for a re-refined structure
        if(RemarkParser.getMatcherRemark(ccRemark, liRaw.get(indexLine)).find()) {
        	SimpleEntry<String, Integer> siIndex = parseRemark(liRaw, indexLine, ccRemark);

            hmNo_Remark.put(ccRemark, siIndex.getKey());

            indexLine = siIndex.getValue();
        }

        ccRemark++;

        if(RemarkParser.getMatcherRemark(ccRemark, liRaw.get(indexLine)) .find()) {
        	SimpleEntry<String, Integer> siIndex = parseRemark(liRaw, indexLine, ccRemark);
            hmNo_Remark.put(ccRemark, siIndex.getKey());
            indexLine = siIndex.getValue();
        }

        // Remark 2 mandatory
        ccRemark++;
        if(RemarkParser.getMatcherRemark(ccRemark, liRaw.get(indexLine)) .find()) {
        	SimpleEntry<String, Integer> siIndex = parseRemark(liRaw, indexLine, ccRemark);
            hmNo_Remark.put(ccRemark, siIndex.getKey());
            indexLine = siIndex.getValue();
        } //else {
       //     throw new RuntimeException("Missing " + PDBFileParser.TAG_REMARK2);
        //}

        // Remark 3 mandatory
        ccRemark++;
        if(RemarkParser.getMatcherRemark(ccRemark, liRaw.get(indexLine)) .find()) {
        	SimpleEntry<String, Integer> siIndex = parseRemark(liRaw, indexLine, ccRemark);
            hmNo_Remark.put(ccRemark, siIndex.getKey());
            indexLine = siIndex.getValue();
        } //else {
        //    throw new RuntimeException("Missing " + PDBFileParser.TAG_REMARK3);
        //}

        Pattern p = getPatternRemark();

        boolean searchRemarks = true;

        while (searchRemarks) {

            String line = liRaw.get(indexLine);

            Matcher matcher = p.matcher(line);

            if(matcher.find()) {

                String [] arr = line.split("[ ]+");

                ccRemark = Integer.parseInt(arr[1]);

                SimpleEntry<String, Integer> siIndex = parseRemark(liRaw, indexLine, ccRemark);

                hmNo_Remark.put(ccRemark, siIndex.getKey());

                indexLine = siIndex.getValue();

            } else {
                searchRemarks = false;
            }
        }

        this.indexLine = indexLine;
    }

    public int getIndexLine() {
        return indexLine;
    }

    public HashMap<Integer, String> getHmNo_Remark() {
        return hmNo_Remark;
    }

    private static SimpleEntry<String, Integer> parseRemark(List<String> liRaw, int indexLine, int indexRemark) throws ParseException {


        Pattern pattern = RemarkParser.getPatternRemark(indexRemark);

        // String titleSub0 = l0.substring(tag.length()).trim();

        StringBuilder sb = new StringBuilder();

        int start = indexLine;

        for (int i = start; i < liRaw.size(); i++) {

            String l = liRaw.get(i);

            Matcher matcher = pattern.matcher(l);

            if(matcher.find()) {

                String [] arr = l.split("[ ]+");
                sb.append(" ");
                for (int j = 2; j < arr.length; j++) {
                    sb.append(arr[j]);

                    if(j < arr.length-1) {
                        sb.append(" ");
                    }
                }

                indexLine++;
            } else {
                break;
            }

        }

        AbstractMap.SimpleEntry<String, Integer> siTextIndex = new AbstractMap.SimpleEntry<>(sb.toString(), indexLine);
 

        return siTextIndex;

    }

    static Pattern getPatternRemark() {

        String p = PDBFileParser.TAG_REMARK + "[ ]+";

        return Pattern.compile(p);
    }

    static Pattern getPatternRemark(int index) {

        String p = PDBFileParser.TAG_REMARK + "[ ]+" + index;

        return Pattern.compile(p);
    }

    static Matcher getMatcherRemark(int indexRemark, String s){

        Pattern p = getPatternRemark(indexRemark);

        Matcher matcher = p.matcher(s);

        return matcher;
    }
}
