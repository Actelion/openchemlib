package com.actelion.research.io;

import java.io.BufferedReader;
import java.io.IOException;

/**
 * Modest v. Korff
 * Idorsia Pharmaceuticals Ltd.
 * 05.05.2022 Start implementation
 **/
public class BufferedReaderHelper {

    public static void skipUntilLineMatchesRegEx(BufferedReader br, String regex) throws NoSuchFieldException, IOException {
        int limit = 10000;

        skipUntilLineMatchesRegEx(br, regex, limit);
    }

    public static String skipUntilLineMatchesRegEx(BufferedReader br, String regex, int limit) throws NoSuchFieldException, IOException {

        String line = br.readLine();

        boolean match = false;

        if(line.matches(regex)){
            match = true;
        }

        int cc=0;

        while(!match){
            line = br.readLine();
            if(line==null) break;
            if(line.matches(regex)){
                match = true;
            }
            cc++;
            if(cc > limit){
                break;
            }
        }

        String lineMatch = null;

        if(match){
            lineMatch = line;
        } else {
            throw new NoSuchFieldException("Regex " + regex + " was not found.");
        }

        return lineMatch;


    }

}
