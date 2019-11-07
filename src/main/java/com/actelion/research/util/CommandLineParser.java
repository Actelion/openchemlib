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
*/

package com.actelion.research.util;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 * CommandLineParser
 *
 * Simple command line parser. Creates a key value table.
 *
 * Created by korffmo1 on 13.07.17.
 */
public class CommandLineParser {

    public static final String HELP = "-h";

    private HashMap<String,String> hmCommandValue;

    public CommandLineParser() {
        hmCommandValue=new HashMap<String,String>();
    }

    public CommandLineParser(String args[]) {
        hmCommandValue=new HashMap<String,String>();
        parse(args);
    }

    public void add(String command, String value) {

        if(hmCommandValue.containsKey(command)){
            throw new RuntimeException("Contains already command " + command);
        }

        hmCommandValue.put(command,value);
    }

    public String get(String command) {
        return hmCommandValue.get(command);
    }

    public File getAsFile(String command) {
        return new File(get(command));
    }

    public int getAsInt(String command) {
        return Integer.parseInt(get(command));
    }

    public boolean contains(String command) {
        return hmCommandValue.containsKey(command);
    }

    public boolean help() {
        return hmCommandValue.containsKey(HELP);
    }

    public boolean checkCommandWithValue(String command) {

        if(!contains(command)){
            throw new RuntimeException("Argument '" + command + "' missing.");
        }

        if(get(command)==null){
            throw new RuntimeException("Value for '" + command + "' missing.");
        }

        return true;
    }


    public int parse(String [] args){

        int index=0;

        while (index<args.length){

            String s0 = args[index];

            if(!s0.startsWith("-")){
                throw new RuntimeException("Wrong command line argument '" + s0 + "'");
            }

            String s1 = null;
            if(index<args.length-1){
                if(!args[index+1].startsWith("-")){
                    s1 = args[index+1];
                    index++;
                }
            }

            index++;

            hmCommandValue.put(s0, s1);

        }

        return hmCommandValue.size();
    }

    public int getNumArguments() {
        return hmCommandValue.size();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();

        List<String> li = new ArrayList<String>(hmCommandValue.keySet());

        Collections.sort(li);

        for (String k : li) {
            sb.append(k);
            sb.append("\t");
            sb.append(hmCommandValue.get(k));
            sb.append("\n");
        }

        return sb.toString();
    }

    public static void main(String[] args) {

        String c = "-f file -D -C command";

        String [] argsCmd = c.split(" ");

        CommandLineParser clp = new CommandLineParser();

        clp.parse(argsCmd);

        System.out.println(clp.toString());

    }
}
