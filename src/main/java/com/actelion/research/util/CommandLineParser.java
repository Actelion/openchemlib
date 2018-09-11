package com.actelion.research.util;

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

    public void add(String command, String value) {

        if(hmCommandValue.containsKey(command)){
            throw new RuntimeException("Contains already command " + command);
        }

        hmCommandValue.put(command,value);
    }

    public String get(String command) {
        return hmCommandValue.get(command);
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
                throw new RuntimeException("Wrong command line argument");
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
