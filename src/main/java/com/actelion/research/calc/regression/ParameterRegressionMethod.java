package com.actelion.research.calc.regression;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

/**
 * ParameterRegressionMethod
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 06.12.18.
 *
 * The Comparable method may order according the parameter.
 *
 */


public abstract class ParameterRegressionMethod implements Comparable<ParameterRegressionMethod>, Serializable {

    // Method name, kNN, PLS etc.
    public static final String TAG_NAME = "Name";

    // For filtering for identical parameter in DataWarrior.
    public static final String TAG_STRING = "ParameterToString";

    protected Properties properties;

    public ParameterRegressionMethod() {
        properties = new Properties();
    }

    public ParameterRegressionMethod(String name) {
        properties = new Properties();
        properties.setProperty(TAG_NAME, name);
    }

    public ParameterRegressionMethod(ParameterRegressionMethod p) {
        properties = new Properties();
        copy(p);
    }

    public void copy(ParameterRegressionMethod orig){
        setName(orig.getName());
    }


    public void setName(String name) {
        properties.setProperty(TAG_NAME, name);
    }

    public String getName() {
        return properties.getProperty(TAG_NAME);
    }


    public static List<String> getHeader(){

        List<String> li = new ArrayList<>();

        li.add(TAG_NAME);
        li.add(TAG_STRING);

        return li;
    }

    public Properties getProperties() {
        return properties;
    }

    /**
     * Decodes the properties to the needed variables.
     */
    protected abstract void decodeProperties2Parameter();

    public void write(File fiProperties) throws IOException {
        properties.store(new FileWriter(fiProperties), "");
    }

    public String write2String() throws IOException {

        StringWriter stringWriter = new StringWriter();

        properties.store(stringWriter, "");

        return stringWriter.toString();
    }

    public void read(File fiProperties) throws IOException {
        properties.load(new FileReader(fiProperties));
        decodeProperties2Parameter();
    }

    public void read(String s) throws IOException {
        StringReader stringReader = new StringReader(s);
        properties.load(stringReader);
        decodeProperties2Parameter();
    }

    @Override
    public boolean equals(Object obj) {

        ParameterRegressionMethod p = (ParameterRegressionMethod)obj;

        if(!getName().equals(p.getName())){
            return false;
        }

        return true;
    }
}
