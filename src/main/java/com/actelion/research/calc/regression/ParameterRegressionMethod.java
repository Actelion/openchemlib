/*
 * Copyright (c) 1997 - 2018
 * Idorsia Pharmaceuticals Ltd.
 * Hegenheimermattweg 91
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
