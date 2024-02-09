package com.actelion.research.chem;

import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.ObjectInputStream;




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

@Deprecated
public class SDFileMolecule extends StereoMolecule implements java.io.Serializable
{
    static final long serialVersionUID = 0x2003CAFE;
    public transient static final String NAME_FIELD = "NAME";
    public transient static final String ID_FIELD = "ID";
    public transient static final String CMNT_FIELD = "CMNT";
    public transient static final String DENSITY_FIELD = "DENSITY";
    public transient static final String PURITY_FIELD = "PURITY";
    
    
    transient java.util.TreeMap map_ = new java.util.TreeMap();
    transient java.util.ArrayList fields_ = new java.util.ArrayList();

    public SDFileMolecule()
    {
    }
    public SDFileMolecule(Molecule m)
    {
        super(m);
    }
    public SDFileMolecule(SDFileMolecule m)
    {
        super(m);
        map_ = new java.util.TreeMap(m.map_);
        fields_ = new java.util.ArrayList(m.fields_);
    }

    public void setFieldData(String field, String data)
    {
        map_.put(field,data);
        if (!fields_.contains(field))
            fields_.add(field);
    }

    public String getFieldData(String field)
    {
        Object o = map_.get(field);
        if (o instanceof String)
            return (String)o;
        return null;
    }

    public boolean moveFieldToStart(String field)
    {
        if (fields_.contains(field)) {
            fields_.remove(field);
            fields_.add(0,field);
            return true;
        }
        return false;
    }
    public String[] getFields()
    {
//        java.util.Set set = map_.keySet();
        java.util.List set = fields_;    
        String res[] = new String[set.size()];
        set.toArray(res);
        return res;
    }

    public String toString()
    {
        return "SDFileMolecule " + super.toString();
    }
        private void writeObject(ObjectOutputStream stream) throws IOException 
        {
            System.out.println("SDFIleMoleucle writeObject");
//            stream.writeObject(map_);
//            stream.writeObject(fields_);
        }

        private void readObject(ObjectInputStream stream) throws IOException 
        {
            System.out.println("SDFIleMoleucle readObject");
/*            try {
                fields_ = (java.util.ArrayList)stream.readObject();
            } catch (ClassNotFoundException e) {
                System.err.println("Error reading Object SDFileMolecule : " + e);
            }
            try {
                map_ = (java.util.TreeMap)stream.readObject();
            } catch (ClassNotFoundException e) {
                System.err.println("Error reading Object SDFileMolecule : " + e);
            }
 */
        }

}
