package com.actelion.research.chem;

import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.ObjectInputStream;




/**
 * <p>Title: Mercury</p>
 * <p>Description: Actelion Electronic Lab Notebook</p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: </p>
 * @author Christian Rufener
 * @version 1.0
 */

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
