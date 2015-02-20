package com.actelion.research.chem;

public class StructureTransferData  implements java.io.Serializable
{
    static final long serialVersionUID = 0x2005CAFE;
//    String name_ , id_ , comment_;
    StructureInfo si_ = null;
    ExtendedMolecule mol_ = new ExtendedMolecule();
    public StructureTransferData(ExtendedMolecule mol, StructureInfo si)
    {
        mol_ = mol;
        si_ = si;
    }

    public StructureTransferData()
    {
    }
    
    public ExtendedMolecule getMolecule()
    {
        return mol_;
    }
    
    public StructureInfo getStructureInfo()
    {
        return si_;
    }
/*
    public String getName()
    {
        return name_;
    }
    public String getId()
    {
        return id_;
    }
    
    public String getComment()
    {
        return comment_;
    }
 */
}
