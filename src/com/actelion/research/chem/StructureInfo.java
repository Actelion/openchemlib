package com.actelion.research.chem;

public class StructureInfo implements java.io.Serializable
{
    static final long serialVersionUID = 0x2005CAFE;
    public static final String FLD_NAME = "FLD_NAME";
    public static final String FLD_ID = "FLD_ID";
    public static final String FLD_CMNT = "FLD_CMNT";
    public static final String FLD_DENSITY = "FLD_DENSITY";
    public static final String FLD_PURITY = "FLD_PURITY";
    
    java.util.Map data_ = new java.util.HashMap();
//    private String name;
//    private String id;
//    private String comment;

    public StructureInfo() {}
    
    public StructureInfo(String name,String id, String comment,float density, String purity)
    {
        data_.put(FLD_NAME,name);
        data_.put(FLD_ID,id);
        data_.put(FLD_CMNT,comment);
		if (!java.lang.Float.isNaN(density))
			data_.put(FLD_DENSITY,Float.toString(density));
		if (purity != null)
			data_.put(FLD_PURITY,purity);
		
    }

/*
    public StructureInfo(String name,String id, String comment)
    {
        this.name = name;
        this.id = id;
        this.comment = comment;
    }
*/
    public String getName()
    {
        Object o = data_.get(FLD_NAME);
        return o != null ? o.toString() : null;
    }

    public String getId()
    {
        Object o = data_.get(FLD_ID);
        return o != null ? o.toString() : null;
    }

    public String getComment()
    {
        Object o = data_.get(FLD_CMNT);
        return o != null ? o.toString() : null;
    }

	public String getDensity()
	{
        Object o = data_.get(FLD_DENSITY);
        return o != null ? o.toString() : null;
	}
	
	public String getPurity()
	{
        Object o = data_.get(FLD_PURITY);
        return o != null ? o.toString() : null;
	}

    public void setData(String field,Object data)
    {
        data_.put(field,data);
    }

    public Object getData(String field)
    {
        return data_.get(field);
    }
}
    
