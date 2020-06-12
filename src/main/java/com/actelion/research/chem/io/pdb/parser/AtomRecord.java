package com.actelion.research.chem.io.pdb.parser;

import com.actelion.research.chem.PeriodicTable;

/**
 * ModelAtom
 * <p>Copyright: Idorsia Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Idorsia Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * Created by korffmo1 on 11.04.18.
 */
public class AtomRecord {

    private int serialId;

    // From field atom name
    private int atomicNo;
    // From field atom name
    private String atomName;

    // Alternate location indicator.
    private String altLoc;

    private String residueName;

    private String chainId;
    
    private String insertionCode;

    private int resNum;

    private double x,y,z;

    private double occupancy;

    private double tempFactor;

    private String element;

    private String charge;

    private String anisou;
    
    private boolean isTerminalC;


    public AtomRecord(int serialId,
                     String atomName,
                     String altLoc,
                     String residueName,
                     String chainId,
                     int resNum,
                     String insertionCode,
                     double x,
                     double y,
                     double z,
                     double occupancy,
                     double tempFactor,
                     String element) {

        this.serialId = serialId;
        this.atomName = atomName;
        this.altLoc = altLoc;
        this.residueName = residueName;
        this.chainId = chainId;
        this.resNum = resNum;
        this.insertionCode = insertionCode;
        this.x = x;
        this.y = y;
        this.z = z;
        this.occupancy = occupancy;
        this.tempFactor = tempFactor;
        this.element = element;
        this.atomicNo = PeriodicTable.number(element);
        isTerminalC = false;
    }
    

    public int getSerialId() {
        return serialId;
    }

    public String getAnisou() {
        return anisou;
    }
    
    public String getAtomName() {
    	return atomName;
    }
    
    public boolean isTerminalC() {
    	return isTerminalC;
    }
    
    public double getX() {
    	return x;
    }
    
    public double getY() {
    	return y;
    }
    
    public double getZ() {
    	return z;
    }
    
    public String getInsertionCode() {
    	return insertionCode;
    }
    
    public int getAtomicNo() {
    	return atomicNo;
    }
    
    public String getResName() {
    	return residueName;
    }
    
    public String getChainID() {
    	return chainId;
    }
    
    public String getString() {
    	StringBuilder sb = new StringBuilder(residueName);
    	sb.append(" ");
    	sb.append(resNum);
    	sb.append(insertionCode);
    	sb.append(" ");
    	sb.append(chainId);
    	return sb.toString();
    }
    
    public int getResNum() {
    	return resNum;
    }
    
    
    public void setAnisou(String anisou) {
        this.anisou = anisou;
    }
    
    public void setTerminalC(boolean isTerminalC) {
    	this.isTerminalC = isTerminalC;
    }

}
