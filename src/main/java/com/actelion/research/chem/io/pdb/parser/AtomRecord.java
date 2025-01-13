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
 * @author Modest v. Korff
 */

package com.actelion.research.chem.io.pdb.parser;

import com.actelion.research.chem.PeriodicTable;

/**
 * ModelAtom
 * Created by korffmo1 on 11.04.18.
 */
public class AtomRecord implements Comparable<AtomRecord> {

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

	public String getAltLoc() {
		return altLoc;
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

	@Override
	public int compareTo(AtomRecord o) {
		return Integer.compare(serialId, o.serialId);
	}
}
