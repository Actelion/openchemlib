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

	private boolean isHetAtom;

	private int serialId;

    // From field atom name
    private int atomicNo;
    // From field atom name
    private String labelAtomName;

    // Alternate location indicator.
    private String labelAltID;

	private int labelSeqID;

	private String authCompID;

    private String authAsymId;
    
    private String insertionCode;

    private int authSeqID;

    private double x,y,z;

    private double occupancy;

    private double tempFactor;

    private String element;

    private String charge;

    private String anisou;
    
    private boolean isTerminalC;

	public AtomRecord() {}

	public void setAtomAndCompName(String labelAtomName, String authCompID) {
		this.labelAtomName = labelAtomName;
		this.authCompID = authCompID;
		this.authAsymId = "A";
	}

	public AtomRecord(boolean isHetAtom,
			         int serialId,
                     String labelAtomName,
                     String labelAltID,
					 int labelSeqID,
                     String authCompID,
                     String authAsymId,
                     int authSeqID,
                     String insertionCode,
                     double x,
                     double y,
                     double z,
                     double occupancy,
                     double tempFactor,
                     String element) {

		this.isHetAtom = isHetAtom;
        this.serialId = serialId;
        this.labelAtomName = labelAtomName;
        this.labelAltID = labelAltID;
		this.labelSeqID = labelSeqID;
        this.authCompID = authCompID;
        this.authAsymId = authAsymId;
        this.authSeqID = authSeqID;
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

	public boolean isHetAtom() {
		return isHetAtom;
	}

    public int getSerialId() {
        return serialId;
    }

    public String getAnisou() {
        return anisou;
    }

	public String getLabelAltID() {
		return labelAltID;
	}

	public int getLabelSeqID() {
		return labelSeqID;
	}

	public String getLabelAtomName() {
    	return labelAtomName;
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

	public double getOccupancy() {
		return occupancy;
	}

    public String getInsertionCode() {
    	return insertionCode;
    }
    
    public int getAtomicNo() {
    	return atomicNo;
    }
    
    public String getResName() {
    	return authCompID;
    }
    
    public String getChainID() {
    	return authAsymId;
    }
    
    public String getString() {
    	StringBuilder sb = new StringBuilder(authCompID);
    	sb.append(" ");
    	sb.append(authSeqID);
    	sb.append(insertionCode);
    	sb.append(" ");
    	sb.append(authAsymId);
    	return sb.toString();
    }

    public int getAuthSeqID() {
    	return authSeqID;
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
