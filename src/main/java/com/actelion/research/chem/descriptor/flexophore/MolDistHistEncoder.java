/*
 * Copyright (c) 2020.
 * Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 *  This file is part of DataWarrior.
 *
 *  DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.
 *
 *  DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License along with DataWarrior.
 *  If not, see http://www.gnu.org/licenses/.
 *
 *  @author Modest v. Korff
 *
 */

package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.descriptor.DescriptorEncoder;
import com.actelion.research.chem.descriptor.DescriptorHandlerFlexophore;
import com.actelion.research.util.datamodel.IntVec;

public class MolDistHistEncoder {
	
	private static MolDistHistEncoder INSTANCE;

	private DistHistEncoder distHistEncoder;
	
	public MolDistHistEncoder() {
		
		distHistEncoder = new DistHistEncoder();
		
	}
	
	public String encode(MolDistHist mdh){
		
		if(!mdh.isFinalized())
			mdh.realize();
				
		String strNodes = encodeByteVec(mdh.getArrNode());
				
		String strDistHist = " " + distHistEncoder.encodeHistograms(mdh);

		String atoms = mdh.getNodeAtoms() == null ? "" : " " + new String(new DescriptorEncoder().encodeIntArray2D(mdh.getNodeAtoms()));
		
		return strNodes + strDistHist + atoms;
	}
	
	public MolDistHist decode(byte [] arr){
		return decode(new String(arr));
	}
	
	public MolDistHist decode(String s){

		if(s.equals(DescriptorHandlerFlexophore.FAILED_STRING)){
			return DescriptorHandlerFlexophore.FAILED_OBJECT;
		}
		
		String[] st = s.split(" ");
		
		byte [] arrNodes = decodeNodes(st[0]);

		int nNodes=0;
		int pos=0;
		while(arrNodes[pos] > 0){
			pos += arrNodes[pos] * PPNode.getNumBytesEntry() + 1;
			
			nNodes++;
			
			if(pos >= arrNodes.length){
				break;
			}
		}

		byte [] arrNodesTrunc = new byte [pos];
		
		System.arraycopy(arrNodes, 0, arrNodesTrunc, 0, arrNodesTrunc.length);
				
		MolDistHist mdh = new MolDistHist(nNodes);
				
		mdh.setArrNode(arrNodesTrunc);
		
		if(st.length >= 2 && st[1].length() != 0)
			distHistEncoder.decodeHistograms(st[1], mdh);

		if (st.length >= 3)
			mdh.setNodeAtoms(new DescriptorEncoder().decodeIntArray2D(st[2].getBytes()));

		return mdh;

	}
	
	/**
	 * 
	 * @param s
	 * @return
	 */
	public static byte [] decodeNodes(String s){
		
		IntVec iv = new IntVec(new DescriptorEncoder().decode(s));
		
		byte [] arr = new byte [iv.sizeBytes()];
		
		int sum=0;
		for (int i = 0; i < arr.length; i++) {
			arr[i] = (byte)(iv.getByte(i) & 0xFF);
			sum += Math.abs(arr[i]);	
		}
		
		if(sum==0){
			throw new RuntimeException("Node vector contains only 0's!");
		}
		
		return arr;
	}

	public static MolDistHistEncoder getInstance(){
		
		if(INSTANCE == null){
			INSTANCE = new MolDistHistEncoder();
		}
		
		return INSTANCE;
	}

	private static String encodeByteVec(byte [] arr){
		
		double ratio = Integer.SIZE / Byte.SIZE;
		int sizeIntVec = (int)((arr.length / ratio) + ((ratio-1) / ratio));
		IntVec iv = new IntVec(sizeIntVec);
		for (int i = 0; i < arr.length; i++) {
			iv.setByte(i, (arr[i] & 0xFF));
		}
		String s= new String (new DescriptorEncoder().encode(iv.get()));
		
		return s;
	}
	


}
