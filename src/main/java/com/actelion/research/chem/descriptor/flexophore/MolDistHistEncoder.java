/*
 * Copyright (c) 1997 - 2022
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
 * 3. Neither the name of the copyright holder nor the
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

package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.descriptor.DescriptorEncoder;
import com.actelion.research.chem.descriptor.DescriptorHandlerFlexophore;
import com.actelion.research.util.datamodel.IntVec;

import java.nio.charset.StandardCharsets;

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

		String atoms = mdh.getNodeAtoms() == null ? "" : " " + new String(new DescriptorEncoder().encodeIntArray2D(mdh.getNodeAtoms()), StandardCharsets.UTF_8);
		
		return strNodes + strDistHist + atoms;
	}
	
	public MolDistHist decode(byte [] arr){
		return decode(new String(arr, StandardCharsets.UTF_8));
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
			mdh.setNodeAtoms(new DescriptorEncoder().decodeIntArray2D(st[2].getBytes(StandardCharsets.UTF_8)));

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
