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

package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.descriptor.DescriptorEncoder;
import com.actelion.research.chem.descriptor.DescriptorHandlerFlexophore;
import com.actelion.research.chem.descriptor.flexophore.generator.ConstantsFlexophoreGenerator;
import com.actelion.research.util.datamodel.IntVec;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

public class DistHistEncoder {

	private static DistHistEncoder INSTANCE;
	
	
	private int iNBitsEntriesCountOneHistogram;
	
	private int iNBitsPos;
	
	private int iNBitsConsequentEntries;
	
	private int iNBitsCountOneField;
	
	/**
	 * 
	 */
	public DistHistEncoder() {
		init();
	}
	
	private void init(){
		
		int iNumConf = DescriptorHandlerFlexophore.NUM_CONFORMATIONS;
		
		int iLenHist = ConstantsFlexophoreGenerator.BINS_HISTOGRAM;

		//
		// So many bits are needed to describe the number of fields > 0 in one histogram.
		//
		int maxNumEntriesOneHist = ConstantsFlexophoreGenerator.BINS_HISTOGRAM / 2;
		iNBitsEntriesCountOneHistogram = 1;
		
		while(Math.pow(2,iNBitsEntriesCountOneHistogram) < maxNumEntriesOneHist){
			iNBitsEntriesCountOneHistogram++;
		}
		
		iNBitsConsequentEntries = iNBitsEntriesCountOneHistogram;
		
		//
		// So many bits are needed to describe the position of one entry in a histogram.
		//
		iNBitsPos = 1;
		while(Math.pow(2,iNBitsPos)<iLenHist){
			iNBitsPos++;
		}

		//
		// So many bits are needed to describe the counts for one field.
		//
		
		iNBitsCountOneField = 1;
		while(Math.pow(2,iNBitsCountOneField) < iNumConf){
			iNBitsCountOneField++;
		}
	}
	
	/**
	 * 
	 * @param dh
	 * @return number_entries_hist pos_in_hist num_consequent_entries count1 count2  pos_in_hist num_consequent_entries count1 count2 count3 ...
	 */
	public String encodeHistograms(DistHist dh){
		
		int numPPNodes = dh.getNumPPNodes();
		
		if(dh.getNumPPNodes()==1){
			return "";
		}
				
		int nApproxBitsNeeded = 0;
		
		for (int i = 0; i < numPPNodes; i++) {
			for (int j = i+1; j < numPPNodes; j++) {
				
				 
				byte [] arrHist = dh.getDistHist(i,j);
				
				int iOccupied = 0;
				for (int k = 0; k < arrHist.length; k++) {
					if(arrHist[k]>0) 
						iOccupied++;
				}
				
				nApproxBitsNeeded += iNBitsEntriesCountOneHistogram + iOccupied * (iNBitsPos + iNBitsCountOneField + iNBitsConsequentEntries);
			}
		}

		boolean [] arr = new boolean [nApproxBitsNeeded];
		
		int posArray=0;
		
		
		for (int i = 0; i < numPPNodes; i++) {
			for (int j = i+1; j < numPPNodes; j++) {
				
				 
				byte [] arrHist = dh.getDistHist(i,j);
				
				int nFieldsOccupied = 0;
				for (int k = 0; k < arrHist.length; k++) {
					if(arrHist[k]>0) 
						nFieldsOccupied++;
				}
				
				// Set number of entries for histogram.
				for (int k = 0; k < iNBitsEntriesCountOneHistogram; k++) {
					if((nFieldsOccupied&1)==1){
						arr[posArray]=true;
					}
					nFieldsOccupied >>>= 1;
					posArray++;
				}
				
				boolean histogramProcessed = false;
				
				int k = 0;
				while(!histogramProcessed){
					
					if(arrHist[k] > 0) {
						int posInHist = k;
						
						// Encode position
						for (int l = 0; l < iNBitsPos; l++) {
							if((posInHist&1)==1){
								arr[posArray]=true;
							}
							posInHist >>>= 1;
							posArray++;
						}
						
						int nConsequentFieldsOcc=1;
						for (int l = k+1; l < arrHist.length; l++) {
							if(arrHist[l] > 0) {
								nConsequentFieldsOcc++;
							} else {
								break;
							}
						}
						
						int nConsequentFieldsOcc2Bit = nConsequentFieldsOcc;
						for (int m = 0; m < iNBitsConsequentEntries; m++) {
							if((nConsequentFieldsOcc2Bit & 1)==1){
								arr[posArray]=true;
							}
							nConsequentFieldsOcc2Bit >>>= 1;
							posArray++;
						}
						
						int startSet = k;
						int endSet = k+nConsequentFieldsOcc;
						
						for (int l = startSet; l < endSet; l++) {
							
							int counts = arrHist[l];
							
							for (int m = 0; m < iNBitsCountOneField; m++) {
								if((counts & 1)==1){
									arr[posArray]=true;
								}
								counts >>>= 1;
								posArray++;
							}
						}
						
						k += nConsequentFieldsOcc;
					} 
					
					// If the last field in the histogram is occupied k becomes here 40.
					// The next field would be 41.
					
					k++;
					if(k >= arrHist.length){
						histogramProcessed=true;
					}
					
				}
				
			}
		}
		
		boolean [] arrTruncated = new boolean [posArray];
		
		System.arraycopy(arr, 0, arrTruncated, 0, arrTruncated.length);
				
		IntVec iv = new IntVec(arrTruncated);
		
		String s = new String(new DescriptorEncoder().encode(iv.get()), StandardCharsets.UTF_8);
				
		return s;
	}
	
	public void decodeHistograms(String s, MolDistHist mdh){
				
		IntVec iv = new IntVec(new DescriptorEncoder().decode(s));
		
		List<byte[]> liHist = new ArrayList<byte[]>();
		
		boolean decodeFinished=false;
		
		int pos = 0;
		while(!decodeFinished){
			
			int nEntriesInHistogram = getDecodedValue(iv, pos, iNBitsEntriesCountOneHistogram);
			
			if(nEntriesInHistogram==0){
				break;
			}
			
			pos += iNBitsEntriesCountOneHistogram;
			
			byte[] arrHist = new byte [ConstantsFlexophoreGenerator.BINS_HISTOGRAM];
			
			int processedEntries = 0;
			
			while(processedEntries != nEntriesInHistogram){
						
				int positionField = getDecodedValue(iv, pos, iNBitsPos);
				pos += iNBitsPos;
				
				int consequentCounts = getDecodedValue(iv, pos, iNBitsConsequentEntries);
				pos += iNBitsConsequentEntries;
				
				for (int j = 0; j < consequentCounts; j++) {
					
					int counts = getDecodedValue(iv, pos, iNBitsCountOneField);
					pos += iNBitsCountOneField;
					
					arrHist[positionField++]=(byte)counts;
				}
				
				processedEntries += consequentCounts;
				
				
			}
			
			liHist.add(arrHist);
			
			if(pos+iNBitsEntriesCountOneHistogram >= iv.sizeBits()){
				break;
			}
		}
		
		int size = getNumNodes(liHist.size());
		
		if(size==0){
			throw new RuntimeException("Number of pharmacophore points is 0.");
		}
				
		int cc=0;
		for (int i = 0; i < size; i++) {
			for (int j = i+1; j < size; j++) {
				mdh.setDistHist(i,j, liHist.get(cc++));
			}
		}
		
	}

	private static int getNumNodes(int nHistogramms){
		
		int nNodes = 0;
		
		int cc=0;
		
		while(cc <= nHistogramms){
			nNodes++;
			cc += nNodes;
		}
		
		return nNodes;
		
	}
	
	private static int getDecodedValue(IntVec iv, int posStart, int widthInBits){
		
		int posEnd = posStart + widthInBits;
		
		int value = 0;
		
		// The first bit is the lowest.
		for (int i = posEnd-1; i >= posStart; i--) {
			value <<=1;
			
			if(iv.isBitSet(i)){
				value |= 1;
			}
		}
		
		return value;
	}
	
	
	public static DistHistEncoder getInstance(){
		
		if(INSTANCE == null){
			INSTANCE = new DistHistEncoder();
		}
		
		return INSTANCE;
	}

}
