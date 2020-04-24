package com.actelion.research.calc.combinatorics;

import com.actelion.research.util.ArrayUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * 
 * 
 * CombinationGenerator
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Oct 12, 2012 MvK: Start implementation
 */
public class CombinationGenerator {
	
	/**
	 * Get all possible index combinations, order independent. 
	 * For a list containing the numbers from  0 to 
	 * nObjects for arrays of the size 'sampleSize'.
	 * @param nObjects so many indices will be permuted.
	 * @param sampleSize Size of the array containing the permutations.
	 * @return
	 */
	public static List<int[]> getAllOutOf(int nObjects, int sampleSize) {

		
		List<int[]> li = new ArrayList<int[]>();
		
		int [] arrCounters = new int [sampleSize];
		
		if(nObjects==sampleSize){
			int [] arr = new int [sampleSize];
			for (int i = 0; i < arr.length; i++) {
				arr[i]=i;
			}
			li.add(arr);
			return li;
		}else if (sampleSize==1){
			for (int i = 0; i < nObjects; i++) {
				int [] arr = new int [1];
				arr[0]=i;
				li.add(arr);
			}
			return li;
		}else if (sampleSize>nObjects){
			return null;
		}
		
		
		// init
		for (int i = 0; i < arrCounters.length; i++) {
			arrCounters[i]=i;
		}
		boolean bProceed = true;
		
		while(bProceed){
			
			int [] arr = new int [sampleSize];
			for (int i = 0; i < sampleSize; i++) {
				arr[i]=arrCounters[i];
			}
			li.add(arr);
			
			
			int depth = sampleSize-1;
			arrCounters[depth]++;
			
			boolean bCounterFieldReset = false;
			
			if(arrCounters[depth]>=nObjects){
				bCounterFieldReset=true;
			}
			
			while(bCounterFieldReset){
				bCounterFieldReset=false;
				depth--;
				arrCounters[depth]++;
				for (int i = depth + 1; i < sampleSize; i++) {
					arrCounters[i] = arrCounters[i-1]+1;
					if(arrCounters[i] >= nObjects){
						bCounterFieldReset=true;
					}
				}
				
				if(depth==0)
					break;
			}
			if(bCounterFieldReset)
				bProceed = false;
		}
		
		return li; 
	}

	public static void main(String[] args) {

		int sizeList = 5;

		for (int i = 2; i < sizeList+1; i++) {

			List<int[]> li = CombinationGenerator.getAllOutOf(sizeList, i);

			for (int j = 0; j < li.size(); j++) {

				int [] a = li.get(j);

				System.out.println(ArrayUtils.toString(a));

			}
		}

	}

}
