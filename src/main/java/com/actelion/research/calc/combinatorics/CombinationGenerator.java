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
		boolean proceed = true;
		
		while(proceed){
			
			int [] arr = new int [sampleSize];
			for (int i = 0; i < sampleSize; i++) {
				arr[i]=arrCounters[i];
			}
			li.add(arr);
			
			
			int depth = sampleSize-1;
			arrCounters[depth]++;
			boolean counterFieldReset = false;
			
			if(arrCounters[depth]>=nObjects){
				counterFieldReset=true;
			}
			
			while(counterFieldReset){
				counterFieldReset=false;
				depth--;
				arrCounters[depth]++;
				for (int i = depth + 1; i < sampleSize; i++) {
					arrCounters[i] = arrCounters[i-1]+1;
					if(arrCounters[i] >= nObjects){
						counterFieldReset=true;
					}
				}
				
				if(depth==0)
					break;
			}
			if(counterFieldReset)
				proceed = false;
		}
		
		return li; 
	}

	public static List<int[]> getCombinations(List<int[]> li){

		int nCombinations = 1;
		for (int[] arr : li) {
			nCombinations *= arr.length;
		}

		List<int[]> liComb = new ArrayList<>(nCombinations);

		for (int i = 0; i < nCombinations; i++) {
			int [] arrComb = new int[li.size()];
			liComb.add(arrComb);
		}

		int nCombCol=1;
		for (int col = 0; col < li.size(); col++) {

			nCombCol *= li.get(col).length;

			int nRepetitions = nCombinations / nCombCol;

			int [] arr = li.get(col);

			int indexArr = 0;

			int row = 0;

			while (row < nCombinations) {

				for (int i = 0; i < nRepetitions; i++) {
					int [] arrComb = liComb.get(row);
					arrComb[col]=arr[indexArr];
					row++;
				}

				indexArr++;
				if(indexArr==arr.length) {
					indexArr=0;
				}
			}
		}

		return liComb;
	}


	public static void main(String[] args) {
		List<int[]> li = new ArrayList<>();

//		int [] a1 = {0,1};
//		int [] a2 = {0};
//		int [] a3 = {1,2};
		int [] a1 = {0,1,2};
		int [] a2 = {3,4};
		int [] a3 = {5,6,7};
		int [] a4 = {8};

		li.add(a1);
		li.add(a2);
		li.add(a3);
		li.add(a4);

		List<int[]> liComb = getCombinations(li);

		for (int[] arr : liComb) {
			System.out.println(ArrayUtils.toString(arr));
		}





	}

//	public static void main(String[] args) {
//
//		int sizeList = 9;
//
//		int sumCombinations=0;
//		for (int i = 4; i < sizeList+1; i++) {
//
//			List<int[]> li = CombinationGenerator.getAllOutOf(sizeList, i);
//
//			sumCombinations += li.size();
//			for (int j = 0; j < li.size(); j++) {
//				int [] a = li.get(j);
//				System.out.println(ArrayUtils.toString(a));
//
//			}
//		}
//
//		System.out.println("sum combinations " + sumCombinations);
//
//
//
//	}

}
