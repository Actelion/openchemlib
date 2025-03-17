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
 */
package com.actelion.research.calc.combinatorics;

import com.actelion.research.util.ListUtils;
import com.actelion.research.util.StringFunctions;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * CombinationGenerator
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

	/**
	 *
	 * @param arrSizeClass, Each value in the array stands for the number of objects.
	 * @return
	 */
	public static List<List<Integer>> getCombinations(int [] arrSizeClass) {

		long capacity = 1;

		for (int sizeClass : arrSizeClass) {
			capacity *= sizeClass;
		}

		if(capacity>Integer.MAX_VALUE){
			throw new RuntimeException("Number of combinations " + capacity + " above Integer.MAX_VALUE!");
		}

		List<List<Integer>> li = new ArrayList<>((int)capacity);

		List<List<Integer>> liTmp = new ArrayList<>((int)capacity);

		int index=0;

		for (int i = 0; i < arrSizeClass[0]; i++) {
			List<Integer> liCombi = new ArrayList<>();
			liCombi.add(i);
			li.add(liCombi);
		}

		index++;

		while (index<arrSizeClass.length) {

			liTmp.clear();

			for (List<Integer> liCombi : li) {
				for (int i = 0; i < arrSizeClass[index]; i++) {
					List<Integer> liCombiNew = new ArrayList<>(liCombi);
					liCombiNew.add(i);
					liTmp.add(liCombiNew);
				}
			}
			li.clear();
			li.addAll(liTmp);

			index++;
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
	

	private static void swap(int[] input, int a, int b) {
		    int tmp = input[a];
		    input[a] = input[b];
		    input[b] = tmp;
	}

	/**
	 *
	 * @param elements
	 * @param n
	 * @return
	 */
	public static List<int[]> getPermutations(int[] elements, int n){
		List<int[]> permutations = new ArrayList<int[]>();
		int[] indexes = new int[n];
		int i = 0;
		permutations.add(Arrays.copyOf(elements, elements.length));
		while (i < n) {
		    if (indexes[i] < i) {
		        swap(elements, i % 2 == 0 ?  0: indexes[i], i);
		        permutations.add(Arrays.copyOf(elements, elements.length));
		        indexes[i]++;
		        i = 0;
		    }
		    else {
		        indexes[i] = 0;
		        i++;
		    }
		}
		return permutations;
	}
	
	/**
	 * https://en.wikipedia.org/wiki/Cartesian_product
	 * generates all possible combinations of elements from a list of lists
	 * @param <T>
	 * @param lists
	 * @return
	 */
	
	public static <T> List<List<T>> cartesianProduct(List<List<T>> lists) {
	    List<List<T>> resultLists = new ArrayList<List<T>>();
	    if (lists.size() == 0) {
	        resultLists.add(new ArrayList<T>());
	        return resultLists;
	    } else {
	        List<T> firstList = lists.get(0);
	        List<List<T>> remainingLists = cartesianProduct(lists.subList(1, lists.size()));
	        for (T condition : firstList) {
	            for (List<T> remainingList : remainingLists) {
	                ArrayList<T> resultList = new ArrayList<T>();
	                resultList.add(condition);
	                resultList.addAll(remainingList);
	                resultLists.add(resultList);
	            }
	        }
	    }
	    return resultLists;
	}

	public static BigInteger getFactorial (int n) {
		BigInteger fact = BigInteger.ONE;
		for (int i = n; i > 1; i--) {
			fact = fact.multiply (new BigInteger (Integer.toString (i)));
		}
		return fact;
	}

	/**
	 * Calculate binomial coefficient or
	 * n choose k
	 *
	 * @param n
	 * @param k
	 * @return
	 */
	public static BigInteger getBinomialCoefficient(int n, int k){
		BigInteger nFac = getFactorial(n);
		BigInteger kFac = getFactorial(k);
		BigInteger nMinus_k_Fac = getFactorial(n-k);
		BigInteger dev = nMinus_k_Fac.multiply(kFac);
		BigInteger bc = nFac.divide(dev);
		return bc;
	}

	public static void main(String[] args) {
		exampleGetAllOutOf2();

	}
	public static void examplePermutations() {
		int [] a = {4,7,9};

		List<int []> liComb = getPermutations(a, 3);
		for(int [] p : liComb) {
			System.out.println(Arrays.toString(p));
		}
	}
	public static void exampleCartesianProduct() {

		List<List<Integer>> lists = new ArrayList<>();

		List<Integer> l0 = new ArrayList<>();

		l0.add(1);
		l0.add(2);
		l0.add(3);

		List<Integer> l1 = new ArrayList<>(l0);
		List<Integer> l2 = new ArrayList<>(l0);

		lists.add(l0);
		lists.add(l1);
		lists.add(l2);

		List<List<Integer>> liComb = cartesianProduct(lists);
		for(List li : liComb) {
			System.out.println(ListUtils.toStringInteger(li));
		}
	}
	public static void exampleCombinations() {

		int [] a = {1,2,3};

		List<List<Integer>> liComb = getCombinations(a);
		for(List li : liComb) {
			System.out.println(ListUtils.toStringInteger(li));
		}
	}

	public static void exampleGetAllOutOf() {

		int object = 4;
		int sampleSize = 3;
		List<int[]>  liComb = getAllOutOf(object, sampleSize);
		for(int [] a : liComb) {
			System.out.println(StringFunctions.toString(a, ","));
		}
	}
	public static void exampleGetAllOutOf2() {

		int object = 4;
		for (int sampleSize = 1; sampleSize < object; sampleSize++) {
			List<int[]>  liComb = getAllOutOf(object, sampleSize);
			for(int [] a : liComb) {
				System.out.println(StringFunctions.toString(a, ","));
			}
		}
	}
}
