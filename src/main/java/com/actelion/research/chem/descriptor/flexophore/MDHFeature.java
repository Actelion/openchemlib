package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.calc.ArrayUtilsCalc;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

/**
 * 
 * MDHFeature
 * @author Modest von Korff
 * @version 1.0
 * 22 Sep 2008 MvK: Start implementation
 */
public class MDHFeature implements Comparable<MDHFeature> {

	private static final NumberFormat NF = new DecimalFormat("0.0000");
	
	private int index;

	private int counter1;

	private int counter2;

	public MolDistHistViz mdhvFeature;
	
	private HashSet<Integer> hsMatchingIndex1;
	
	private HashSet<Integer> hsMatchingIndex2;

	public MDHFeature(MolDistHistViz mdhvFeature, int index) {
		init();
		this.mdhvFeature = mdhvFeature;
		this.index = index;
	}

	public MDHFeature(MolDistHistViz mdhvFeature) {
		init();
		this.mdhvFeature = mdhvFeature;
	}

	private void init() {
		
		hsMatchingIndex1 = new HashSet<Integer>();
		
		hsMatchingIndex2 = new HashSet<Integer>();
	}
	
	public int compareTo(MDHFeature m) {

		if (counter1 > m.counter1)
			return 1;
		if (counter1 < m.counter1)
			return -1;
		else
			return 0;
	}

	public int getIndex() {
		return index;
	}

	public int getCounter1() {
		return counter1;
	}

	public int getCounter2() {
		return counter2;
	}

	public MolDistHistViz getMdhvFeature() {
		return mdhvFeature;
	}

	public void increaseClass1(int index) {
		counter1++;
		
		if(!hsMatchingIndex1.add(index))
			throw new RuntimeException("Index " + index + " is already in active list.");
	}

	public void increaseClass2(int index) {
		counter2++;
		
		if(!hsMatchingIndex2.add(index))
			throw new RuntimeException("Index " + index + " is already in inactive list.");
	}

	public double getRatioClass1Class2(){
		
		if(counter2==0){
			return Double.POSITIVE_INFINITY;
		}
		
		return ((double) counter1) / ((double) counter2);
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		
		sb.append(index);
		
		double ratioRisk = -1;
		if(counter1==0 && counter2==0) {
			ratioRisk = -1;
		} else if(counter2==0) {
			ratioRisk = counter1;
		} else {
			ratioRisk = getRatioClass1Class2();
		}
		
		sb.append(" risk\t" + counter1 + ", no risk\t" + counter2 + ", ratio\t" + NF.format(ratioRisk));
		
		// sb.append(mdhvFeature.toString());

		return sb.toString();
	}
	
	public String toStringListClass1() {
		List<Integer> li = new ArrayList<Integer>(hsMatchingIndex1);
		Collections.sort(li);
		return ArrayUtilsCalc.toString(li);
	}
	
	public String toStringListClass2() {
		List<Integer> li = new ArrayList<Integer>(hsMatchingIndex2);
		Collections.sort(li);
		return ArrayUtilsCalc.toString(li);
	}

}
