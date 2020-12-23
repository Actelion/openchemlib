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
*/

package com.actelion.research.chem;

import com.actelion.research.calc.DataProcessor;

import java.util.Arrays;
import java.util.Comparator;

public class DiversitySelector extends DataProcessor {
	private int			mNoOfFeatures,mExistingSetCount;
	private long[][] 	mFeatureList;
	private boolean		mAddingToExistingSet;
	private double[]	mCentroidVector;

	public DiversitySelector() {
		super();
		}

	public void initializeExistingSet(int noOfKeys) {
		mFeatureList = new long[1][];
		mCentroidVector = new double[noOfKeys];
		mExistingSetCount = 0;
	    }

	public void addToExistingSet(long[] featureList) {
		mFeatureList[0] = featureList;
		DiversitySelectorRecord record = new DiversitySelectorRecord(0);
		record.addToCentroidVector();
		mAddingToExistingSet = true;
		mExistingSetCount++;
    	}


	public void setExistingSet(long[][] featureList) {
		// create centroid vector based on already selected set of compounds
		mCentroidVector = new double[64*featureList[0].length];
		mFeatureList = featureList;
		for (int compound=0; compound<featureList.length; compound++) {
			DiversitySelectorRecord record = new DiversitySelectorRecord(compound);
			record.addToCentroidVector();
			}

		mAddingToExistingSet = true;
		mExistingSetCount = featureList.length;
		}


	public int[] select(long[][] featureList, int compoundsToSelect) {
		int compoundsAvailable = featureList.length;
		mNoOfFeatures = 64*featureList[0].length;
		mFeatureList = featureList;

		if (compoundsToSelect > compoundsAvailable)
			compoundsToSelect = compoundsAvailable;

		startProgress("Creating Key Lists...", 0, compoundsAvailable);

		DiversitySelectorRecord[] recordList = new DiversitySelectorRecord[compoundsAvailable];
		for (int compound=0; compound<compoundsAvailable; compound++) {
			recordList[compound] = new DiversitySelectorRecord(compound);
			if ((compound & 255) == 255) {
				if (threadMustDie()) {
				    stopProgress("Selection cancelled");
					return null;
					}
				updateProgress(compound);
				}
			}

	    startProgress("Locating Starting Compound...", 0, compoundsAvailable);

		if (!mAddingToExistingSet) {
			// find most similar compound as starting point
			mCentroidVector = new double[mNoOfFeatures];
			for (int compound=0; compound<compoundsAvailable; compound++) {
				recordList[compound].addToCentroidVector();
				if ((compound & 255) == 255) {
					if (threadMustDie()) {
					    stopProgress("Selection cancelled");
						return null;
						}
					updateProgress(compound/2);
					}
				}
			double maxDotProduct = 0.0;
			int maxCompoundIndex = 0;
			for (int compound=0; compound<compoundsAvailable; compound++) {
				// dot product must be based on the centroid vector of the complete set minus the compound
				// under investigation.
				double dotProduct = 0.0;
				for (int keyIndex=0; keyIndex<recordList[compound].mKeyList.length; keyIndex++) {
					int key = recordList[compound].mKeyList[keyIndex];
					dotProduct += (mCentroidVector[key] - recordList[compound].mWeight) * recordList[compound].mWeight;
					}
				if (maxDotProduct < dotProduct) {
					maxDotProduct = dotProduct;
					maxCompoundIndex = compound;
					}
				if ((compound & 255) == 255) {
					if (threadMustDie()) {
					    stopProgress("Selection cancelled");
						return null;
						}
					updateProgress((compoundsAvailable+compound)/2);
					}
				}
			DiversitySelectorRecord startCompound = recordList[maxCompoundIndex];
			recordList[maxCompoundIndex] = recordList[0];
			recordList[0] = startCompound;

			mCentroidVector = new double[mNoOfFeatures];
			startCompound.addToCentroidVector();
			}

/*
try {
BufferedWriter writer = new BufferedWriter(new FileWriter("d:\\test.txt"));
writer.write("run\tposition\tindex\tdotProduct\n");
*/

	    startProgress("Selecting Compounds...", 0, compoundsToSelect);

		int selectionCycle = 0;
		for (int compound=(mAddingToExistingSet)?0:1; compound<compoundsToSelect; compound++) {
			// set noOfCompoundsToSort larger than typical maximum of sort position change considering
			// no of compounds already in centroid vector and size of remaining compound set
			int noOfCompoundsToSort = (int)(10.0*(mExistingSetCount+compoundsAvailable)/(mExistingSetCount+compound));

			// increase every forth sorting noOfCompoundsToSort by factor of 4
			// increase every 16th sorting noOfCompoundsToSort by factor of 16
			// increase every 64th sorting noOfCompoundsToSort by factor of 64 and so on
			int mask = 0x00000003;
			while (noOfCompoundsToSort < compoundsAvailable - compound) {
				if ((selectionCycle & mask) != 0)
					break;
				mask = (mask << 2) | 0x00000003;
				noOfCompoundsToSort *= 4;
				}

			int lastCompoundToConsider = compound + noOfCompoundsToSort;
			if (lastCompoundToConsider > compoundsAvailable)
				lastCompoundToConsider = compoundsAvailable;

			for (int i=compound; i<lastCompoundToConsider; i++)
				recordList[i].calculateDotProduct();

			Arrays.sort(recordList, compound, lastCompoundToConsider,
                        new DiversitySelectorComparator<>());
			recordList[compound].addToCentroidVector();
			selectionCycle++;
/*
for (int i=compound; i<compoundsAvailable; i++)
if (recordList[i].mCompoundIndex < 20)
writer.write(compound+"\t"+i+"\t"+(1+recordList[i].mCompoundIndex)+"\t"+recordList[i].mDotProduct+"\n");
*/
			if (threadMustDie()) {
			    stopProgress("Selection cancelled");
				return null;
				}
			updateProgress(compound);
			}

/*
writer.close();
} catch (Exception e) {
e.printStackTrace();
}
*/
		int[] selected = new int[compoundsToSelect];
		for (int compound=0; compound<compoundsToSelect; compound++)
			selected[compound] = recordList[compound].mCompoundIndex;

	    stopProgress("Compound Selection Done");

			// the selected array contains compound indices of selected compounds only
			// sorted by the selection order (first selected is first in array)
		return selected;
		}


	protected class DiversitySelectorRecord {
		private int mCompoundIndex;
		private int[] mKeyList;
		private double mWeight;
		public double mDotProduct;

		private DiversitySelectorRecord(int index) {
			int count = 0;
			for (int feature=0; feature<mNoOfFeatures; feature++)
				if ((mFeatureList[index][feature/64] & (1 << (63-feature%64))) != 0)
					count++;
			mKeyList = new int[count];
			count = 0;
			for (int feature=0; feature<mNoOfFeatures; feature++)
				if ((mFeatureList[index][feature/64] & (1 << (63-feature%64))) != 0)
					mKeyList[count++] = feature;

			mWeight = 1.0 / Math.sqrt((double)count);
			mCompoundIndex = index;
			}

		private void addToCentroidVector() {
			for (int i=0; i<mKeyList.length; i++)
				mCentroidVector[mKeyList[i]] += mWeight;
			}

		private void calculateDotProduct() {
			mDotProduct = 0.0;
			for (int i=0; i<mKeyList.length; i++)
				mDotProduct += mWeight * mCentroidVector[mKeyList[i]];
			}
		}

    class DiversitySelectorComparator<T> implements Comparator<T> {
        public int compare(T o1, T o2) {
            double d1 = ((DiversitySelectorRecord)o1).mDotProduct;
            double d2 = ((DiversitySelectorRecord)o2).mDotProduct;
            return (d1 < d2) ? -1 : (d1 == d2) ? 0 : 1;
            }
        }
    }
