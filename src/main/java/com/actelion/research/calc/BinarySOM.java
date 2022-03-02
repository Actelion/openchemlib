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
 * @author Thomas Sander
 */

package com.actelion.research.calc;

import com.actelion.research.chem.SSSearcherWithIndex;

import java.awt.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class BinarySOM extends SelfOrganizedMap {
    private static final int KEY_COUNT = SSSearcherWithIndex.getNoOfKeys();
    private static final int MASK_COUNT = (SSSearcherWithIndex.getNoOfKeys()+31)/32;
	private double[]	mKeyFrequency;
	private int			mMaxKeyCount;
	private int[]		mMask;
	private int[][]		mKeyMaskMap;
	private ArrayList	mRandomizedKeyIndexList;
	private byte[]		mBitCount;
	private Random		mRandom;

	public BinarySOM() {
			// constructor to be used if SOM interna are read from a SOM file with read()
		init();
		}
	
	public BinarySOM(int nx, int ny, int mode) {
		super(nx, ny, mode);
		init();
		}

	private void init() {
		mRandom = new Random();

		mMask = new int[32];
		mMask[0] = 0x80000000;
		for (int i=0; i<31; i++)
		    mMask[i+1] = mMask[i] >>> 1;

		mRandomizedKeyIndexList = new ArrayList();
		for (int i=0; i<KEY_COUNT; i++)
		    mRandomizedKeyIndexList.add(new Integer(i));

		mBitCount = new byte[0x10000];
		for (int i=0; i<0x10000; i++)
		    mBitCount[i] = (byte)Integer.bitCount(i);
		}
	
	protected void initializeNormalization() {
		startProgress("Calculating key frequencies...", 0, mController.getInputVectorCount());

		mKeyFrequency = new double[KEY_COUNT];
		
		for (int row=0; row<mController.getInputVectorCount(); row++) {
			if (threadMustDie())
				break;
			updateProgress(row);

			int[] keyList = (int[])mController.getInputVector(row);

			for (int i=0; i<MASK_COUNT; i++)
			    for (int j=0; j<32; j++)
			        if ((keyList[i] & mMask[j]) != 0)
			            mKeyFrequency[i*32+j]++;
			}
		if (!threadMustDie())
			for (int i=0; i<KEY_COUNT; i++)
			    mKeyFrequency[i] /= (double)mController.getInputVectorCount();
		}

	public void write(BufferedWriter writer) throws IOException {
		super.write(writer);

		writer.write("<keyFrequency=\""+VectorSOM.doubleArrayToString(mKeyFrequency)+"\">");
		writer.newLine();
		}

	public void read(BufferedReader reader) throws Exception {
		super.read(reader);

		String theLine = reader.readLine();
		boolean error = !theLine.startsWith("<keyFrequency=");
		if (!error) {
		    mKeyFrequency = VectorSOM.stringToDoubleArray(extractValue(theLine));
			}

		if (error)
			throw new IOException("Invalid SOM file format");
		}

	protected String referenceVectorToString(int x, int y) {
		return SSSearcherWithIndex.getHexStringFromIndex((int[])mReferenceVector[x][y]);
		}

	protected void setReferenceVector(int x, int y, String ref) throws Exception {
		mReferenceVector[x][y] = SSSearcherWithIndex.getIndexFromHexString(ref);
		}

	public double getDissimilarity(Object vector1, Object vector2) {
	    int[] k1 = (int[])vector1;
	    int[] k2 = (int[])vector2;
	    int sharedKeys = 0;
	    int allKeys = 0;
	    for (int i=0; i<MASK_COUNT; i++) {
	        int sk = k1[i] & k2[i];
	        int ak = k1[i] | k2[i];
	        sharedKeys += mBitCount[0xFFFF & sk] + mBitCount[sk >>> 16];
	        allKeys += mBitCount[0xFFFF & ak] + mBitCount[ak >>> 16];
	    	}
	    return 1.0 - (double)sharedKeys/(double)allKeys;

//	    return 1.0 - SSSearcherWithIndex.getSimilarityTanimoto((int[])vector1, (int[])vector2);
		}

	protected void updateReference(Object inputVector, Object referenceVector, double influence) {
	    int[] mask = mKeyMaskMap[(int)influence];

	    int[] inputKeyList = (int[])inputVector;
	    int[] referenceKeyList = (int[])referenceVector;
		for (int i=0; i<MASK_COUNT; i++) {
		    referenceKeyList[i] ^= mask[i] & (inputKeyList[i] ^ referenceKeyList[i]);
// alternative with same result:
//		    referenceKeyList[i] &= ~mask[i];
//		    referenceKeyList[i] |= (mask[i] & inputKeyList[i]);
			}
		}

	protected Object getRandomVector() {
	    int[] keyList = new int[MASK_COUNT];
	    for (int i=0; i<MASK_COUNT; i++)
	        for (int j=0; j<32; j++)
	            if (mRandom.nextDouble() < mKeyFrequency[32*i+j])
	                keyList[i] |= mMask[j];
	    return keyList;
		}

	protected Object getMeanVector(Object vector1, Object vector2) {
			// must be overridden if input vectors aren't double arrays
		int[] v1 = (int[])vector1;
		int[] v2 = (int[])vector2;
		int[] mv = new int[v1.length];

        for (int i=0; i<MASK_COUNT; i++) {
            mv[i] = v1[i];
            for (int j=0; j<32; j++) {
                if ((v1[i] & mMask[j]) != (v2[i] & mMask[j])
                 && Math.random() < 0.5) {
                    mv[i] &= ~mMask[j];
                    mv[i] |= (v2[i] & mMask[j]);
                    }
                }
            }

		return mv;
		}

	public Object normalizeVector(Object vector) {
		return vector;
		}

	protected void calculateInfluences(double time) {
	    	// modify mInfluence[][] to contain noOfBits to adapt
	    	// and allocate empty mask for used bitCounts
	    super.calculateInfluences(time);
	    mMaxKeyCount = 0;
		mKeyMaskMap = new int[KEY_COUNT][];
		for (int i=0; i<mInfluence.length; i++) {
		    for (int j=0; j<mInfluence[i].length; j++) {
		        if (mInfluence[i][j] != 0.0) {
		            int keyCount = (int)(KEY_COUNT * mInfluence[i][j] + 0.5);
		            if (mMaxKeyCount < keyCount)
		                mMaxKeyCount = keyCount;
	                mInfluence[i][j] = keyCount;
		            if (mKeyMaskMap[keyCount] == null)
		                mKeyMaskMap[keyCount] = new int[MASK_COUNT];
		    		}
		    	}
			}
		}

	protected void applyInfluences(Object inputVector, Point location) {
	    randomizeKeyMasks();
		super.applyInfluences(inputVector, location);
		}

	private void randomizeKeyMasks() {
	    int[] mask = new int[MASK_COUNT];
	    Collections.shuffle(mRandomizedKeyIndexList);
	    for (int keyCount=1; keyCount<=mMaxKeyCount; keyCount++) {
            int keyIndex = ((Integer)mRandomizedKeyIndexList.get(keyCount-1)).intValue();
            mask[keyIndex >> 5] |= mMask[keyIndex & 0x1F];
            if (mKeyMaskMap[keyCount] != null)
                for (int i=0; i<MASK_COUNT; i++)
                    mKeyMaskMap[keyCount][i] = mask[i];
	    	}
		}
	}