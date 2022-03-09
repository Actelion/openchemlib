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

package com.actelion.research.util.datamodel;

import com.actelion.research.calc.Matrix;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * ModelXYCrossValidation
 * @author Modest von Korff
 * @version 1.0
 * Aug 14, 2015 MvK Start implementation

 */
@Deprecated // has to be deleted
public class ModelXYCrossValidation extends ModelXY {

	
	private double fractionLeaveOut;
	
	private int nTest;
	
	private int nTrain;
	
	private Matrix Xtest;
	
	private Matrix Ytest;
	
	private Matrix Xtrain;
	
	private Matrix Ytrain;
	
	private List<Integer> liIndex;
	
	
	/**
	 * 
	 */
	public ModelXYCrossValidation(ModelXY modelXY) {
		super(modelXY);
	}

	/**
	 * @param fractionLeaveOut the fractionLeaveOut to set
	 */
	public void setFractionLeaveOut(double fractionLeaveOut) {
		this.fractionLeaveOut = fractionLeaveOut;
		init();
	}
	
	private void init() {
		
		nTest = (int)(X.rows() * fractionLeaveOut);
		
		nTrain = X.rows()-nTest;
		
		Xtest = new Matrix(nTest, X.cols());
		
		Ytest = new Matrix(nTest, Y.cols());
		
		Xtrain = new Matrix(nTrain, X.cols());
		
		Ytrain = new Matrix(nTrain, Y.cols());
		
		liIndex = new ArrayList<Integer>(X.rows());
		for (int i = 0; i < X.rows(); i++) {
			liIndex.add(i);
		}
	}
	
	public void next(){
		
		Collections.shuffle(liIndex);
		
		List<Integer> liIndexTest = liIndex.subList(0, nTest); 
		
		copy(X, liIndexTest, Xtest);
		
		copy(Y, liIndexTest, Ytest);
		
		List<Integer> liIndexTrain = liIndex.subList(0, nTrain); 
		
		copy(X, liIndexTrain, Xtrain);
		
		copy(Y, liIndexTrain, Ytrain);

	}
	
	private static void copy(Matrix maSource, List<Integer> liIndex, Matrix maDest) {
		
		for (int i = 0; i < liIndex.size(); i++) {
			maDest.setRow(i, maSource.getRow(liIndex.get(i)));
		}
		
	}

	/**
	 * @return the xtest
	 */
	public Matrix getXtest() {
		return Xtest;
	}

	/**
	 * @return the ytest
	 */
	public Matrix getYtest() {
		return Ytest;
	}

	/**
	 * @return the xtrain
	 */
	public Matrix getXtrain() {
		return Xtrain;
	}

	/**
	 * @return the ytrain
	 */
	public Matrix getYtrain() {
		return Ytrain;
	}
	
	


}
