package com.actelion.research.util.datamodel;

import com.actelion.research.calc.Matrix;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * ModelXYCrossValidation
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
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
