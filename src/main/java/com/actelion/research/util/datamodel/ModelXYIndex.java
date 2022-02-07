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
 * 
 * 
 * ModelDataXY
 * @author Modest von Korff
 * @version 1.0
 * Aug 4, 2011 MvK: Start implementation
 */
public class ModelXYIndex extends ModelXY {
			
	
	// Index of the X row origin in the original file.
	// Becomes important when rows in X and Y were removed.
	public List<Integer> liIndex;
		
	public ModelXYIndex() {
		
	}
	public ModelXYIndex(Matrix x, Matrix y) {
		super(x,y);
	}

	public ModelXYIndex(int rows, int colsX, int colsY) {
		super(rows, colsX, colsY);
	}


	/**
	 * Deep copy constructor
	 * @param dataXY
	 */
	public ModelXYIndex(ModelXYIndex dataXY) {
		super(dataXY);
		
		liIndex = new ArrayList<Integer>(dataXY.X.rows());
				
		for (int index : dataXY.liIndex) {
			liIndex.add(index);
		}
	}


	
	public ModelXYIndex(ModelXY dataXY) {
		super(dataXY);
	}

	/**
	 *
	 * @param indexStart inclusive
	 * @param indexEnd exclusive
	 * @return
	 */
	public ModelXYIndex sub(int indexStart, int indexEnd){

		List<Integer> liIndexRow = new ArrayList<>(indexEnd-indexStart);

		for (int i = indexStart; i < indexEnd; i++) {
			liIndexRow.add(i);
		}

		return sub(liIndexRow);
	}

	
	public ModelXYIndex sub(List<Integer> liIndexRow){

		ModelXYIndex modelDataXY = new ModelXYIndex();
		
		modelDataXY.X = X.getSubMatrix(liIndexRow);
		
		modelDataXY.Y = Y.getSubMatrix(liIndexRow);
		
		modelDataXY.liIndex = new ArrayList<>(liIndexRow);

		return modelDataXY;
	}

	/**
	 * Rows are copied from this into target.
	 * @param liIndexRow
	 * @param modelXYIndexTarget
	 */
	public void sub(List<Integer> liIndexRow, ModelXYIndex modelXYIndexTarget){

		if((modelXYIndexTarget.X.rows()!=liIndexRow.size()) || (modelXYIndexTarget.Y.rows()!=liIndexRow.size())) {
			throw new RuntimeException("Wrong number of rows!");
		}

		for (int i = 0; i < liIndexRow.size(); i++) {

			int index = liIndexRow.get(i);

			double [] a = X.getRow(index);

			modelXYIndexTarget.X.setRow(i, a);

			double [] b = Y.getRow(index);

			modelXYIndexTarget.Y.setRow(i, b);
		}
	}

	public ModelXYIndex getVariableSelectionOnX(IntVec iv) {

		int colsNew = iv.getBitsSet();

		Matrix XNew = new Matrix(X.rows(), colsNew);

		int n = iv.sizeBits();

		int colNew = 0;
		for (int i = 0; i < n; i++) {

			if(iv.isBitSet(i)){
				XNew.copyColumn(X, i, colNew++);
			}
		}

		ModelXYIndex modelXYIndexSelectedCols = new ModelXYIndex();

		modelXYIndexSelectedCols.X = XNew;
		modelXYIndexSelectedCols.Y = new Matrix(Y);

		if(liIndex!=null) {
			modelXYIndexSelectedCols.liIndex = new ArrayList<>(X.rows());

			for (int index : liIndex) {
				modelXYIndexSelectedCols.liIndex.add(index);
			}
		}

		return modelXYIndexSelectedCols;
	}

	@Override
	public ModelXYIndex getDeepClone() {
		return new ModelXYIndex(this);
	}


	/**
	 *
	 * @return Deep copy.
	 */
	public List<XYIndex> getAsListWithIndex(){

		int rows = X.rows();

		List<XYIndex> li = new ArrayList<>(rows);

		for (int i = 0; i < rows; i++) {
			XYIndex XY = new XYIndex(X.getRowCopy(i), Y.getRowCopy(i));

			if(liIndex!=null){
				XY.index = liIndex.get(i);
			}
			li.add(XY);
		}

		return li;
	}

	public ModelXYIndex getSortedByY(int colY){

		List<XYIndex>  li = getAsListWithIndex();

		Collections.sort(li, XY.getComparatorY(colY));

		double [][] arrX = new double[li.size()][];
		double [][] arrY = new double[li.size()][];

		List<Integer> liIndex = new ArrayList<>(li.size());
		for (int i = 0; i < li.size(); i++) {
			arrX[i]=li.get(i).x;
			arrY[i]=li.get(i).y;
			liIndex.add(li.get(i).index);
		}

		ModelXYIndex modelSorted = new ModelXYIndex();

		modelSorted.X = new Matrix(arrX);
		modelSorted.Y = new Matrix(arrY);
		modelSorted.liIndex=liIndex;

		return modelSorted;
	}






}