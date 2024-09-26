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

package com.actelion.research.calc;

import com.actelion.research.util.DoubleVec;
import com.actelion.research.util.convert.String2DoubleArray;
import com.actelion.research.util.datamodel.IntegerDouble;
import com.actelion.research.util.datamodel.ScorePoint;

import java.awt.*;
import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.List;
import java.util.*;

public class Matrix {

    public static final double TINY04 = 0.0001;
    public static final double TINY08 = 0.00000001;
    public static final double TINY16 = Math.pow(10, -16);
    public static String OUT_SEPARATOR_COL = "\t";
    
    public static String OUT_SEPARATOR_ROW = "\n";
    
    private final static String FORMAT = "0.############";

    public static final double TINY = 10e-12;

    private double [][] data;

    private int identifier;

    public Matrix() {
        data = new double[1][1];
    }

    public Matrix(int rows, int cols) {
        data = new double[rows][cols];
    }

    public Matrix(Matrix ma) {
        this(ma.getArray());
    }

    /**
     *
     * @param row if true the matrix has one row. If false the matrix has one column.
     * @param arr
     */
    public Matrix(boolean row, double [] arr) {
        if (row) {
            data = new double[1][];
            data[0]=arr;
        }
        else {
            data = new double[arr.length][1];
            for (int i = 0; i < getRowDim(); i++) {
                data[i][0] = arr[i];
            }
        }
    }
    
    public Matrix(boolean row, byte [] arr) {
        if (row) {
            data = new double[1][arr.length];
            for (int i = 0; i < arr.length; i++) {
                data[0][i]=arr[i];
            }
        }
        else {
            data = new double[arr.length][1];
            for (int i = 0; i < getRowDim(); i++) {
                data[i][0] = arr[i];
            }
        }
    }

    public Matrix(boolean row, int [] dArray) {
        if (row) {
            data = new double[1][dArray.length];
            for (int jj = 0; jj < getColDim(); jj++) {
                data[0][jj] = dArray[jj];
            }
        }
        else {
            data = new double[dArray.length][1];
            for (int ii = 0; ii < getRowDim(); ii++) {
                data[ii][0] = dArray[ii];
            }
        }
    }
    /**
     * Deep copy
     * @param arrArr
     */
    public Matrix(double [][]arrArr) {
        data = new double[arrArr.length][];
        int rows = arrArr.length;
        int cols = arrArr[0].length;
        for (int i = 0; i < rows; i++) {
        	double [] arr = new double [cols];
        	System.arraycopy(arrArr[i], 0, arr, 0, cols);
            data[i] = arr;
        }
    }

    public Matrix(double [][]arrArr, boolean flat) {
        if(!flat){
            throw new RuntimeException("Only flat constructor!");
        }
        data = arrArr;
    }

    public Matrix(float [][]arrArr) {
        data = new double[arrArr.length][arrArr[0].length];
        int rows = arrArr.length;
        int cols = arrArr[0].length;
        for (int i = 0; i < rows; i++) {
        	double [] arr = new double [cols];
        	for (int j = 0; j < arr.length; j++) {
        		data[i][j] = arrArr[i][j];
			}
        }
    }

    public Matrix(int [][]arrArr) {
        
        data = new double[arrArr.length][];
        
        int rows = arrArr.length;
        
        int cols = arrArr[0].length;
        
        // 29.06.2012 System.arrayCopy gave ArrayStoreException.
        for (int i = 0; i < rows; i++) {
        	double [] arr = new double [cols];
        	
        	for (int j = 0; j < cols; j++) {
        		arr[j]=arrArr[i][j];
			}
        	
            
            data[i] = arr;
        }
    }
    
    public Matrix(byte [][]arrArr) {
    	
        data = new double[arrArr.length][];
        
        int rows = arrArr.length;
        
        int cols = arrArr[0].length;
        
        for (int i = 0; i < rows; i++) {
        	double [] arr = new double [cols];
        	
        	System.arraycopy(arrArr[i], 0, arr, 0, cols);
            
            data[i] = arr;
        }
    }

    /**
     *
     * @param vecDoubleVec Vector with DoubleVec
     */
    public Matrix(List<DoubleVec> vecDoubleVec) {
        DoubleVec dv = (DoubleVec) vecDoubleVec.get(0);
        data = new double[vecDoubleVec.size()][dv.size()];
        for (int ii = 0; ii < getRowDim(); ii++) {
            dv = (DoubleVec) vecDoubleVec.get(ii);
            for (int jj = 0; jj < getColDim(); jj++) {
                data[ii][jj] = dv.get(jj);
            }
        }
    }

    public Matrix(boolean bRow, List<Double> liDoubles) {
        if (bRow) {
            data = new double[1][liDoubles.size()];
            for (int i = 0; i < getColDim(); i++) {
                data[0][i] = liDoubles.get(i);
            }
        }
        else {
            data = new double[liDoubles.size()][1];
            for (int i = 0; i < getRowDim(); i++) {
                data[i][0] = liDoubles.get(i);
            }
        }
    }


    public Matrix add(double v) {

        int r = rows();
        int c = cols();

        Matrix A = new Matrix(r,c);

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                double s = get(i,j) + v;
                A.set(i,j,s);
            }
        }
        return A;
    }

    /**
     *
     * @param arr adds a row at the end of the matrix.
     */
    public void addRow(double [] arr) {
        if(getColDim() != arr.length) {
            throw new RuntimeException("Matrices have wrong dimensions.");
        }

        resize(getRowDim() + 1, getColDim());

        for (int ii = 0; ii < arr.length; ii++)
          data[getRowDim() - 1][ii] = arr[ii];
    }

    public void addRow(int [] arr) {
        if(getColDim() != arr.length) {
            throw new RuntimeException("Matrices have wrong dimensions.");
        }

        resize(getRowDim() + 1, getColDim());

        for (int ii = 0; ii < arr.length; ii++)
          data[getRowDim() - 1][ii] = arr[ii];
    }
    /**
     * @param row where the data are added on
     * @param ma2 matrix from where the data are taken
     * @param row2 data row
     */
    public void add2Row(int row, Matrix ma2, int row2) {
        for (int ii = 0; ii < getColDim(); ii++)
          data[row][ii] += ma2.data[row2][ii];
    }

    public void addToElement(int i, int j, double value) {
        data[i][j]+=value;

    }

    /**
     * The value in the col from the input matrix is
     * added to all values in the corresponding col in the matrix.
     * @param maRow matrix with one row.
     */
    public Matrix add2CompleteCol(Matrix maRow) {
        Matrix m = new Matrix(this);
        for (int ii = 0; ii < getRowDim(); ii++)
            for (int jj = 0; jj < getColDim(); jj++)
                m.data[ii][jj] += maRow.data[0][jj];

        return m;
    }

    public void assignCol(int iCol, Matrix ma) {
        if(getRowDim() != ma.getRowDim()) {
            throw new RuntimeException("Matrices have wrong dimensions.");
        }
        for (int ii = 0; ii < data.length; ii++) {
          data[ii][iCol] = ma.data[ii][0];
        }
    }

    public void assignCol(int iCol, Matrix ma, int iColMa) {
        if(getRowDim() != ma.getRowDim()) {
            throw new RuntimeException("Matrices have wrong dimensions.");
        }
        for (int ii = 0; ii < data.length; ii++) {
          data[ii][iCol] = ma.data[ii][iColMa];
        }
    }

    public void assignRow(int iRow, Vector<Double> vecRow) {
        if(getColDim() != vecRow.size()) {
            throw new RuntimeException("Matrix and Vector have wrong dimensions.");
        }
        for (int jj = 0; jj < data[0].length; jj++) {
          data[iRow][jj] = vecRow.get(jj);
        }
    }

    public void assignRow(int iRow, double [] arr) {
        if(getColDim() != arr.length) {
            throw new RuntimeException("Matrix and Vector have wrong dimensions.");
        }
        for (int jj = 0; jj < arr.length; jj++) {
          data[iRow][jj] = arr[jj];
        }
    }

    public void assignRow(int iRow, int [] arr) {
        if(getColDim() != arr.length) {
            throw new RuntimeException("Matrix and Vector have wrong dimensions.");
        }
        for (int jj = 0; jj < arr.length; jj++) {
          data[iRow][jj] = arr[jj];
        }
    }

    public void assignRow(int iRow, DoubleVec dv) {
        if(getColDim() != dv.size()) {
            throw new RuntimeException("Matrix and Vector have wrong dimensions.");
        }
        for (int jj = 0; jj < dv.size(); jj++) {
          data[iRow][jj] = dv.get(jj);
        }
    }

    public boolean containsRow(DoubleVec dv) {
        boolean bOK = false;

      int iRows = getRowDim();
      int iCols = getColDim();

      for(int ii=0; ii < iRows; ii++) {
          bOK = true;
          for (int jj = 0; jj < iCols; jj++)
              if (data[ii][jj] != dv.get(jj)) {
                  bOK = false;
                  break;
              }
          if(bOK)
              break;
      }
        return bOK;
    }

    /**
     * Copies a matrix ma into this, the pointer to ma1 is not changed.
     * @param maSource
     */
    public void copy(Matrix maSource) {
		int rows = maSource.rows();
		int cols = maSource.cols();

		if (rows() == rows && cols() == cols) {
			for (int i = 0; i < rows; i++) {
				System.arraycopy(maSource.data[i], 0, data[i], 0, cols);
			}
		} else {
			data = new double[rows][];
			for (int i = 0; i < rows; i++) {

				double[] a = new double[cols];

				System.arraycopy(maSource.data[i], 0, a, 0, cols);

				data[i] = a;
			}
		}
    }

    /**
     * This (target matrix) must have appropriate dimensions.
     * @param offsetRows
     * @param maSource
     */
    public void copy(int offsetRows, Matrix maSource) {

		if(cols() != maSource.cols()){
		    throw new RuntimeException("Matrix col dimension error. Number of cols differ!");
        }

		int cols = cols();
		int r = maSource.rows();
        for (int i = 0; i < r; i++) {
            System.arraycopy(maSource.data[i], 0, data[offsetRows+i], 0, cols);
        }
    }

    public void copyColumn(Matrix maSource, int colSource, int colDestination) {
		int rows = maSource.rows();
        for (int i = 0; i < rows; i++) {
            data[i][colDestination] = maSource.get(i, colSource);
        }
    }

    final public double get(final int row, final int col) {
        return data[row][col];
    }

    public Matrix getAbs() {
        Matrix ma = new Matrix(getRowDim(), getColDim());
        for (int ii = 0; ii < getRowDim(); ii++) {
            for (int jj = 0; jj < getColDim(); jj++) {
                ma.data[ii][jj] = Math.abs(data[ii][jj]);
            }
        }
        return ma;
    }

    final public double [][] getArray() {
        return data;
    }

    public int [][] getArrayAsInt() throws NumberFormatException {
        int [][] arr = new int[data.length][data[0].length];
        for (int ii = 0; ii < getRowDim(); ii++) {
            for (int jj = 0; jj < getColDim(); jj++) {
                if(data[ii][jj] > Integer.MAX_VALUE ||
                   data[ii][jj] < Integer.MIN_VALUE) {
                    String err = "Double value instead of int " + data[ii][jj] + ".";
                    throw new NumberFormatException(err);
                }
                arr[ii][jj] = (int)data[ii][jj];
            }
        }
        return arr;
    }

    public double [][] getArrayCopy() {
        double [][] arr = new double[data.length][data[0].length];
        for (int ii = 0; ii < getRowDim(); ii++) {
            for (int jj = 0; jj < getColDim(); jj++) {
                arr[ii][jj] = data[ii][jj];
            }
        }
        return arr;
    }

    /**
     *
     * @return column wise centered matrix.
     */
    public Matrix getCenteredMatrix() {
    	
    	final int cols = cols();
    	
    	final int rows = rows();

    	double [][] arr = new double[rows][];
    	
        Matrix maMean = getMeanCols();
        
        double [] arrMean = maMean.getRow(0);
        
        for (int i = 0; i < rows; i++) {
        	
        	double [] arrRowSrc = getRow(i);
        	
        	double [] arrRow = new double [cols];
            
            System.arraycopy(arrRowSrc, 0, arrRow, 0, arrRowSrc.length);
        	
            for (int j = 0; j < cols; j++) {
            	arrRow[j] = arrRow[j] - arrMean[j];
            }
            
            arr[i]=arrRow;
            
        }
        
        Matrix ma = new Matrix();
        
        ma.setFlat(arr);
        
        return ma;
    }

    public Matrix getCenteredMatrix(Matrix maMean) {
        Matrix ma = new Matrix(getRowDim(), getColDim());
        for (int i = 0; i < rows(); i++) {
            for (int j = 0; j < cols(); j++) {
                ma.data[i][j] = data[i][j] - maMean.data[0][j];
            }
        }
        return ma;
    }

    public Matrix getCol(int col) {
        Matrix ma = new Matrix(rows(), 1);
        int iRows = getRowDim();
        for(int i=0; i < iRows; i++){
            ma.data[i][0] = data[i][col];
        }
        return ma;
    }

    public Matrix removeCol(int col2Remove) {

        int r = rows();
        int c = cols();

        Matrix ma = new Matrix(r,c-1);

        for(int i=0; i < r; i++){

            int ccColNew=0;

            for (int j = 0; j < c; j++) {
                if(j!=col2Remove) {
                    ma.data[i][ccColNew++] = data[i][j];
                }
            }
        }

        return ma;
    }

    public double [] getColAsDouble(int col) {
    	double [] arr = new double [rows()];
        for(int i=0; i < rows(); i++){
            arr[i] = data[i][col];
        }
        return arr;
    }
    
    public List<Double> getColAsList(int col) {
    	List<Double> li = new ArrayList<Double>(rows());
        for(int i=0; i < rows(); i++){
        	li.add(data[i][col]);
        }
        return li;
    }

    public float [] getColAsFloat(int iCol) {
        float [] arr = new float [rows()];
        for(int i=0; i < rows(); i++){
            arr[i] = (float)data[i][iCol];
        }
        return arr;
    }
    
    public int [] getColAsInt(int iCol) {
        int [] arr = new int [rows()];
        for(int i=0; i < rows(); i++){
            arr[i] = (int)(data[i][iCol]+0.5);
        }
        return arr;
    }

    /**
     * 
     * @param vecIndices
     * @return
     * @deprecated use getColumns(List<Integer> vecIndices) instead
     */
    public Matrix getColumns(Vector<Integer> vecIndices) {
        Matrix maReduced = new Matrix(getRowDim(), vecIndices.size());
        for (int i = 0; i < vecIndices.size(); i++) {
            int col= vecIndices.get(i);
            maReduced.assignCol(i, this, col);
        }

        return maReduced;
    }

    public Matrix getColumns(List<Integer> liIndex) {
        Matrix maReduced = new Matrix(getRowDim(), liIndex.size());
        for (int i = 0; i < liIndex.size(); i++) {
            int col = liIndex.get(i);
            maReduced.assignCol(i, this, col);
        }
        return maReduced;
    }

    public Matrix getColumns(int [] arrIndex) {
        Matrix maReduced = new Matrix(getRowDim(), arrIndex.length);
        for (int i = 0; i < arrIndex.length; i++) {
            int col = arrIndex[i];
            maReduced.assignCol(i, this, col);
        }
        return maReduced;
    }

    public int getColDim() {
        return data[0].length;
    }
    public int cols() {
        return data[0].length;
    }
    /**
     * Deletes columns with zero variance.
     * @param vecIndices the integer vector is filled with the indices from the input
     * matrix which are taken for the new matrix.
     * @return input matrix reduced by all cols with zero variance.
     */
    public Matrix getDeleteColsZeroVar(List<Integer> vecIndices) {
        vecIndices.clear();
        Matrix maReduced = new Matrix(0,0);

        Matrix maVar = getVarianceCols();
        vecIndices.clear();
        for(int ii = 0; ii < maVar.getColDim(); ii++) {
          if(maVar.get(0,ii) > 0) {
            vecIndices.add(new Integer(ii));
          }
        }
        maReduced = getColumns(vecIndices);

        return maReduced;

    }

    /** Generate a covariance matrix, each column contains values of a pulling.
     @return     An n-by-n matrix.
     */
    /*T
     <instructions>Matrix A = Matrix.random(4,4);Matrix B = A.covariance();</instructions>
     <result>A</result>
     <result>B</result>
     */

    public Matrix covariance() {

        int n = getColDim();
        int m = getRowDim();
        Matrix X = new Matrix(n, n);
        int degrees = (m - 1);
        double c;
        double s1;
        double s2;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                c = 0;
                s1 = 0;
                s2 = 0;
                for (int k = 0; k < m; k++) {
                    s1 += this.get(k, i);
                    s2 += this.get(k, j);
                }
                s1 = s1 / m;
                s2 = s2 / m;
                for (int k = 0; k < m; k++) {
                    c += (this.get(k, i) - s1) * (this.get(k, j) - s2);
                }
                X.set(i, j, c / degrees);
            }
        }
        return X;
    }

    /** Generate a correlation matrix, each column contains values of a pulling.
     @return     An n-by-n matrix.
     */
    /*T
     <instructions>Matrix A = Matrix.random(4,4);Matrix B = A.correlation();</instructions>
     <result>A</result>
     <result>B</result>
     */

    public Matrix correlation() {
        int n = getColDim();
        int m = getRowDim();
        Matrix X = new Matrix(n, n);
        int degrees = (m - 1);
        Matrix V = new Matrix(n, n);
        double c;
        double s1;
        double s2;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                c = 0;
                s1 = 0;
                s2 = 0;
                for (int k = 0; k < m; k++) {
                    s1 += this.get(k, i);
                    s2 += this.get(k, j);
                }
                s1 = s1 / m;
                s2 = s2 / m;
                for (int k = 0; k < m; k++) {
                    c += (this.get(k, i) - s1) * (this.get(k, j) - s2);
                }
                V.set(i, j, c / degrees);
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                X.set(i, j, V.get(i, j) / Math.sqrt(V.get(i, i) * V.get(j, j)));
            }
        }
        return X;
    }

    /**
     * Matrix diagonal extraction.
     * @return An d*1 Matrix of diagonal elements, d = min(m,n).
     */
    public Matrix getDiagonal() {
        int n = Math.min(getRowDim(), getColDim());
        Matrix ma = new Matrix(n,1);

        for (int i = 0; i < n; i++) {
          ma.set(i,0, data[i][i]);
        }


        return ma;
    }

    public Matrix getLinedCol() {
        Matrix ma = new Matrix(getNumElements(), 1);

        int ind = 0;
        for (int i = 0; i < getRowDim(); i++) {
            for (int j = 0; j < getColDim(); j++) {
                ma.data[ind][0] = data[i][j];
                ind++;
            }
        }

        return ma;
    }
    
    public List<Double> getList() {
    	
    	List<Double> li = new ArrayList<Double>(getRowDim()*getColDim());

        for (int i = 0; i < getRowDim(); i++) {
            for (int j = 0; j < getColDim(); j++) {
            	li.add(data[i][j]);
               
            }
        }

        return li;
    }

    public Matrix getMatrix (int i0, int i1, int j0, int j1) {
        Matrix X = new Matrix(i1-i0+1,j1-j0+1);
        double[][] B = X.getArray();
        try {
            for (int i = i0; i <= i1; i++) {
                for (int j = j0; j <= j1; j++) {
                    B[i-i0][j-j0] = data[i][j];
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }

    /** Get a submatrix.
     @param r    Array of row indices.
     @param c    Array of column indices.
     @return     A(r(:),c(:))
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */

    public Matrix getMatrix (int[] r, int[] c) {
        Matrix X = new Matrix(r.length,c.length);
        double[][] B = X.getArray();
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = 0; j < c.length; j++) {
                    B[i][j] = data[r[i]][c[j]];
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }

    /** Get a submatrix.
     @param i0   Initial row index
     @param i1   Final row index
     @param c    Array of column indices.
     @return     A(i0:i1,c(:))
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */

    public Matrix getMatrix (int i0, int i1, int[] c) {
        Matrix X = new Matrix(i1-i0+1,c.length);
        double[][] B = X.getArray();
        try {
            for (int i = i0; i <= i1; i++) {
                for (int j = 0; j < c.length; j++) {
                    B[i-i0][j] = data[i][c[j]];
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }

    /** Get a submatrix.
     @param r    Array of row indices.
     @param j0   Initial column index
     @param j1   Final column index
     @return     A(r(:),j0:j1)
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */

    public Matrix getMatrix (int[] r, int j0, int j1) {
        Matrix X = new Matrix(r.length,j1-j0+1);
        double[][] B = X.getArray();
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = j0; j <= j1; j++) {
                    B[i][j-j0] = data[r[i]][j];
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }


    public int getNumElements() {
        return cols() * rows();
    }
    
    public int getNumElementsLarger(double val2Compare) {
    	
        int ind = 0;
        for (int i = 0; i < rows(); i++) {
            for (int j = 0; j < cols(); j++) {
                if(data[i][j]>val2Compare)
                	ind++;
            }
        }
    	
        return ind;
    }
    
    public int getNumElementsEqual(double val2Compare) {
    	
        int ind = 0;
        for (int i = 0; i < rows(); i++) {
            for (int j = 0; j < cols(); j++) {
                if(data[i][j]==val2Compare)
                	ind++;
            }
        }
    	
        return ind;
    }

    /**
     * Log function from the class Math .
     * @return matrix with the natural logarithm (base e) of a double value.
     */
    public Matrix log() {
        Matrix maResult = new Matrix(getRowDim(), getColDim());
        for (int i = 0; i < rows(); i++) {
            for (int j = 0; j < cols(); j++) {
                if(data[i][j] > 0) {
                    maResult.data[i][j] = Math.log(data[i][j]);
                } else {
                    maResult.data[i][j] = 0.0;

                }
            }
        }
        return maResult;
    }

    public Matrix diagonalize() {
        Matrix maResult = new Matrix();

        if((getRowDim() > 1) && (cols() == 1)) {
            maResult.resize(rows(), rows());
            for (int i = 0; i < maResult.rows(); i++) {
                maResult.data[i][i] = data[i][0];
            }
        } else if((cols() > 1) && (rows() == 1)) {
            maResult.resize(cols() , cols());
            for (int i = 0; i < maResult.cols(); i++) {
                maResult.data[i][i] = data[0][i];
            }
        }

        return maResult;
    }

    public Matrix devide(double dDivisor) {
        Matrix maResult = new Matrix(getRowDim(), getColDim());

        for (int i = 0; i < getRowDim(); i++) {
            for (int j = 0; j < getColDim(); j++) {
                maResult.data[i][j] = data[i][j] / dDivisor;
            }
        }
        return maResult;
    }

    public Matrix devide(Matrix maDivisor) {
        Matrix maResult = new Matrix(getRowDim(), getColDim());

        for (int i = 0; i < getRowDim(); i++) {
            for (int j = 0; j < getColDim(); j++) {
                maResult.data[i][j] = data[i][j] / maDivisor.data[i][j];
            }
        }
        return maResult;
    }

    public Matrix devideDivisorBigger(Matrix maDivisor) {
        Matrix maResult = new Matrix(getRowDim(), getColDim());

        for (int i = 0; i < getRowDim(); i++) {
            for (int j = 0; j < getColDim(); j++) {
                if(data[i][j] < maDivisor.data[i][j]) {
                    maResult.data[i][j] = data[i][j] / maDivisor.data[i][j];
                } else {
                    maResult.data[i][j] = maDivisor.data[i][j] / data[i][j] ;
                }
            }
        }
        return maResult;
    }

    public Matrix devideCols(Matrix maDivisor) {
        Matrix maResult = new Matrix(getRowDim(), getColDim());

        for (int i = 0; i < getRowDim(); i++) {
            for (int j = 0; j < getColDim(); j++) {
                maResult.data[i][j] = data[i][j] / maDivisor.get(0,j);
            }
        }
        return maResult;
    }

    public void devideRow(int row, double denominator) {
        for (int j = 0; j < getColDim(); j++) {
            data[row][j] /= denominator;
        }
    }

/*
    public static void getEigenvectorBig(Matrix a, int n, Matrix d, Matrix e)
    {
        mt.DenseMatrix A = new mt.DenseMatrix(a.mData);
        // mt.Matrix B = new mt.DenseMatrix(ma.mData);
        // mt.Matrix C = new mt.DenseMatrix(getRowDim(), ma.getColDim());
        // maResult = new Matrix(mt.util.Matrices.getArray(C));

        mt.fact.EigenvalueComputer evc = new mt.fact.EigenvalueComputer(a.getRowDim());
        try {
            EigenvalueDecomposition evd = evc.factor(A);
            double [] arrEV = evd.getRealEigenvalues();
            d.resize(arrEV.length,1);
            for (int ii = 0; ii < arrEV.length; ii++) {
              d.set(ii, 0, arrEV[ii]);
            }
            a.set(mt.util.Matrices.getArray(evd.getLeftEigenvectors()));
            e.set(mt.util.Matrices.getArray(evd.getRightEigenvectors()));
        }
        catch (NotConvergedException ex) {
            ex.printStackTrace();
        }
   }
*/
    /**
     * Householder reduction according num rec 11.2.
     * Followed from an Eigen value decomposition with the calculation of the
     * Eigen vectors according num rec 11.3.
     * @param A intercept input matrix
     * @param n number of considered eigen values
     * @param D
     * @param E
     */

    final public static void getEigenvector(Matrix A, int n, Matrix D, Matrix E) {
        int l,k,j, i;
        double scale,hh,h,g,f;

       double [] a = A.toArray();
       double [] d = new double [n]; 
       double [] e = new double [n];
       
       int cols = n;
       
       // Householder reduction of a real, symmetric matrix
       for (i = n-1; i >= 1; i--) {
           l = i;
           h = scale = 0.0;
           if (l > 0) {
               for (k = 0; k < l; k++)
                   scale += Math.abs(a[i * cols + k]);
               
               if (scale == 0.0)
                   e[i] = a[i * cols + l-1];
               
               else {
            	   
                   for (k = 0; k < l; k++) {
                       a[i * cols + k] /= scale;
                       h += a[i * cols + k] * a[i * cols + k];
                   }
                   
                   f = a[i * cols + l-1];
                   g = (f >= 0.0 ? - Math.sqrt(h) : Math.sqrt(h));
                   e[i] = scale * g;
                   h -= f * g;
                   a[i * cols + l-1] = f - g;
                   f = 0.0;
                   
                   for (j = 0; j < l; j++) {
                       a[j * cols + i] = a[i * cols + j] / h;
                       g = 0.0;
                       
                       for (k = 0; k < j+1; k++)
                           g += a[j * cols + k] * a[i * cols + k];
                       
                       for (k = j+1; k < l; k++)
                           g += a[k * cols + j] * a[i * cols + k];
                       
                       e[j] = g / h;
                       f += e[j] * a[i * cols + j];
                   }
                   
                   hh = f / (h + h);
                   for (j = 0; j < l; j++) {
                       f=a[i * cols + j];
                       e[j]= g = e[j] - hh * f;
                       
                       for (k = 0; k < j+1; k++)
                           a[j * cols + k] -= (f * e[k] + g * a[i * cols + k]);
                   }
               }
           } else
               e[i] = a[i * cols + l-1];
           
           d[i] = h;
       }
       
       d[0] = 0.0;
       e[0] = 0.0;
       
       // Contents of this loop can be omitted if Eigen vectors not
       // wanted except for statement d[i]=a[i * cols + i]
       for (i = 0; i < n; i++) {
           l = i;
           if (d[i] != 0) {
               for (j = 0; j < l; j++) {
                   g = 0.0;
                   for (k = 0; k < l; k++)
                       g += a[i * cols + k] * a[k * cols + j];
                   for (k = 0; k < l; k++)
                       a[k * cols + j] -= g * a[k * cols + i];
               }
           }
           
           d[i] = a[i * cols + i];
           
           a[i * cols + i] = 1.0;
           
           for (j = 0; j < l; j++)
               a[j * cols + i] = a[i * cols + j] = 0.0;
       }

       // Eigen values and Eigen vectors of a trigiagonal matrix
       // void tqli(float d[], float e[], int n, float **z)

        int m,iter;
        double s,r,p,dd,c,b;

        for (i = 1; i < n; i++)
            e[i - 1] = e[i];
        
        e[n-1] = 0.0;
        
        for (l = 0; l < n; l++) {
        	
            iter = 0;
            
            do {
                for (m = l; m <= n-2; m++) {
                    dd = Math.abs(d[m]) + Math.abs(d[m + 1]);
                    
                    if ((double)(Math.abs(e[m]) + dd) == dd)
                        break;
                }
                
                if (m != l) {
                	
                    if (iter++ == 30)
                        System.err.println("Too many iterations in tqli");
                    
                    g=(d[l+1] - d[l]) / (2.0 * e[l]);
                    
                    r = getPythag(g,1.0);
                    
                    g=d[m] - d[l] + e[l] / (g + Sign(r,g)); // ??? Sign correct
                    
                    s=c=1.0;
                    
                    p=0.0;
                    
                    for (i = m-1; i >= l; i--) {
                    	
                        f = s * e[i];
                        
                        b = c * e[i];
                        
                        e[i+1]= (r = getPythag(f,g));
                        
                        if (r == 0.0) {
                            d[i+1] -= p;
                            e[m] = 0.0;
                            break;
                        }
                        s = f / r;
                        c = g / r;
                        
                        g = d[i+1] - p;
                        
                        r = (d[i] - g) * s + 2.0 * c * b;
                        d[i+1] = g + (p = s * r);
                        g = c * r - b;
                        for (k = 0; k < n; k++) {
                            f = a[k * cols + i + 1];
                            a[k * cols + i + 1] = s * a[k * cols + i] + c * f;
                            a[k * cols + i] = c * a[k * cols + i] - s * f;
                        }
                    }
                    if (r == 0.0 && i >= l)
                        continue;
                    
                    d[l] -= p;
                    e[l] = g;
                    e[m] = 0.0;
                }
            } while (m != l);
        }

        for (int o = 0; o < A.rows(); o++) {
			for (int q = 0; q < A.cols(); q++) {
				A.set(o,q, a[o*n+q]);
			}
		}
        
        D.resize(n, 1);
        
        D.setCol(0, d);
        
        E.resize(n, 1);

        E.setCol(0, e);
        
        return;
    }
    
//    public static void getEigenvector(MatrixK a, int n, MatrixK d, MatrixK e) {
//        int l,k,j, i;
//        double scale,hh,h,g,f;
//
//       a.resize((n + 1), (n + 1));
//       d.resize((n + 1), 1);
//       e.resize((n + 1), 1);
//
//       // Shift since routine works from indices 1..n and not from 0..n-1
//       for (i = n; i > 0; i--)
//           for (j = n; j > 0; j--)
//              a.data[i][j] = a.data[i - 1][j - 1];
//
//       // Householder reduction of a real, symmetric matrix
//        for (i = n; i >= 2; i--) {
//            l = i - 1;
//            h = scale = 0.0;
//            if (l > 1) {
//                for (k = 1; k <= l; k++)
//                    scale += Math.abs(a.data[i][k]);
//                if (scale == 0.0)
//                    e.data[i][0] = a.data[i][l];
//                else {
//                    for (k = 1; k <= l; k++) {
//                        a.data[i][k] /= scale;
//                        h += a.data[i][k] * a.data[i][k];
//                    }
//                    f = a.data[i][l];
//                    g = (f >= 0.0 ? - Math.sqrt(h) : Math.sqrt(h));
//                    e.data[i][0] = scale * g;
//                    h -= f * g;
//                    a.data[i][l] = f - g;
//                    f = 0.0;
//                    for (j = 1; j <= l; j++) {
//                        a.data[j][i] = a.data[i][j] / h;
//                        g = 0.0;
//                        for (k = 1; k <= j; k++)
//                            g += a.data[j][k] * a.data[i][k];
//                        for (k = j + 1; k <= l; k++)
//                            g += a.data[k][j] * a.data[i][k];
//                        e.data[j][0] = g / h;
//                        f += e.data[j][0] * a.data[i][j];
//                    }
//                    hh = f / (h + h);
//                    for (j = 1; j <= l; j++) {
//                        f=a.data[i][j];
//                        e.data[j][0]= g = e.data[j][0] - hh * f;
//                        for (k = 1; k <= j; k++)
//                            a.data[j][k] -= (f * e.data[k][0] + g * a.data[i][k]);
//                    }
//                }
//            } else
//                e.data[i][0] = a.data[i][l];
//            d.data[i][0] = h;
//        }
//        d.data[1][0] = 0.0;
//        e.data[1][0] = 0.0;
//        // Contents of this loop can be omitted if eigenvectors not
//        // wanted except for statement d[i]=a[i][i]
//        for (i = 1; i <= n; i++) {
//            l = i - 1;
//            if (d.data[i][0] != 0) {
//                for (j = 1; j <= l; j++) {
//                    g = 0.0;
//                    for (k = 1; k <= l; k++)
//                        g += a.data[i][k] * a.data[k][j];
//                    for (k = 1; k<= l; k++)
//                        a.data[k][j] -= g * a.data[k][i];
//                }
//            }
//            d.data[i][0] = a.data[i][i];
//            a.data[i][i] = 1.0;
//            for (j = 1; j <= l; j++)
//                a.data[j][i] = a.data[i][j] = 0.0;
//        }
//
//       // Eigenvalues and Eigenvectors of a trigiagonal matrix
//       // void tqli(float d[], float e[], int n, float **z)
//
//        int m,iter;
//        double s,r,p,dd,c,b;
//
//        for (i = 2; i <= n; i++)
//            e.data[i - 1][0] = e.data[i][0];
//        e.data[n][0] = 0.0;
//        
//        for (l = 1; l <= n; l++) {
//            iter = 0;
//            do {
//                for (m = l; m <= n-1; m++) {
//                    dd = Math.abs(d.data[m][0]) + Math.abs(d.data[m + 1][0]);
//                    
//                    if ((double)(Math.abs(e.data[m][0]) + dd) == dd)
//                        break;
//                }
//                if (m != l) {
//                    if (iter++ == 30)
//                        System.err.println("Too many iterations in tqli");
//                    g=(d.data[l+1][0] - d.data[l][0]) / (2.0 * e.data[l][0]);
//                    r = getPythag(g,1.0);
//                    g=d.data[m][0] - d.data[l][0] + e.data[l][0] / (g + Sign(r,g)); // ??? Sign correct
//                    
//                    s=c=1.0;
//                    p=0.0;
//                    for (i=m-1;i>=l;i--) {
//                        f = s * e.data[i][0];
//                        b = c * e.data[i][0];
//                        e.data[i+1][0] = (r = getPythag(f,g));
//                        if (r == 0.0) {
//                            d.data[i+1][0] -= p;
//                            e.data[m][0] = 0.0;
//                            break;
//                        }
//                        s = f / r;
//                        c = g / r;
//                        g = d.data[i+1][0] - p;
//                        r = (d.data[i][0] - g) * s + 2.0 * c * b;
//                        d.data[i+1][0] = g + (p = s * r);
//                        g = c * r - b;
//                        for (k = 1; k <= n; k++) {
//                            f = a.data[k][i + 1];
//                            a.data[k][i + 1] = s * a.data[k][i] + c * f;
//                            a.data[k][i] = c * a.data[k][i] - s * f;
//                        }
//                    }
//                    if (r == 0.0 && i >= l)
//                        continue;
//                    d.data[l][0] -= p;
//                    e.data[l][0] = g;
//                    e.data[m][0] = 0.0;
//                }
//            } while (m != l);
//        }
//
//       // Shift back;
//       for(i = 0; i < n; i++)
//           for (j = 0; j < n; j++)
//              a.data[i][j] = a.data[i + 1][j + 1];
//
//       for(j = 0; j < n; j++)
//            d.data[j][0] = d.data[j + 1][0];
//
//       for(j = 0; j < n; j++)
//           e.data[j][0] = e.data[j + 1][0];
//
//       // Resize arrays properly
//       a.resize(n, n);
//       d.resize(n, 1);
//       e.resize(n, 1);
//       
//       return;
//    }
    /**
     *
     * @param row row
     * @return smallest value in specified row.
     */
    public double getMinRow(int row) {
        double val = Double.MAX_VALUE;
        for (int i = 0; i < getColDim(); i++) {
          if(data[row][i] < val) {
              val = data[row][i];
          }
        }
        return val;
    }

    /**
     *
     * @param row row
     * @return the column index of the smallest value in the row.
     */
    public int getMinRowIndex(int row) {
        int index = 0;
        double min = Double.MAX_VALUE;
        for (int i = 0; i < getColDim(); i++) {
          if(data[row][i] < min) {
              min = data[row][i];
              index = i;
          }
        }
        return index;
    }

    /**
     * 03.10.04 MvK
     * Shuffles the indices before searching the minimum. This randomizes the
     * search order.
     * @param row row
     * @return the column index of the smallest value in the row.
     */
    public int getMinRowIndexRND(int row) {
        int indexMin = 0;
        List<Integer> list = new ArrayList<Integer>();
        for (int i = 0; i < getColDim(); i++) {
          list.add(new Integer(i));
        }
        Collections.shuffle(list);

        int [] arrIndex = new int [getColDim()];
        for (int i = 0; i < arrIndex.length; i++) {
          arrIndex[i] = ((Integer) list.get(i)).intValue();
        }

        double min = Double.MAX_VALUE;
        for (int i = 0; i < getColDim(); i++) {
            int ind = arrIndex[i];
            if (data[row][ind] < min) {
                min = data[row][ind];
                indexMin = ind;
            }
        }
        return indexMin;
    }

    public int getID() {
        return identifier;
    }

    /**
     * Skipping the sqrt.
     * @param iRow
     * @param A
     * @param iRowA
     * @return
     */
    public double getEuclideanDistanceFastRows(int iRow, Matrix A, int iRowA) {
        double distance = 0;

        int cols = cols();

        for (int i = 0; i < cols; i++) {
          distance += (data[iRow][i] - A.data[iRowA][i]) * (data[iRow][i] - A.data[iRowA][i]);
        }
        return distance;
    }

    /**
     * Skipping the sqrt.
     * @param row1
     * @param row2
     * @return
     */
    public double getEuclideanDistanceFastRows(int row1, int row2) {
        double distance = 0;

        int cols = cols();

        for (int i = 0; i < cols; i++) {
          distance += (data[row1][i] - data[row2][i]) * (data[row1][i] - data[row2][i]);
        }

        return distance;
    }

    public double getMaximumValue() {
        double max = Double.MIN_VALUE;

        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length; j++) {
                if(max < data[i][j])
                    max = data[i][j];
            }
        }
        return max;
    }

    public double getMean() {
        int cols = getColDim();
        int rows = getRowDim();
        double mean = getSum() / (cols * rows);
        return mean;
    }

    public Matrix getMeanCols() {
        
    	final int cols = cols();
    	
    	final int rows = rows();
    	
    	final double [] a = new double [cols];

    	for (int i = 0; i < rows; i++) {
    		final double [] arrRow = data[i];
    		
    		for (int j = 0; j < arrRow.length; j++) {
    			a[j] += arrRow[j];
			}
		}

    	for (int i = 0; i < cols; i++) {
    		a[i] /= rows;
    	}

        return new Matrix(true, a);
    }

    public Matrix getMeanCols(int [] rowIndex) {
        Matrix ma = new Matrix(1, getColDim());
        int cols = getColDim();
        for (int j = 0; j < cols; j++) {
            ma.data[0][j] = getMeanCol(j, rowIndex);
        }
        return ma;
    }

    public double getMeanCol(int col) {

        int rows = rows();
        
        double sum = 0;
        
        for (int i = 0; i < rows; i++) {
            sum += data[i][col];
        }

        return sum / rows;
    }

    public double getMeanRow(int row) {

        int iCols = getColDim();
        double dMean = 0;
        for (int i = 0; i < iCols; i++) {
            dMean += data[row][i];
        }

        return dMean / iCols;
    }
    
    public Matrix getMeanRows() {
        Matrix ma = new Matrix(1, getRowDim());
        int iRows = getRowDim();
        for (int i = 0; i < iRows; i++) {
        	ma.set(0,i,getMeanRow(i));
        }
        return ma;
    }

    public double getMeanCol(int iColIndex, int [] rowIndex) {

        int iRows = getRowDim();
        double dMean = 0;
        for (int ii = 0; ii < rowIndex.length; ii++) {
            dMean += data[rowIndex[ii]][iColIndex];
        }

        return dMean / iRows;
    }
    
    public double getMedian(int row) {
    	
    	double [] arr = getRowCopy(row);
    	
    	Arrays.sort(arr);
    	
		double median = 0;
		int len = arr.length;
		if(len % 2 != 0) {
			median = arr[len / 2];
		} else {
			int ind = (int)(((double)len / 2.0)+0.5);
			median = (arr[ind] + arr[ind-1]) / 2.0;
		}
		
		return median;
    }

    public double getMedian() {
    	
    	int nRows = rows();
    	int nCols = cols();
    	
    	double [] arr = new double [nRows*nCols]; 
    	
    	int k = 0;
    	for (int i = 0; i < nRows; i++) {
			for (int j = 0; j < nCols; j++) {
				arr[k++]=get(i, j);
			}
		}
    	
    	Arrays.sort(arr);
    	
		double median = 0;
		int len = arr.length;
		if(len % 2 != 0) {
			median = arr[len / 2];
		} else {
			int ind = (int)(((double)len / 2.0)+0.5);
			median = (arr[ind] + arr[ind-1]) / 2.0;
		}
		
		return median;
    }

    public Matrix getMedianCols() {

        Matrix maMedian = new Matrix(1, cols());

    	int rows = rows();
    	int cols = cols();

    	double [] arr = new double [rows];

        for (int i = 0; i < cols; i++) {

            for (int j = 0; j < rows; j++) {
                arr[j]=get(j,i);
            }

            Arrays.sort(arr);

            double median = 0;
            int len = arr.length;
            if(len % 2 != 0) {
                median = arr[len / 2];
            } else {
                int ind = (int)(((double)len / 2.0)+0.5);
                median = (arr[ind] + arr[ind-1]) / 2.0;
            }

            maMedian.set(0, i, median);

        }

		return maMedian;
    }

    public Matrix getMergeRows(Matrix ma) {

        Matrix maMerge = new Matrix(getRowDim() + ma.getRowDim(), getColDim());

        for (int i = 0; i < getRowDim(); i++) {
            for (int j = 0; j < getColDim(); j++) {
                maMerge.data[i][j] = data[i][j];
            }
        }

        int iRow = getRowDim();
        for (int ii = 0; ii < ma.getRowDim(); ii++) {
            for (int jj = 0; jj < getColDim(); jj++) {
                maMerge.data[iRow + ii][jj] = ma.data[ii][jj];
            }
        }

        return maMerge;
    }

    /**
     * max and min vals for all cols  
     * @return matrix with two rows.
     */
    public Matrix getMaxMin() {

        Matrix maMaxMin = new Matrix(2, getColDim());

        int cols = getColDim();
        
        for (int i = 0; i < cols; i++) {
            maMaxMin.set(1,i, getMin(i));
            maMaxMin.set(0,i, getMax(i));
        }
        return maMaxMin;
    }

    /**
     *
     * @return
     */
    public double getMax() {

        int rows = getRowDim();
        int cols = getColDim();
        double dMax = -Double.MAX_VALUE;
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                if(data[i][j] > dMax)
                    dMax = data[i][j];

        return dMax;
    }

    /**
     * get max value for that col.
     * @param col
     * @return
     */
    public double getMax(int col) {

        int iRows = getRowDim();
        
        double max = -Double.MAX_VALUE;
        
        for (int i = 0; i < iRows; i++) {
            if(data[i][col] > max)
                max = data[i][col];
        }
        
        return max;
    }
    
    /**
     * 
     * @return Point(row,col) for the largest value.
     */
    public Point getMaxIndex() {

    	Point p = new Point(0,0);
    	
        int rows = rows();
        
        int cols = cols();
        
        double max = data[0][0];
        
        for (int i = 0; i < rows; i++) {
        	for (int j = 0; j < cols; j++) {
                if(data[i][j] > max) {
                    max = data[i][j];
                    p.y = i;
                    p.x = j;
                }
			}
        }
        
        return p;
    }

    public double getMaxRow(int row) {
        double max = -Double.MAX_VALUE;
        for (int i = 0; i < getColDim(); i++) {
            if(data[row][i] > max)
                max = data[row][i];
        }
        return max;
    }

    /**
     * Shuffles the indices before searching the minimum.
     * @param row row
     * @return the column index of the largest value in the row.
     */
    public int getMaxRowIndexRND(int row) {
        int indexMax = 0;
        
        List<Integer> list = new ArrayList<Integer>();
        
        for (int i = 0; i < cols(); i++) {
          list.add(i);
        }
        Collections.shuffle(list);
        int [] arrIndex = new int [getColDim()];
        for (int i = 0; i < arrIndex.length; i++) {
          arrIndex[i] = ((Integer) list.get(i)).intValue();
        }
        double max = -Double.MAX_VALUE;
        for (int i = 0; i < getColDim(); i++) {
            int ind = arrIndex[i];
            if (data[row][ind] > max) {
                max = data[row][ind];
                indexMax = ind;
            }
        }
        return indexMax;
    }
    
    /**
     * @param col
     * @return the row index of the largest value in the col.
     */
    public int getMaxRowIndex(int col) {
        int row = -1;
        
        double max = -Double.MAX_VALUE;
        
        for (int i = 0; i < rows(); i++) {
            
            if (data[i][col] > max) {
                max = data[i][col];
                row = i;
            }
        }
        
        return row;
    }

    public double getMin() {

        int rows = getRowDim();
        int cols = getColDim();
        double min = Double.MAX_VALUE;
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                if(data[i][j] < min)
                    min = data[i][j];

        return min;
    }
/**
 * 
 * @return ValPos
 */
    public ScorePoint getMinPos() {

        int rows = getRowDim();
        int cols = getColDim();
        
        ScorePoint td = new ScorePoint();
        
        td.setScore(Double.MAX_VALUE);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                if(data[i][j] < td.getScore()) {
                	td.y = i;
                	td.x = j;
                	td.setScore(data[i][j]);
                }

        return td;
    }
/**
 * 
 * @return matrix with minimum value from each row.
 */
    public Matrix getMinRows() {
        int rows = getRowDim();
        int cols = getColDim();
        Matrix ma = new Matrix(rows, 1);
        for (int i = 0; i < rows; i++) {
        	double dMin = Double.MAX_VALUE;
            for (int j = 0; j < cols; j++)
                if(data[i][j] < dMin)
                    dMin = data[i][j];
            ma.set(i, 0, dMin);
        }

        return ma;
    }
/**
 * 
 * @return matrix dimension n,2. First row contains column position of min val,
 * second col the min value.
 */
    public Matrix getMinRowsPosCol() {
        int rows = getRowDim();
        int cols = getColDim();
        Matrix ma = new Matrix(rows, 2);
        for (int i = 0; i < rows; i++) {
        	double dMin = Double.MAX_VALUE;
            for (int j = 0; j < cols; j++)
                if(data[i][j] < dMin) {
                	ma.set(i, 0, j);
                	ma.set(i, 1, data[i][j]);
                }
        }

        return ma;
    }

    public double getMin(int col) {

        int rows = getRowDim();
        double min = Double.MAX_VALUE;
        for (int i = 0; i < rows; i++) {
            if(data[i][col] < min)
                min = data[i][col];
        }
        return min;
    }
	/**
	 * 
	 * @param row
	 * @return column index for the field with the lowest value for the given row.
	 */
    public int getMinColIndex(int row) {

        int cols = getColDim();
        double min = Double.MAX_VALUE;
        int ind = -1;
        for (int i = 0; i < cols; i++) {
            if(data[row][i] < min) {
                min = data[row][i];
                ind=i;
            }
        }
        return ind;
    }
    
	/**
	 * 
	 * @param row
	 * @return column index for the field with the largest value for the given row.
	 */
    public int getColIndexContainingMaxVal(int row) {

        int cols = getColDim();
        double max = -Double.MAX_VALUE;
        int col = -1;
        for (int i = 0; i < cols; i++) {
            if(data[row][i] > max) {
                max = data[row][i];
                col=i;
            }
        }
        return col;
    }
    
    /**
     * 
     * @param rowEnd the values in the matrix will be considered until this row (exclusively).
     * @param colEnd the values in the matrix will be considered until this col (exclusively).
     * @return
     */
    public ScorePoint getColIndexContainingMaxVal(int rowEnd, int colEnd) {
    	
    	int cols = colEnd;
    	
    	int rows = rowEnd;
    	        
        double max = -Double.MAX_VALUE;
        
        int rowMax = -1;
        int colMax = -1;
        
        for (int i = 0; i < cols; i++) {
        	
        	double maxInCol = -Double.MAX_VALUE;
        	
        	int rowMaxInCol=-1;
        	for (int j = 0; j < rows; j++) {
        		if(data[j][i] > maxInCol) {
        			maxInCol = data[j][i];
        			rowMaxInCol=j;
                }
			}
        	
        	if(maxInCol>max){
        		max = maxInCol;
        		rowMax = rowMaxInCol;
        		colMax = i;
        	}
        }
        
        ScorePoint sc = new ScorePoint(rowMax, colMax);
        
        sc.setScore(max);
        
        return sc;
    }

    public Matrix getNormalizedMatrix() {

        int iRows = getRowDim();
        int iCols = getColDim();

        Matrix maNorm = new Matrix(iRows, iCols);
        Matrix maStandardDeviation = getStandardDeviationCols();
        for (int ii = 0; ii < iRows; ii++)
            for (int jj = 0; jj < iCols; jj++) {
                double dStandardDeviation = maStandardDeviation.get(0, jj);
                if (dStandardDeviation == 0) {
                    dStandardDeviation = Float.MIN_VALUE;
                }
                double dVal = get(ii, jj) / dStandardDeviation;
                maNorm.set(ii, jj, dVal);
            }
        return maNorm;
    }

    public Matrix getNormalizedMatrix(Matrix maStandardDeviation) {

        int iRows = getRowDim();
        int iCols = getColDim();

        Matrix maNorm = new Matrix(iRows, iCols);

        for (int ii = 0; ii < iRows; ii++)
            for (int jj = 0; jj < iCols; jj++) {
                double dStandardDeviation = maStandardDeviation.get(0, jj);
                if (dStandardDeviation == 0) {
                    dStandardDeviation = Float.MIN_VALUE;
                }
                double dVal = get(ii, jj) / dStandardDeviation;
                maNorm.set(ii, jj, dVal);
            }
        return maNorm;
    }

    /**
     * For a quadratic matrix only.
     * @return
     */
    public double [] getUpperTriangle(){

        int r = rows();
        int c = cols();

        if(r != c){
            throw new RuntimeException("Not a quadratic matrix.");
        }

        int n = ((r * r) - r)/2;

        double [] a = new double[n];

        int cc = 0;
        for (int i = 0; i < r; i++) {

            for (int j = i+1; j < r; j++) {
                a[cc++] = get(i,j);
            }
        }

        return a;
    }

    /**
     *pythag computes sqrt(a^2 + b^2) without destructive underflow or overflow.
     * @param a length a
     * @param b length b
     * @return double
     */
    public static double getPythag(double a, double b)
    {
      double absa = Math.abs(a);
      double absb = Math.abs(b);
      if (absa > absb) {
          double dScale = (absb / absa) * (absb / absa);
          return absa * Math.sqrt(1.0 + dScale);
      }
      else if(absb == 0.0)
          return 0.0;
      else {
          double dScale = (absa / absb) * (absa / absb);
          return (absb * Math.sqrt(1.0 + dScale));
      }
    }

    public static double getPythag2(double a, double b)
    {
        return Math.sqrt((a * a) + (b * b));
    }


    public boolean areRowsEqual(int row1, int row2) {

        boolean equal = true;

        int cols = cols();

        for (int i = 0; i < cols; i++) {

            double diff = Math.abs(data[row1][i]-data[row2][i]);

            if(diff>TINY){
                equal=false;
                break;
            }
        }

        return equal;
    }

    public boolean equal(Matrix ma) {
        boolean bEQ = true;

        if(equalDimension(ma)) {
            for (int i = 0; i < getRowDim(); i++) {
                for (int j = 0; j < getColDim(); j++) {
                    if(data[i][j] != ma.data[i][j]) {
                        bEQ = false;
                        break;
                    }
                }
            }
        } else {
            bEQ = false;
        }

        return bEQ;
    }

    public boolean equal(Matrix ma, double dLimit) {
        boolean bEQ = true;

        if(equalDimension(ma)) {
            for (int i = 0; i < getRowDim(); i++) {
                for (int j = 0; j < getColDim(); j++) {
                    double dDiff = Math.abs(data[i][j] - ma.data[i][j]);
                    if(dDiff > dLimit) {
                        bEQ = false;
                        break;
                    }
                }
            }
        } else {
            bEQ = false;
        }

        return bEQ;
    }

    /**
     * Checks two matrices for equal dimensions.
     * @param ma matrix to compare with.
     * @return true if the number of getColDim and the number of getRowDim are
     * corresponding.
     */
    public boolean equalDimension(Matrix ma) {
        boolean bEqual = false;

        if((getColDim() == ma.getColDim()) && (getRowDim() == ma.getRowDim()))
            bEqual = true;

        return bEqual;
    }

    public boolean hasOnlyFinite() {

        boolean finite = true;

        int r = rows();

        int c = cols();

        ma:
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                if(!Double.isFinite(data[i][j])){
                    finite = false;
                    break ma;
                }
            }
        }

        return finite;

    }
    /**
     * Value by value multiplication
     * Matrices must have the same dimensions.
     * @param ma
     * @return
     */
    public Matrix multiplyValByVal(Matrix ma) {
    	
    	Matrix maProd = new Matrix(rows(), cols());
    	
    	for (int i = 0; i < rows(); i++) {
			for (int j = 0; j < cols(); j++) {
				maProd.set(i, j, get(i,j) * ma.get(i, j));
			}
		}
    	
    	return maProd;
    	
    }

    
    public Matrix multCols(Matrix maMuiltiplicant) {
        Matrix maResult = new Matrix(getRowDim(), getColDim());

        for (int i = 0; i < getRowDim(); i++) {
            for (int j = 0; j < getColDim(); j++) {
                maResult.data[i][j] = data[i][j] * maMuiltiplicant.get(0,j);
            }
        }
        return maResult;
    }

    public Matrix multiply(double dScalar) {
        Matrix maResult = new Matrix(getRowDim(), getColDim());

        for (int i = 0; i < getRowDim(); i++) {
            for (int j = 0; j < getColDim(); j++) {
                maResult.data[i][j] = data[i][j] * dScalar;
            }
        }

        return maResult;
    }

    public Matrix multiply(Matrix ma) {
        return multiply(false,false,ma);
    }

    public double multiply(int row1, Matrix ma2, int row2) {
        double sum = 0;
        for (int ii = 0; ii < getColDim(); ii++) {
            sum += data[row1][ii] * ma2.data[row2][ii];
        }
        return sum;
    }


    /**
     * Multiplication of two matrices. The algorithm checks for the correct
     * dimensions of the matrices. If the dimensions are not corresponding a
     * runtime exception is thrown. The left matrix is the instantiated object.
     * @param transA if true the left matrix is transposed before the
     * multiplication
     * @param transB if true the right matrix is transposed before the
     * multiplication
     * @param ma right matrix
     * @return matrix with dimension getRowDim(), ma.getColDim().
     */
//    public MatrixK multiply(boolean transA, boolean transB, MatrixK ma) {
//        MatrixK maResult = null;
//
//
//        // Tasks:
//        // Computes
//        // C = A * B
//        // C = A' * B
//        // C = A * B'
//        // C = A' * B'
//        //
//        
//        // C = A * B
//        if ( (transA == false) && (transB == false)) {
//        	MatrixK B = ma.getTranspose();
//            
//        	maResult = multATransB(B);
//
//        }
//
//        //
//        // C = A' * B
//        //
//        // Matrix A is transposed
//        if ( (transA == true) && (transB == false)) {
//
//        	MatrixK A = getTranspose();
//        	MatrixK B = ma.getTranspose();
//            
//        	maResult = A.multATransB(B);
//            
//        } 
//
//        //
//        // C = A * B'
//        //
//        // Matrix B is transposed
//        if ( (transA == false) && (transB == true)) {
//        	maResult = multATransB(ma);
//        } 
//
//        //
//        // C = A' * B'
//        //
//        // Matrix A + B are transposed
//        if ( (transA == true) && (transB == true)) {
//        	MatrixK A = getTranspose();
//        	maResult = A.multATransB(ma);
//        }
//
//        return maResult;
//    }
    
    
    
//    private MatrixK multATransB(MatrixK ma){
//        //
//        // C = A * B'
//        //
//        int n = cols();
//        
//        int m = rows();
//        
//        int maRows = ma.rows();
//        
//        if (n != ma.cols()) {
//            throw new RuntimeException("Error in Routine SMatrix::Mult(). Attempt to calculate the product of two incompatible matrices. Do nothing and return.");
//        }
//
//        MatrixK maResult = new MatrixK(rows(), maRows);
//        
//        double[][] C = maResult.getArray();
//        
//        double[] Bcolj = new double[n];
//        
//        for (int j = 0; j < maRows; j++) {
//           for (int k = 0; k < n; k++) {
//              Bcolj[k] = ma.data[j][k];
//           }
//           for (int i = 0; i < m; i++) {
//              double[] Arowi = data[i];
//              double s = 0;
//              for (int k = 0; k < n; k++) {
//                 s += Arowi[k]*Bcolj[k];
//              }
//              C[i][j] = s;
//           }
//        }
//        return maResult;
//    }
    
    
    /**
     * Multiplication of two matrices. The algorithm checks for the correct
     * dimensions of the matrices. If the dimensions are not corresponding a
     * runtime exception is thrown. The left matrix is the instantiated object.
     * Algorithm taken from Sedgewick.
     * @param transA if true the left matrix is transposed before the
     * multiplication
     * @param transB if true the right matrix is transposed before the
     * multiplication
     * @param ma right matrix
     * @return matrix with dimension getRowDim(), ma.getColDim().
     */
    
    final public Matrix multiply(boolean transA, boolean transB, Matrix ma) {
        Matrix maResult = null;


        // Tasks:
        // Computes
        // C = A * B
        // C = A' * B
        // C = A * B'
        // C = A' * B'
        //
        
        // C = A * B
        if (!transA && !transB) {
        	
            int n = cols();
            
            int m = rows();
            
            int maCols = ma.cols();
            
            if (n != ma.rows()) {
                throw new RuntimeException("Error in Routine Matrix.multiply(...). Attempt to calculate the product of two incompatible matrices. Do nothing and return.");
            }

            maResult = new Matrix(rows(), ma.cols());
            
            double[][] C = maResult.getArray();
            
            Matrix Btrans = ma.getTranspose(); 
            
            double[][] bTrans = Btrans.getArray(); 
            
            double [] Bcolj = null;
            
            double [] Arowi = null;
            
            int c4 = n/4 * 4;
            
            for (int j = 0; j < maCols; j++) {
               
               Bcolj = bTrans[j];
               
               for (int i = 0; i < m; i++) {
                  Arowi = data[i];
                  
                  double s = 0;
                  
                  for (int k = 0; k < c4; k+=4) {
                     
                	  final double s1 = Arowi[k]*Bcolj[k];
                	  final double s2 = Arowi[k+1]*Bcolj[k+1];
                	  final double s3 = Arowi[k+2]*Bcolj[k+2];
                	  final double s4 = Arowi[k+3]*Bcolj[k+3];
                	  
                	  s += s1+s2+s3+s4;
                  }
                  
                  for (int k = c4; k < n; k++) {
                      s += Arowi[k]*Bcolj[k];
                  }
                                    
                  C[i][j] = s;
               }
            }

        }

        //
        // C = A' * B
        //
        // Matrix A is transposed
        if (transA && !transB) {

            // Reverse n
            int n = rows();
            
            int m = cols();
            
            int maCols = ma.cols();
            
            if (n != ma.rows()) {
                throw new RuntimeException("Error in Routine SMatrix::Mult(). Attempt to calculate the product of two incompatible matrices. Do nothing and return.");
            }
            
            // Perform C = A' * B
            maResult = new Matrix(cols(), ma.cols());
            
            double[][] C = maResult.getArray();
            
            double [] Bcolj = null;
            
            int c4 = n/4 * 4;
            
            Matrix Atrans = getTranspose(); 
            
            double [][] aTrans = Atrans.getArray();
            
            Matrix Btrans = ma.getTranspose(); 
            
            double[][] bTrans = Btrans.getArray();
            
            double [] Arowi = null;
            
            for (int j = 0; j < maCols; j++) {
               
               Bcolj = bTrans[j];
               
               for (int i = 0; i < m; i++) {
            	   
            	  Arowi = aTrans[i];
                  double s = 0;
                  
                  for (int k = 0; k < c4; k+=4) {
                      
                	  final double s1 = Arowi[k]*Bcolj[k];
                	  final double s2 = Arowi[k+1]*Bcolj[k+1];
                	  final double s3 = Arowi[k+2]*Bcolj[k+2];
                	  final double s4 = Arowi[k+3]*Bcolj[k+3];
                 	  
                 	  s += s1+s2+s3+s4;
                   }
                   
                   for (int k = c4; k < n; k++) {
                       s += Arowi[k]*Bcolj[k];
                   }
                   
                  C[i][j] = s;
               }
            }
            
        } 

        //
        // C = A * B'
        //
        // Matrix B is transposed
        if (!transA && transB) {
            int n = cols();
            
            int m = rows();
            
            int maRows = ma.rows();
            
            if (n != ma.cols()) {
                throw new RuntimeException("Error in Routine SMatrix::Mult(). Attempt to calculate the product of two incompatible matrices. Do nothing and return.");
            }

            maResult = new Matrix(rows(), maRows);
            
            double[][] C = maResult.getArray();
            
            double [] Bcolj = null;
            
            int c4 = n/4 * 4;
            
            double [] Arowi = null;
            
            for (int j = 0; j < maRows; j++) {
               
               Bcolj = ma.data[j];
               
               for (int i = 0; i < m; i++) {
            	   
                  Arowi = data[i];
                  
                  double s = 0;
                  
                  for (int k = 0; k < c4; k+=4) {
                      
                	  final double s1 = Arowi[k]*Bcolj[k];
                	  final double s2 = Arowi[k+1]*Bcolj[k+1];
                	  final double s3 = Arowi[k+2]*Bcolj[k+2];
                	  final double s4 = Arowi[k+3]*Bcolj[k+3];
                 	  
                 	  s += s1+s2+s3+s4;
                   }
                   
                   for (int k = c4; k < n; k++) {
                       s += Arowi[k]*Bcolj[k];
                   }
                  
                  C[i][j] = s;
               }
            }
        } 

        //
        // C = A' * B'
        //
        // Matrix A + B are transposed
        if (transA && transB) {
        	
            // Reverse n
            int n = rows();
            
            int m = cols();
            
            int maRows = ma.rows();
            
            if (n != ma.cols()) {
                throw new RuntimeException("Error in Routine SMatrix::Mult(). Attempt to calculate the product of two incompatible matrices. Do nothing and return.");
            }
            
            // Perform C = A' * B'
            maResult = new Matrix(cols(), maRows);
            
            double[][] C = maResult.getArray();
            
            double [] Bcolj = null;
            
            int c4 = n/4 * 4;
            
            Matrix Atrans = getTranspose(); 
            
            double [][] aTrans = Atrans.getArray();
            
            double [] Arowi = null;
            
            for (int j = 0; j < maRows; j++) {
               
               Bcolj = ma.data[j];
               
               for (int i = 0; i < m; i++) {
            	   
             	  Arowi = aTrans[i];
             	  
                  double s = 0;
                  
                  for (int k = 0; k < c4; k+=4) {
                      
                	  final double s1 = Arowi[k]*Bcolj[k];
                	  final double s2 = Arowi[k+1]*Bcolj[k+1];
                	  final double s3 = Arowi[k+2]*Bcolj[k+2];
                	  final double s4 = Arowi[k+3]*Bcolj[k+3];
                 	  
                 	  s += s1+s2+s3+s4;
                   }
                   
                   for (int k = c4; k < n; k++) {
                       s += Arowi[k]*Bcolj[k];
                   }
                   
                  C[i][j] = s;
               }
            }
        }

        return maResult;
    }
    
/*
    public Matrix multiplyBig(boolean bTransA, boolean bTransB, Matrix ma) {
        Matrix maResult = new Matrix(getRowDim(), ma.getColDim());


        // Tasks:
        // Computes
        // C = A * B
        // C = A' * B
        // C = A * B'
        // C = A' * B'
        //
        // ==>> where *this.Val == C
        //
        // Input:
        // A() n x m matrix
        // B() s x r matrix, where m == s
        //
        // Output:
        // C() n x r matrix (*this.Val is changed)
        //
        // Requirements:
        // m = s
        //
        // Resulting dimensions of matrix product:
        // Matrix C => n x r (C == *this.Val)
        //
        // Remark:
        // For Var/Covar-Matrix the routine requires a centered Matrix
        // => more options could be added here

        int i, j, k, n, m, s, r;
        //
        // C = A * B
        //
        // No matrix should be transposed
        if ( (bTransA == false) && (bTransB == false)) {
            mt.Matrix A = new mt.DenseMatrix(mData);
            mt.Matrix B = new mt.DenseMatrix(ma.mData);
            mt.Matrix C = new mt.DenseMatrix(getRowDim(), ma.getColDim());
            C = A.mult(B,C);
            maResult = new Matrix(mt.util.Matrices.getArray(C));

        } // End if

        //
        // C = A' * B
        //
        // Matrix A should be transposed
        if ( (bTransA == true) && (bTransB == false)) {

            mt.Matrix A = new mt.DenseMatrix(mData);
            mt.Matrix B = new mt.DenseMatrix(ma.mData);
            mt.Matrix C = new mt.DenseMatrix(getColDim(), ma.getColDim());
            C = A.transAmult(B,C);
            maResult = new Matrix(mt.util.Matrices.getArray(C));

        } // End if

        //
        // C = A * B'
        //
        // Matrix B should be transposed
        if ( (bTransA == false) && (bTransB == true)) {
            mt.Matrix A = new mt.DenseMatrix(mData);
            mt.Matrix B = new mt.DenseMatrix(ma.mData);
            mt.Matrix C = new mt.DenseMatrix(getRowDim(), ma.getRowDim());
            C = A.transBmult(B,C);
            maResult = new Matrix(mt.util.Matrices.getArray(C));
        } // End if

        //
        // C = A' * B'
        //
        // Matrix A + B should be transposed
        if ( (bTransA == true) && (bTransB == true)) {
            mt.Matrix A = new mt.DenseMatrix(mData);
            mt.Matrix B = new mt.DenseMatrix(ma.mData);
            mt.Matrix C = new mt.DenseMatrix(getColDim(), ma.getRowDim());
            C = A.transABmult(B,C);
            maResult = new Matrix(mt.util.Matrices.getArray(C));
        }

        return maResult;
    }
*/

/**
 * Parses an vector with Strings and converts it to a Matrix.
 * @param vecStringMatrix vector with strings.
 * @return Matrix
 */
    public static Matrix getParsed(Vector<String> vecStringMatrix) throws Exception {

        Matrix ma = new Matrix();
        String sRow = vecStringMatrix.get(0);

        int iLen = 0;

        double[] dArr = String2DoubleArray.convert(sRow);
        iLen = dArr.length;

        ma.resize(vecStringMatrix.size(), dArr.length);

        for (int ii = 0; ii < dArr.length; ii++) {
            ma.set(0, ii, dArr[ii]);
        }

        for (int ii = 1; ii < vecStringMatrix.size(); ii++) {
            sRow = (String) vecStringMatrix.get(ii);
            dArr = String2DoubleArray.convert(sRow);

            if (iLen != dArr.length) {
                System.err.println(
                    "Vectors for matrix generation differ in length");
                (new RuntimeException()).printStackTrace();
                break;
            }

            for (int jj = 0; jj < dArr.length; jj++) {
                ma.set(ii, jj, dArr[jj]);
            }
        }
        return ma;
    }


    public Matrix plus(Matrix ma) {
        Matrix maResult = new Matrix(getRowDim(), getColDim());
        if(!equalDimension(ma)) {
            throw new RuntimeException("Matrices have wrong dimensions.");
        }

        for (int ii = 0; ii < getRowDim(); ii++) {
            for (int jj = 0; jj < getColDim(); jj++) {
                maResult.data[ii][jj] = data[ii][jj] + ma.data[ii][jj];
            }
        }
        return maResult;
    }

    public Matrix pow(double exp) {
        Matrix ma = new Matrix(getRowDim(), getColDim());
        for (int i = 0; i < getRowDim(); i++) {
            for (int j = 0; j < getColDim(); j++) {
                ma.data[i][j] = Math.pow(data[i][j], exp);
            }
        }
        return ma;
    }

    /** Matrix determinant
     @return     determinant
     */

    public double det () {
        return new LUDecomposition(this).det();
    }

    /** Matrix rank
     @return     effective numerical rank, obtained from SVD.
     */



    /** Matrix condition (2 norm)
     @return     ratio of largest to smallest singular value.
     */



    /** Matrix trace.
     @return     sum of the diagonal elements.
     */


    /**
     * Resizes the matrix, the old data are written to the new matrix, if
     * possible.
     * @param rowsNew new number of getRowDim.
     * @param colsNew new number of getColDim.
     */
    public void resize(int rowsNew, int colsNew) {
        double [][] arrTmp = new double [rowsNew][];

        int rowsMin = Math.min(rows(), rowsNew);
        
        int colsMin = Math.min(cols(), colsNew);

        for (int i = 0; i < rowsMin; i++) {
			double [] a = new double [colsNew];
			
			System.arraycopy(data[i], 0, a, 0, colsMin);
			
			arrTmp[i]=a;
			
		}
        
        for (int i = rowsMin; i < rowsNew; i++) {
			double [] a = new double [colsNew];
			
			arrTmp[i]=a;
			
		}
        
        data = arrTmp;
    }
    
    
    
    
//    public void resize(int iNumberRowsNew, int iNumberColsNew) {
//        double [][] dTmp = new double [iNumberRowsNew][iNumberColsNew];
//
//        int iRows = iNumberRowsNew;
//        int iCols = iNumberColsNew;
//
//        if(iNumberRowsNew > getRowDim())
//            iRows = getRowDim();
//        if(iNumberColsNew > getColDim())
//            iCols = getColDim();
//        for (int ii = 0; ii < iRows; ii++) {
//            for (int jj = 0; jj < iCols; jj++) {
//                dTmp[ii][jj] = data[ii][jj];
//            }
//        }
//        data = dTmp;
//    }

    /**
     * Flat copy.
     * @param row
     * @return
     */
    public double [] getRow(int row) {
        return data[row];
    }

    public double [] getRowCopy(int row) {
    	double [] arr = new double[data[0].length];
    	for (int i = 0; i < arr.length; i++) {
			arr[i] = data[row][i];
		}
        return arr;
    }

    public List<Double> getRowAsList(int row) {
        List<Double> list = new ArrayList<Double>();
        for (int ii = 0; ii < data[0].length; ii++) {
          list.add(data[row][ii]);
        }
        return list;
    }
    
    public float [] getRowAsFloat(int row) {
    	float [] arr = new float [data[0].length];
        for (int i = 0; i < data[0].length; i++) {
        	arr[i]=(float)data[row][i];
        }
        return arr;
    }

    public int getRowDim() {
        return data.length;
    }
    public int rows() {
        return data.length;
    }
    /**
     * Sorts the row of the matrix according the compareTo(...) function in
     * DoubleVec
     * @return sorted matrix
     */
    public Matrix getSorted() {
    	
        Matrix ma = new Matrix(getRowDim(), getColDim());
        
        List<DoubleVec> list = new ArrayList<DoubleVec>();
        
        for (int i = 0; i < data.length; i++) {
          list.add(new DoubleVec(data[i]));
        }
        Collections.sort(list);
        for (int i = 0; i < data.length; i++) {
          ma.assignRow(i, (DoubleVec) list.get(i));
        }
        return ma;
    }

    public Matrix getSQRT() {
        Matrix ma = new Matrix(getRowDim(), getColDim());
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length; j++) {
                ma.data[i][j] = Math.sqrt(data[i][j]);
            }
        }
        return ma;
    }

    public double getSquaredSum() {
        double s = 0;
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length; j++) {
                s += data[i][j] * data[i][j];
            }
        }
        return s;
    }

    public Matrix getStandardDeviationCols() {
        
        Matrix maMean = getMeanCols();

        final double [] arrMeanCols = maMean.getRow(0);
    	
    	final int cols = cols();
    	
    	final int rows = rows();
    	
    	final double [] arrStdvCols = new double [cols];
        
        for (int i = 0; i < rows; i++) {
            
            final double [] arrRow = data[i];
            for (int j = 0; j < cols; j++) {
            	arrStdvCols[j] += (arrRow[j] - arrMeanCols[j]) * (arrRow[j] - arrMeanCols[j]);
            }
        }
        
        for (int i = 0; i < cols; i++) {
        	double sdv = Math.sqrt(arrStdvCols[i] / (rows - 1));
        	
        	arrStdvCols[i] = sdv;
        }
        
        return new Matrix(true, arrStdvCols);
    }
    
    public double getStandardDeviation() {
        
    	double sdv=0;
    	
    	double n = rows()*cols();
    	
    	double mean = getMean();
    	
    	double sum=0;
    	for (int i = 0; i < data.length; i++) {
    		for (int j = 0; j < data[0].length; j++) {
    			sum += (data[i][j]-mean)*(data[i][j]-mean);
    		}
		}
    	
    	sdv = Math.sqrt(sum / (n - 1));
        
        return sdv;
    }

    public double getVariance() {
        double var = 0;

        double mean = getMean();
        double sum = 0;
        for (int j = 0; j < data[0].length; j++) {
            for (int i = 0; i < data.length; i++) {
                sum += (data[i][j] - mean) * (data[i][j] - mean);
            }
        }

        var = sum / (getNumElements() - 1);

        return var;
    }
    /**
     * Coefficient of variation (CV) is a normalized measure of dispersion of a probability distribution. 
     * It is defined as the ratio of the standard deviation sigma to the mean mu.
     * @return
     */
    public double getCoefficientVariation() {

        double mean = getMean();
        double sum = 0;
        for (int j = 0; j < data[0].length; j++) {
            for (int i = 0; i < data.length; i++) {
                sum += (data[i][j] - mean) * (data[i][j] - mean);
            }
        }

        double var = sum / (getNumElements() - 1);

        double sdv = Math.sqrt(var);
        
        double varCoeff=sdv/mean;
        
        return varCoeff;
    }

    public double getVarianceCol(int col) {
        double var = 0;

        double mean = getMeanCol(col);
        double dSum = 0;
        for (int i = 0; i < data.length; i++) {
        	dSum += (data[i][col] - mean) * (data[i][col] - mean);
        }

        var = dSum / (data.length - 1.0);

        return var;
    }

    public double getVarianceRow(int row) {
        double var = 0;

        int cols = cols();
        double mean = getMeanRow(row);
        double dSum = 0;
        for (int i = 0; i < cols; i++) {
        	dSum += (data[row][i] - mean) * (data[row][i] - mean);
        }

        var = dSum / (cols - 1.0);

        return var;
    }

    public double getVarianceCentered() {
        double var = 0;

        double sum = 0;
        for (int j = 0; j < data[0].length; j++) {
            for (int i = 0; i < data.length; i++) {
                sum += (data[i][j]) * (data[i][j]);
            }
        }

        var = sum / (getNumElements() - 1);

        return var;
    }

    /**
     * 
     * @return row with variance
     */
    public Matrix getVarianceCols() {
        Matrix ma = new Matrix(1, getColDim());
        Matrix maMean = getMeanCols();

        for (int j = 0; j < data[0].length; j++) {
        	
            double sum = 0;
            
            for (int i = 0; i < data.length; i++) {
                sum += (data[i][j] - maMean.data[0][j]) * (data[i][j] - maMean.data[0][j]);
            }
            
            double dVariance = sum / (getRowDim() - 1);
            
            ma.data[0][j] = dVariance;
        }
        return ma;
    }

    public void increase(int row, int col, double v){
    	data[row][col]+=v;
    }
    
    final public void set(final int row, final int col, final double v) {
        data[row][col] = v;
    }

    public void set(Matrix ma) {
        resize(ma.getRowDim(), ma.getColDim());
        for (int ii = 0; ii < data.length; ii++) {
            for (int jj = 0; jj < data[0].length; jj++) {
                data[ii][jj] = ma.data[ii][jj];
            }
        }
    }

    public void set(double [][] arr) {
        resize(arr.length, arr[0].length);
        for (int ii = 0; ii < data.length; ii++) {
            for (int jj = 0; jj < data[0].length; jj++) {
                data[ii][jj] = arr[ii][jj];
            }
        }
    }
    
    public void setFlat(double [][] arr) {
    	data = arr;
    }

    public void set(double v) {
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length; j++) {
                data[i][j] = v;
            }
        }
    }

    public void setID(int id) {
        identifier = id;
    }

    public void setCol(int col, double v) {
        for (int i = 0; i < data.length; i++)
            data[i][col] = v;
    }
    
    public void setCol(int col, double [] arr) {
        for (int i = 0; i < data.length; i++)
            data[i][col] = arr[i];
    }
    
    public void setRow(int iRow, double v) {
        for (int i = 0; i < data[0].length; i++)
            data[iRow][i] = v;
    }

    /**
     * Deep copy
     * @param iRow
     * @param arr
     */
    public void setRow(int iRow, double [] arr) {
        for (int i = 0; i < data[0].length; i++)
            data[iRow][i] = arr[i];
    }
    
    public void setRow(int iRow, int [] arr) {
        for (int i = 0; i < data[0].length; i++)
            data[iRow][i] = arr[i];
    }

    public static void setSeparatorCol(String s) {
        OUT_SEPARATOR_COL = s;
    }

    public static void setSeparatorRow(String s) {
        OUT_SEPARATOR_ROW = s;
    }


    public void shuffleRows() {
        int r = rows();
        List<Integer> li = new ArrayList<>(r);
        for (int i = 0; i < r; i++) {
            li.add(i);
        }
        Collections.shuffle(li);
        for (int i = 0; i < li.size(); i++) {
            swapRows(i, li.get(i));
        }
    }

    public void swapRows(int a, int b) {
        double [] t = data[a];
        data[a]=data[b];
        data[b]=t;
    }

    public void sortRows(int col){
        int r= rows();
        List<IntegerDouble> li = new ArrayList<>(r);
        for (int i = 0; i < r; i++) {
            li.add(new IntegerDouble(i, get(i, col)));
        }
        Collections.sort(li, IntegerDouble.getComparatorDouble());
        double [][] dataSorted = new double[r][];
        for (int i = 0; i < r; i++) {
            dataSorted[i]=data[li.get(i).getInt()];
        }
        data = dataSorted;
    }

    /**
     * Get the standard scores, also known as z-scores.
     * @return matrux with standardized values.
     */
    public Matrix getStandardized() {

        Matrix ma = new Matrix(getRowDim(), getColDim());
        Matrix maStandardDeviation = getStandardDeviationCols();
        Matrix Xc = getCenteredMatrix();

        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length; j++) {
                //double d = Xc.mData[ii][jj] / maStandardDeviation.mData[0][jj];
                ma.data[i][j] = Xc.data[i][j] / maStandardDeviation.data[0][j];
            }
        }
        return ma;
    }

    /**
     *
     * @param indexRowStart start row included
     * @param indexRowEnd end row included
     * @param indexColStart start col included
     * @param indexColEnd end col included
     * @return matrix
     */

    public Matrix getSubMatrix(int indexRowStart, int indexRowEnd, int indexColStart, int indexColEnd) {

        int rows = indexRowEnd - indexRowStart + 1;
        int cols = indexColEnd - indexColStart + 1;

        Matrix ma = new Matrix();

        double [][] arr = new double [rows][cols];
        
        for (int i = 0; i < rows; i++) {
        	System.arraycopy(data[indexRowStart+i], indexColStart, arr[i], 0, cols);
        }
        
        ma.data = arr;
        
        return ma;
    }
    
    public Matrix getSubMatrix(List<Integer> liIndexRow) {

        int cols = cols();

        Matrix ma = new Matrix();

        double [][] arr = new double [liIndexRow.size()][cols];
        
        int row=0;
        for (int i = 0; i < liIndexRow.size(); i++) {
        	System.arraycopy(data[liIndexRow.get(i)], 0, arr[row], 0, cols);
        	row++;
        }
        
        ma.data = arr;
        
        return ma;
    }

    public double getSum() {
        double sum = 0;
        for (int i = 0; i < getRowDim(); i++) {
            for (int j = 0; j < getColDim(); j++) {
                sum += data[i][j];
            }
        }
        return sum;
    }

    /**
     * Matrix has to be quadratic.
     * @return
     */
    public double getSumUpperTriangle() {
    	
    	if(rows()!=cols()){
    		throw new RuntimeException("Not a quadratic matrix.");
    	}
    	
        int rows = cols();

		double sum=0;
		for (int i = 0; i < rows; i++) {
			for (int j = i+1; j < rows; j++) {
				if(!Double.isNaN(get(i, j))) {
					sum += get(i, j);
				}
			}
		}

		return sum;
    }
    
    
    public Matrix getSumCols() {
        Matrix ma = new Matrix(1,getColDim());
        for (int i = 0; i < getRowDim(); i++) {
            for (int j = 0; j < getColDim(); j++) {
                ma.data[0][j] += data[i][j];
            }
        }
        return ma;
    }
    
    public Matrix getSumRows() {
        Matrix ma = new Matrix(1,getRowDim());
        for (int i = 0; i < getRowDim(); i++) {
           ma.set(0, i, getSumRow(i));
        }
        return ma;
    }

    public double getSumCol(int col) {
        double sum = 0;
        for (int i = 0; i < getRowDim(); i++) {
                sum += data[i][col];
            
        }
        return sum;
    }

    public double getSumRow(int row) {
        double sum = 0;
        for (int i = 0; i < getColDim(); i++) {
                sum += data[row][i];
            
        }
        return sum;
    }
    
    public double getSumSquared() {
        double sum = 0;
        for (int i = 0; i < getRowDim(); i++) {
            for (int j = 0; j < getColDim(); j++) {
                sum += (data[i][j] * data[i][j]);
            }
        }
        return sum;
    }

    public double [] getNext4NeighboursTorus(int row, int col) {
        double [] arr = new double [4];

        // value above
        arr[0] = getTorus(row - 1, col);
        // left
        arr[1] = getTorus(row, col - 1);
        // right
        arr[2] = getTorus(row, col + 1);
        // down
        arr[3] = getTorus(row + 1, col);

        return arr;
    }

    public double [] getNext8NeighboursTorus(int row, int col) {
        double [] arr = new double [8];

        int cc = 0;
        for (int i = -1; i < 2; i++) {
            for (int j = -1; j < 2; j++) {
                if(i != 0 || j != 0) {
                    arr[cc] = getTorus(row + i, col + j);
                    cc++;
                }
            }
        }
        return arr;
    }

    public Matrix row2Matrix(int row, int len) {
    	Matrix m = new Matrix(getColDim()/len, len);
    	
		int r = 0;
		int c = 0;
    	for (int i = 0; i < getColDim(); i++) {
    		if((i % len == 0) && (i > 0)) {
    			r++;
    			c = 0;
    		}
    		m.set(r,c, get(0,i));
			c++;
		}    	
    	
    	return m;
    }
    
    public double getTorus(int row, int col) {

        int torRow = row;
        int torCol = col;

        while(torRow < 0) {
            torRow += getRowDim();
        }

        while(torRow >= getRowDim()) {
            torRow -= getRowDim();
        }

        while(torCol < 0) {
            torCol += getColDim();
        }

        while(torCol >= getColDim()) {
            torCol -= getColDim();
        }

        return get(torRow, torCol);
    }

    /** Matrix trace.
       @return     sum of the diagonal elements.
     */
    public double getTrace() {
        double t = 0;

        int m = getRowDim();
        int n = getColDim();

        for (int ii = 0; ii < Math.min(m, n); ii++) {
            t = t + get(ii, ii);
        }
        return t;
    }

    public Matrix getTranspose() {
        int rows = getRowDim();
        int cols = getColDim();
        Matrix ma = new Matrix(cols, rows);
        for (int ii = 0; ii < rows; ii++) {
            for (int jj = 0; jj < cols; jj++) {
                ma.data[jj][ii] = data[ii][jj];
            }
        }
        return ma;
    }
    /**
     * Opposite direction to transpose.
     */ 
    public Matrix getFlipped() {
        int rows = getRowDim();
        int cols = getColDim();
        Matrix ma = new Matrix(cols, rows);
        for (int ii = 0; ii < rows; ii++) {
            for (int jj = 0; jj < cols; jj++) {
                ma.data[cols-jj-1][ii] = data[ii][jj];
            }
        }
        return ma;
    }

    public Matrix toRow() {
    	Matrix m = new Matrix(1, getRowDim() * getColDim());
    	for (int i = 0; i < getRowDim(); i++) {
        	for (int j = 0; j < getColDim(); j++) {
    			m.set(0,i * getColDim() + j, get(i,j));
    		}
		}
    	return m;
    }
    
    /**
     * 
     * @return deep copy.
     */
    public double [] toArray() {
    	
    	final int rows = rows();
    	final int cols = cols();
    	
    	double [] a = new double [rows*cols];
    	
    	for (int i = 0; i < rows; i++) {
        	for (int j = 0; j < cols; j++) {
    			a[i * cols + j] = get(i,j);
    		}
		}
    	
    	return a;
    }
    
    public String toString() {

        int iRequireDigits = 20;
        int len = getRowDim() * getColDim() * iRequireDigits;
        StringBuilder sb = new StringBuilder(len);

        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length - 1; j++) {
                sb.append(data[i][j] + OUT_SEPARATOR_COL);
            }
            sb.append(data[i][data[0].length - 1]);
            // No line break after the last line
            if(i < data.length - 1)
                sb.append(OUT_SEPARATOR_ROW);
        }


        return sb.toString();
    }
    
    public String toStringBinary() {

        int len = getRowDim() * getColDim();
        StringBuilder sb = new StringBuilder(len);

        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length; j++) {
            	if(data[i][j]==0)
            		sb.append(0);
            	else
            		sb.append(1);
            }
            // No line break after the last line
            if(i < data.length - 1)
                sb.append(OUT_SEPARATOR_ROW);
        }


        return sb.toString();
    }
    public String toString(int digits) {
        return toString(rows(), cols(), digits, 0);
    }
    public String toString(int digits, int width) {
        return toString(rows(), cols(), digits, width);
    }

    /**
     *
     * @param rowEnd exclusive
     * @param colEnd exclusive
     * @param digits
     * @return
     */
    public String toString(int rowEnd, int colEnd, int digits, int width) {

        int iRequireDigits = 20;

        String sFormat = "";

        sFormat += "0";
        int iCounter = 0;
        if(digits > 0)
            sFormat += ".";

        while(iCounter < digits) {
          sFormat += "0";
          iCounter++;
        }

        DecimalFormat nf = new DecimalFormat(sFormat, new DecimalFormatSymbols(Locale.US));

        int len = getRowDim() * getColDim() * iRequireDigits;
        StringBuilder sb = new StringBuilder(len);

        
        for (int i = 0; i < rowEnd; i++) {
            for (int j = 0; j < colEnd; j++) {
            	
            	String sVal = nf.format(data[i][j]);
            	if(data[i][j]==Double.MAX_VALUE)
            		sVal = "Max";

                StringBuilder sbVal = new StringBuilder(sVal);
                while (sbVal.length()<width){
                    sbVal.insert(0, " ");
                }
                sb.append(sbVal.toString());
                
                if(j<data[0].length-1)
                	sb.append(OUT_SEPARATOR_COL);
            }
            
            if(i < data.length - 1)
                sb.append(OUT_SEPARATOR_ROW);
        }

        return sb.toString();
    }

    public String toStringRow(int row, int iDigits) {
        int iRequireDigits = 20;

        String sFormat = "";

        sFormat += "0";
        int iCounter = 0;
        if(iDigits > 0)
            sFormat += ".";

        while(iCounter < iDigits) {
          sFormat += "0";
          iCounter++;
        }

        DecimalFormat nf = new DecimalFormat(sFormat);

        int len = getRowDim() * getColDim() * iRequireDigits;
        StringBuilder sb = new StringBuilder(len);


        for (int j = 0; j < data[0].length; j++) {

            String sVal = nf.format(data[row][j]);
            if(data[row][j]==Double.MAX_VALUE)
                sVal = "Max";

            sb.append(sVal);

            if(j<data[0].length-1)
                sb.append(OUT_SEPARATOR_COL);
        }


        return sb.toString();
    }

    public String toStringRowNumber(int iDigits, String sSeparatorCol) {
        int iRequireDigits = 20;

        String sFormat = "##,";

        int iCounter = 0;
        while(iCounter < iDigits) {
          sFormat += "#";
          iCounter++;
        }
        sFormat += "0.";
        iCounter = 0;
        while(iCounter < iDigits) {
          sFormat += "0";
          iCounter++;
        }

        DecimalFormat nf = new DecimalFormat(sFormat);

        int len = getRowDim() * getColDim() * iRequireDigits;
        StringBuilder sb = new StringBuilder(len);

        for (int i = 0; i < data.length; i++) {
            sb.append(i + sSeparatorCol);
            for (int j = 0; j < data[0].length - 1; j++) {
                sb.append(nf.format(data[i][j]) + sSeparatorCol);
            }
            sb.append(nf.format(data[i][data[0].length - 1]));
            if(i < data.length - 1)
                sb.append(OUT_SEPARATOR_ROW);
        }

        return sb.toString();
    }

    public Matrix subtract(Matrix ma) {
        Matrix maResult = new Matrix(getRowDim(), getColDim());
        if(!equalDimension(ma)) {
            throw new RuntimeException("Matrices have wrong dimensions.");
        }

        for (int ii = 0; ii < getRowDim(); ii++) {
            for (int jj = 0; jj < getColDim(); jj++) {
                maResult.data[ii][jj] = data[ii][jj] - ma.data[ii][jj];
            }
        }
        return maResult;
    }
    
    public Matrix subtract(double val) {
        Matrix maResult = new Matrix(getRowDim(), getColDim());

        for (int ii = 0; ii < getRowDim(); ii++) {
            for (int jj = 0; jj < getColDim(); jj++) {
                maResult.data[ii][jj] = data[ii][jj] - val;
            }
        }
        return maResult;
    }

    public Matrix subtractFromCols(Matrix maDivisor) {
        Matrix maResult = new Matrix(getRowDim(), getColDim());

        for (int ii = 0; ii < getRowDim(); ii++) {
            for (int jj = 0; jj < getColDim(); jj++) {
                maResult.data[ii][jj] = data[ii][jj] - maDivisor.get(0,jj);
            }
        }
        return maResult;
    }

    public static double Sign(double a, double b) {
        // ( (b) >= 0.0 ? Math.abs(a) : -Math.abs(a))

        if(b >= 0.0)
            return Math.abs(a);
        else
            return -Math.abs(a);

    }
    
	public static Matrix getRND(int rows, int cols){
		
		Matrix ma = new Matrix(rows, cols);
		
		for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                ma.set(i,j,Math.random());
		
		return ma;
	}

    public void write(String sFile, boolean apppend) {
    	
    	DecimalFormat df = new DecimalFormat(FORMAT);
    	
    	write(new File(sFile), df, apppend);
    }
    
    public void write(String sFile) {
    	DecimalFormat df = new DecimalFormat(FORMAT);
    	
    	write(new File(sFile), df, false);
    }
    
    public void write(File fiMa) {
    	
    	DecimalFormat df = new DecimalFormat(FORMAT);
    	
    	write(fiMa, df, false);
    }
    
    public void write(File fiMa, boolean apppend, int digits) throws IOException{
        String sFormat = "##,";
        int iCounter = 0;
        while(iCounter < digits) {
          sFormat += "#";
          iCounter++;
        }
        sFormat += "0.";
        iCounter = 0;
        while(iCounter < digits) {
          sFormat += "0";
          iCounter++;
        }

        DecimalFormat nf = new DecimalFormat(sFormat);

        write(fiMa,  nf, apppend);
    }
    
    public void write(String sFiMa, boolean apppend, int digits) throws IOException{
    	write(new File (sFiMa), apppend, digits);
    }

    public void write(File fiMa, DecimalFormat nf, boolean apppend) {

        try {
            FileOutputStream os = new FileOutputStream(fiMa, apppend);
            write(os, nf);
            os.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
    }

    public void write(OutputStream os) {
        NumberFormat nf = new DecimalFormat("#.###############");
        write(os, nf);
    }

    public String writeAsLineBase64Encoded() {

        NumberFormat nf = new DecimalFormat("#.###############");

        ByteArrayOutputStream baos = new ByteArrayOutputStream();

        write(baos, nf);

        Base64.Encoder encoder = Base64.getEncoder();

        byte [] arr64 = encoder.encode(baos.toByteArray());

        return new String(arr64);
    }


    public void write(OutputStream os, NumberFormat nf) {

        try {
			for (int i = 0; i < data.length; i++) {

			    StringBuilder sb =  new StringBuilder();

			    for (int j = 0; j < data[0].length; j++) {

			    	sb.append(nf.format(data[i][j]));

			    	if(j<data[0].length-1){
			    		sb.append(OUT_SEPARATOR_COL);
			    	}
			    }

			    os.write(sb.toString().getBytes());
			    // No line break after the last line
			    if(i < data.length - 1) {
                    os.write(OUT_SEPARATOR_ROW.getBytes());
			    }
			}

			os.flush();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
    }


    public String toStringWithColTags(List<String> liColTags, DecimalFormat nf, String separator) {
    	    	
    	if(cols()!=liColTags.size()){
    		throw new RuntimeException("Number of cols and col tags differ.");
    	}
    	
    	int [] arrWidth = new int [cols()];
    	
    	for (int i = 0; i < arrWidth.length; i++) {
			
    		arrWidth[i] = liColTags.get(i).length();
    		
		}
    	
    	for (int i = 0; i < arrWidth.length; i++) {
    		
    		for (int j = 0; j < rows(); j++) {
    			arrWidth[i] = Math.max(arrWidth[i], nf.format(get(j,i)).length());
			}
    	}
    	
    	StringBuilder sbAll = new StringBuilder();
    	
    	for (int i = 0; i < arrWidth.length; i++) {
    		
    		int w = arrWidth[i];
    		
    		StringBuilder sb = new StringBuilder(liColTags.get(i));
    		
    		int l = w-sb.length();
    		
    		for (int j = 0; j < l; j++) {
    			sb.append(" ");
			}
    		
    		sbAll.append(sb.toString());
    		sbAll.append(" ");
    		
    	}
    	
    	sbAll.append("\n");
    	
    	
    	for (int i = 0; i < rows(); i++) {
    		
    		for (int j = 0; j < arrWidth.length; j++) {
        		
        		int w = arrWidth[j];
        		
        		StringBuilder sb = new StringBuilder(nf.format(get(i,j)));
        		
        		int l = w-sb.length();
        		
        		for (int k = 0; k < l; k++) {
        			sb.append(" ");
    			}
        		
        		sbAll.append(sb.toString());
        		sbAll.append(" ");
        		
        	}
    		
    		sbAll.append("\n");
		}
    	
    	return sbAll.toString();
    }
    
    public String toStringWithRowTags(List<String> liRowTags, DecimalFormat nf, String separator) {
    	
    	if(rows()!=liRowTags.size()){
    		throw new RuntimeException("Number of rows and row tags differ.");
    	}
    	
    	int [] arrWidth = new int [cols()+1];
    	
    	for (int i = 0; i < liRowTags.size(); i++) {
			
    		arrWidth[0] = Math.max(arrWidth[0], liRowTags.get(i).length());
    		
		}
    	
    	for (int i = 0; i < cols(); i++) {
    		
    		for (int j = 0; j < rows(); j++) {
    			arrWidth[i+1] = Math.max(arrWidth[i+1], nf.format(get(j,i)).length());
			}
    	}
    	
    	StringBuilder sbAll = new StringBuilder();
    	
    	for (int i = 0; i < arrWidth.length; i++) {
    		
    		int w = arrWidth[i];
    		
    		StringBuilder sb = new StringBuilder(liRowTags.get(i));
    		
    		int l = w-sb.length();
    		
    		for (int j = 0; j < l; j++) {
    			sb.append(" ");
			}
    		
    		sbAll.append(sb.toString());
    		sbAll.append(" ");
    		
    	}
    	
    	sbAll.append("\n");
    	
    	
    	for (int i = 0; i < rows(); i++) {
    		
    		for (int j = 0; j < arrWidth.length; j++) {
        		
        		int w = arrWidth[j];
        		
        		StringBuilder sb = new StringBuilder(nf.format(get(i,j)));
        		
        		int l = w-sb.length();
        		
        		for (int k = 0; k < l; k++) {
        			sb.append(" ");
    			}
        		
        		sbAll.append(sb.toString());
        		sbAll.append(" ");
        		
        	}
    		
    		sbAll.append("\n");
		}
    	
    	return sbAll.toString();
    }

	public void writeSerialized(File fiOut) throws IOException {
		FileOutputStream fos = new FileOutputStream(fiOut);
		ObjectOutputStream oos = new ObjectOutputStream(fos);

		
		oos.writeObject(data);
		
		
		oos.close();
	}
	
	public static Matrix readSerialized(File fiIn) throws FileNotFoundException, IOException, ClassNotFoundException {
		
		Matrix ma = null;
		
		FileInputStream fos = new FileInputStream(fiIn);
		
		ObjectInputStream ois = new ObjectInputStream(fos);
		
		double [][] data = (double [][]) ois.readObject();
		
		/**
		 * We are using a deep copy because of memory handling problems with using the de-serialized object directly.
		 */
		ma = new Matrix(data);
		
		ois.close();
		
		return ma;
	}

    public static DecimalFormat format(int digits) {
        String sFormat = "##,###";
        
        
        if(digits > 0) {
        	sFormat += ".";
	        int iCounter = 0;
	        while (iCounter < digits) {
	            sFormat += "0";
	            iCounter++;
	        }
        }
        return new DecimalFormat(sFormat, new DecimalFormatSymbols(Locale.US));
    }

    public void write(String sFile, boolean bApppend, int digits, int totalWidth) {

        try {
            DecimalFormat nf = format(digits);
            BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
                sFile), bApppend));

            StringBuffer sVal;
            for (int ii = 0; ii < data.length; ii++) {
                sVal = new StringBuffer();
                for (int jj = 0; jj < data[0].length; jj++) {
                    sVal.append(format(data[ii][jj], nf, totalWidth) +
                                OUT_SEPARATOR_COL);
                }
                writer.write(sVal.toString());
                // No line break after the last line
                if (ii < data.length - 1) {
                    writer.write(OUT_SEPARATOR_ROW);
                }
            }
            writer.flush();
            writer.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    public static String format(double val, DecimalFormat nf, int totalWidth) {
        String str = "";
        str = nf.format(val);
        while(str.length() < totalWidth) {
            str = " " + str;
        }
        return str;
    }

}