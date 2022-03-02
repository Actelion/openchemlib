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

package com.actelion.research.calc;


import com.actelion.research.chem.descriptor.SimilarityCalculatorDoubleArray;
import com.actelion.research.util.DoubleVec;
import com.actelion.research.util.Formatter;
import com.actelion.research.util.datamodel.*;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.Base64.Decoder;
import java.util.List;

/**
 * <p>Title: MatrixFunctions</p>
 * <p>Description: Matrix operations for which a direct access to the matrix
 * member variables is not necessary</p>
 *
 * References
 * http://www.ini.uzh.ch/~fred/java/Matrix.java
 * https://introcs.cs.princeton.edu/java/95linear/Cholesky.java.html
 *
 * @author Modest von Korff
 * @version 1.0
 * 20.11.2003 MvK: Start implementation
 * 10.06.2004 MvK read matrix.
 */

public class MatrixFunctions {

    private static final double TINY = 0.0000001;




    public static Matrix convert2Binary(Matrix A, double thresh) {
        int rows = A.rows();
        int cols = A.cols();
        Matrix B = new Matrix(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if(A.get(i,j)>thresh){
                    B.set(i,j,1);
                }
            }
        }
        return B;
    }

    /**
     * Splits the matrix row wise into two matrices.
     * @param A
     * @param row
     * @return
     */
    public static Matrix [] split(Matrix A, int row) {

        int rows = A.rows();
        int cols = A.cols();
        int rowsC = rows-row;

        Matrix B = new Matrix(row, cols);
        Matrix C = new Matrix(rowsC, cols);

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < cols; j++) {
                B.set(i, j, A.get(i,j));
            }
        }

        for (int i = 0; i < rowsC; i++) {
            for (int j = 0; j < cols; j++) {
                C.set(i, j, A.get(i+row,j));
            }
        }

        Matrix [] arr = new Matrix[2];

        arr[0]=B;
        arr[1]=C;

        return arr;
    }

    public static Matrix [] splitCol(Matrix A, int col) {

        int rows = A.rows();
        int cols = A.cols();
        int colsC = cols-col;

        Matrix B = new Matrix(rows, col);
        Matrix C = new Matrix(rows, colsC);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < col; j++) {
                B.set(i, j, A.get(i,j));
            }
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < colsC; j++) {
                C.set(i, j, A.get(i, j+col));
            }
        }

        Matrix [] arr = new Matrix[2];

        arr[0]=B;
        arr[1]=C;

        return arr;
    }

    public static boolean isEqualCol(Matrix A, int colA, Matrix B, int colB, double thresh) {

        boolean equal = true;

        int r = A.rows();

        if(r != B.rows()){
            throw new RuntimeException("Row number differs!");
        }

        for (int i = 0; i < r; i++) {

            double abs = Math.abs(A.get(i,colA)-B.get(i,colB));

            if(abs>thresh){
                equal=false;
                break;
            }
        }

        return equal;
    }


    public static boolean equals(double d1, double d2){
        return (Math.abs(d1-d2)<TINY)?true:false;
    }

    /**
     * Returns a (n x n) identity matrix.
     * (The elements on the diagonal have value <code>1.0</code>, all others are set to <code>0.0</code>).
     *
     * @param n the size of the matrix
     * @return the identity matrix
     */
    public static double[][] id(int n){
        double[][] A = new double [n][n];
        for(int i=0; i<n; i++)
            for(int j=0; j<n; j++){
                if (i==j) A[i][j] = 1.0;
                else A[i][j] = 0.0;
            }
        return A;
    }

    /**
     * Returns the inverse of a matrix. The method uses a LU decomposition with pivot.
     * There is no check on the size (is it a square matrix) or the rank (is it inversible, ie det doesn't equal 0).
     * @param A
     * @return A^(-1)
     */
    public static Matrix inv(Matrix A){
// based on the chapter IV of the course "Analyse numerique" taugth by Pr Ernst Hairer at the University of Geneva

        int n = A.rows();
        double[][] invertedA = new double[n][n];
        double[][] R = new double[n][n];

        double[][] basis = id(n);

        // R = copy of A (the matrix on which we make the operations, so we leave A unchanged)
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                R[i][j] = A.get(i,j);
            }
        }
        // Triangulation of R
        for(int i=0; i<n-1; i++){
            // Find pivot and swap lines
            double a = Math.abs(R[i][i]);
            for (int j = i+1; j < n; j++) {
                if(Math.abs(R[j][i]) > a){
                    a = Math.abs(R[j][i]);
                    double[] tempLine = R[j]; R[j] = R[i] ; R[i] = tempLine;
                    tempLine = basis[j];
                    basis[j] = basis[i] ;
                    basis[i] = tempLine;
                }
            }
            // Elimination of (i+i)th element of each line >i in R
            for(int j = i+1; j<n; j++){
                double l = R[j][i]/R[i][i];
                for(int k = i; k<n; k++){  // CAUTION we start at k = i !!
                    R[j][k] = R[j][k] - l*R[i][k];
                }
                for(int k = 0; k<n; k++){  // CAUTION we start at k = 0 !!
                    basis[j][k] = basis[j][k] - l*basis[i][k];
                }
            }
        }

        for (int k = 0; k < basis.length; k++) {
            for (int i = n - 1; i > -1; i--) {
                double sum = 0.0;
                for (int j = i + 1; j < n; j++)
                    sum += R[i][j] * invertedA[j][k];
                invertedA[i][k] = (basis[i][k] - sum) / R[i][i];
            }
        }
        return new Matrix(invertedA);
    }


    public static boolean isSymmetric(Matrix A) {
        int rows = A.rows();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < i; j++) {
                if (A.get(i,j) != A.get(j,i)) return false;
            }
        }
        return true;
    }

    public static boolean isSquare(Matrix A) {
        int rows = A.rows();
        int cols = A.cols();
        return (rows==cols)?true:false;
    }


    /**
     *
     * @param A
     * @return Cholesky factor L of psd matrix A = L L^T
     */
    public static Matrix cholesky(Matrix A) {
        if (!isSquare(A)) {
            throw new RuntimeException("Matrix is not square");
        }
        if (!isSymmetric(A)) {
            throw new RuntimeException("Matrix is not symmetric");
        }

        int rows  = A.rows();
        double[][] L = new double[rows][rows];

        for (int i = 0; i < rows; i++)  {
            for (int j = 0; j <= i; j++) {
                double sum = 0.0;
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * L[j][k];
                }
                if (i == j)
                    L[i][i] = Math.sqrt(A.get(i,i) - sum);
                else
                    L[i][j] = 1.0 / L[j][j] * (A.get(i,j) - sum);
            }
            if (L[i][i] <= 0) {
                throw new RuntimeException("Matrix not positive definite");
            }
        }

        return new Matrix(L);
    }

    public static Matrix calculateSimilarityMatrixRowWise(Matrix X1, Matrix X2){

        SimilarityCalculatorDoubleArray similarityCalculatorDoubleArray = new SimilarityCalculatorDoubleArray();


        int r1 = X1.rows();
        int r2 = X2.rows();

        Matrix maSim = new Matrix(r1, r2);

        for (int i = 0; i < r1; i++) {

            double [] x1 = X1.getRow(i);

            for (int j = 0; j < r2; j++) {

                double [] x2 = X2.getRow(j);

                double v = similarityCalculatorDoubleArray.getSimilarity(x1, x2);

                maSim.set(i,j,v);
            }
        }

        return maSim;
    }


    public static double [] calculateMaxSimilarity(Matrix XTrain, Matrix Xtest){

        double [] arrSim = new double[Xtest.rows()];

        int r = Xtest.rows();

        for (int i = 0; i < r; i++) {

            double [] arrDescriptor = Xtest.getRow(i);

            IntegerDouble maxSim = MatrixFunctions.calculateMaxSimilarity(XTrain, arrDescriptor);

            arrSim[i] = maxSim.getDouble();

        }

        return arrSim;
    }


    public static IntegerDouble[] calculateMaxSimilarity(Matrix X, double [] arr, int nMostSimilar){

        List<IntegerDouble> li = new ArrayList<>();

        int rows = X.rows();

        for (int i = 0; i < rows; i++) {
            double similarity = DoubleVec.getTanimotoSimilarity(X.getRow(i), arr);
            li.add(new IntegerDouble(i, similarity));
        }

        Collections.sort(li, IntegerDouble.getComparatorDouble());
        Collections.reverse(li);

        int n = Math.min(li.size(), nMostSimilar);

        IntegerDouble [] arrMax = new IntegerDouble[n];

        for (int i = 0; i < n; i++) {

            arrMax[i]=li.get(i);

        }

        return arrMax;
    }

    public static IntegerDouble calculateMaxSimilarity(Matrix X, double [] arr){

        int rows = X.rows();

        IntegerDouble maxSim = new IntegerDouble(-1, 0);
        for (int i = 0; i < rows; i++) {
            double similarity = DoubleVec.getTanimotoSimilarity(X.getRow(i), arr);

            if(similarity>maxSim.getDouble()){
                maxSim.setInteger(i);
                maxSim.setDouble(similarity);
            }
        }

        return maxSim;
    }


    /**
     *
     * @param maQuadratic
     * @return upper triangle from a quadratic matrix as an array.
     */
    public static double [] upperTriangle(Matrix maQuadratic){

        int r = maQuadratic.rows();
        int n = ((r * r)-r) / 2;

        double [] arr= new double[n];

        int cc=0;

        for (int i = 0; i < r; i++) {

            for (int j = i+1; j < r; j++) {
                arr[cc++]=maQuadratic.get(i,j);
            }
        }

        return arr;
    }

    public static Matrix appendRows(List<Matrix> liMatrix) {

        if(liMatrix==null || liMatrix.size()==0){
            return null;
        }

        int cols = liMatrix.get(0).cols();

        int rows = 0;

        for (Matrix m : liMatrix) {
            rows+=m.rows();
        }

        Matrix mAll = new Matrix(rows, cols);

        int offsetRow=0;
        for (Matrix m : liMatrix) {
            mAll.copy(offsetRow, m);
            offsetRow+=m.rows();
        }


        return mAll;
    }


    public static Matrix appendRows(Matrix ma0, Matrix ma1) {

        Matrix ma = new Matrix(ma0.rows() + ma1.rows(), ma0.cols());

        for (int i = 0; i < ma0.rows(); i++) {
            for (int j = 0; j < ma0.cols(); j++) {
                ma.set(i,j, ma0.get(i,j));
            }
        }

        int offsetRows = ma0.rows();

        for (int i = 0; i < ma1.rows(); i++) {
            for (int j = 0; j < ma1.cols(); j++) {
                ma.set(i+offsetRows,j, ma1.get(i,j));
            }
        }

        return ma;
    }

    public static Matrix appendCols(Matrix ma0, Matrix ma1) {

        if(ma0.rows() != ma1.rows()){
            throw new RuntimeException("Number rows differ!");
        }

        Matrix ma = new Matrix(ma0.rows(), ma0.cols()+ma1.cols());

        for (int i = 0; i < ma0.rows(); i++) {
            for (int j = 0; j < ma0.cols(); j++) {
                ma.set(i,j, ma0.get(i,j));
            }
        }

        int offsetCols = ma0.cols();

        for (int i = 0; i < ma1.rows(); i++) {
            for (int j = 0; j < ma1.cols(); j++) {
                ma.set(i,j+offsetCols, ma1.get(i,j));
            }
        }

        return ma;
    }

	/**
	 * List is already used as a constructor for Matrix. So we have to place this method here.
	 * @param li
	 * @return
	 */
	public static Matrix create(List<double []> li){

		double [][] a = new double[li.size()][];

		for (int i = 0; i < li.size(); i++) {
			a[i]=li.get(i);
		}

		return new Matrix(a);

	}

    public static int countFieldsBiggerThan(Matrix ma, int row, double thresh) {
    	int cc = 0;
    	for (int i = 0; i < ma.getColDim(); i++) {
			if(ma.get(row,i) > thresh) {
				cc++;
			}
		}
    	return cc;
    }
    
    public static int countFieldsBiggerThanThreshColWise(Matrix ma, int col, double thresh) {
    	int cc = 0;
    	for (int i = 0; i < ma.rows(); i++) {
			if(ma.get(i,col) > thresh) {
				cc++;
			}
		}
    	return cc;
    }
    
    public static Matrix countFieldsBiggerThanThreshColWise(Matrix ma, double thresh) {
    	Matrix maCounts = new Matrix(1, ma.cols());
    	
    	for (int i = 0; i < ma.cols(); i++) {
			int cc = countFieldsBiggerThanThreshColWise(ma, i, thresh);
			
			maCounts.set(0, i, cc);
		}
    	
    	return maCounts;
    }
    
    /**
     * Calculates the inverse Tanimoto coefficient from row wise comparison of the two input matrices.
     * 
     * @param ma1
     * @param ma2
     * @return complete distance matrix calculated between all rows from  the two input matrices.
     */
    public static Matrix getDistanceMatrixTanimotoInv(Matrix ma1, Matrix ma2) {
        Matrix maDist = new Matrix(ma1.getRowDim(), ma2.getRowDim());
        for (int i = 0; i < ma1.getRowDim(); i++) {
            for (int j = 0; j < ma2.getRowDim(); j++) {
                double dist = getDistanceTanimotoInv(ma1, i, ma2, j);
                maDist.set(i,j, dist);
            }
        }
        return maDist;
    }
    
    public static Matrix getDistanceMatrix(List<Point> li) {
        Matrix maDist = new Matrix(li.size(), li.size());
        for (int i = 0; i <li.size(); i++) {
            for (int j = 0; j < li.size(); j++) {
                double dist = li.get(i).distance(li.get(j));
                maDist.set(i,j, dist);
            }
        }
        return maDist;
    }

    public static Matrix getDistTanimotoInvReduced(Matrix ma1, Matrix ma2) {
        Matrix maDist = new Matrix(ma1.getRowDim(), ma2.getRowDim());
        for (int i = 0; i < ma1.getRowDim(); i++) {
            for (int j = 0; j < ma2.getRowDim(); j++) {
                double dist = getDistTanimotoInvReduced(ma1, i, ma2, j);
                maDist.set(i,j, dist);
            }
        }
        return maDist;
    }
    /**
     *
     * @param ma1 Matrix
     * @param row1 row
     * @param ma2 Matrix
     * @param row2 row
     * @return maximum distance = 0, minimum distance = 1
     */
    public static double getDistanceTanimotoInv(Matrix ma1, int row1, Matrix ma2, int row2) {
        double dist = 0;
        double dAtB = ma1.multiply(row1, ma2, row2);
        double dAtA = ma1.multiply(row1, ma1, row1);
        double dBtB = ma2.multiply(row2, ma2, row2);
        dist = dAtB / (dAtA + dBtB - dAtB);
        return 1-dist;
    }

    public static double getSumSquaredDiff(Matrix ma1, int row1, Matrix ma2, int row2) {

        double sum = 0;

        for (int i = 0; i < ma1.cols(); i++) {
            double diff = ma1.get(row1, i) - ma2.get(row2, i);
            sum += diff*diff;
        }

        return sum;
    }

    public static double getTotalSumSquaredDiffRowWise(Matrix ma1, Matrix ma2) {

        int rows = ma1.rows();

        double sum = 0;
        for (int i = 0; i < rows; i++) {
            double dist = getSumSquaredDiff(ma1, i, ma2, i);

            sum += dist;

            System.out.println(Formatter.format1(dist));
        }

        return sum;
    }




    /**
     * Only fields are considered which are not nut 0 in both or in one of the
     * rows we from where the distance is calculated.
     * @param ma1 Matrix
     * @param row1 row
     * @param ma2 Matrix
     * @param row2 row
     * @return maximum distance = 0, minimum distance = 1
     */
    public static double getDistTanimotoInvReduced(Matrix ma1, int row1, Matrix ma2, int row2) {
        double dist = 0;

        List<Double> li1 = new ArrayList<Double>();
        List<Double> li2 = new ArrayList<Double>();
        for (int ii = 0; ii < ma1.getColDim(); ii++) {
            if((ma1.get(row1, ii) != 0) || (ma2.get(row2, ii) != 0)) {
                li1.add(new Double(ma1.get(row1, ii)));
                li2.add(new Double(ma2.get(row2, ii)));
            }
        }

        Matrix maRow1 = new Matrix(true, li1);
        Matrix maRow2 = new Matrix(true, li2);

        double dAtB = maRow1.multiply(0, maRow2, 0);
        double dAtA = maRow1.multiply(0, maRow1, 0);
        double dBtB = maRow2.multiply(0, maRow2, 0);
        dist = dAtB / (dAtA + dBtB - dAtB);
        return 1-dist;
    }


    /**
     * A square-shaped neighborhood that can be used to define a set of cells surrounding a given point.
     * 
     * @param p
     * @param ma
     * @return
     */
    public static List<Point> getMooreNeighborhood(Point p, Matrix ma) {
    	List<Point> li = new ArrayList<Point>();
    	
    	int startX = Math.max(0, p.x - 1);
    	int endX = Math.min(ma.cols(), p.x+2);
    	
    	int startY = Math.max(0, p.y - 1);
    	int endY = Math.min(ma.rows(), p.y+2);
    	
    	for (int i = startY; i < endY; i++) {
			for (int j = startX; j < endX; j++) {
				if(i!=p.y || j!=p.x) {
					if(ma.get(i,j)>0) {
						li.add(new Point(j,i));
					}
				}
			}
		}
    	return li;
    }
    
    public static List<Point> getMooreNeighborhood(Point p, int r, Matrix ma) {
    	List<Point> li = new ArrayList<Point>();
    	
    	int startX = Math.max(0, p.x - r);
    	int endX = Math.min(ma.cols(), p.x+r+1);
    	
    	int startY = Math.max(0, p.y - r);
    	int endY = Math.min(ma.rows(), p.y+r+1);
    	
    	for (int i = startY; i < endY; i++) {
			for (int j = startX; j < endX; j++) {
				if(i!=p.y || j!=p.x) {
					if(ma.get(i,j)>0) {
						li.add(new Point(j,i));
					}
				}
			}
		}
    	return li;
    }


    /**
	 * 
	 * @param ma
	 *            Matrix with ma.rows and 1 col. Containing the min value from
	 *            each row.
	 * @return
	 */
    public static Matrix getRowMinUnique(Matrix ma) {
    	
    	Matrix maTmp = new Matrix(ma);
    	
    	Matrix maMin = new Matrix(maTmp.getRowDim(), 1);
    	
    	for (int i = 0; i < maTmp.getRowDim(); i++) {
        	ScorePoint td = maTmp.getMinPos();
        	maMin.set(td.y,0,td.getScore());
        	maTmp.setRow(td.y, Double.MAX_VALUE);
        	maTmp.setCol(td.x, Double.MAX_VALUE);
		}
    	return maMin;
    }
    
    /**
     * 
     * @param ma
     * @return list with indices Point(row,col) for values>0.
     */
	public static List<Point> getPoints(Matrix ma){
		List<Point> li = new ArrayList<Point>();
		
		for (int i = 0; i < ma.rows(); i++) {
			for (int j = 0; j < ma.cols(); j++) {
				if(ma.get(i, j)>0){
					li.add(new Point(j,i));
				}
			}
		}
		
		return li;
	}
	
	public static List<Point> getIndicesUniqueMaxRowWise(Matrix maIn){
		List<Point> li = new ArrayList<Point>();
		
		Matrix ma = new Matrix(maIn);
		
		for (int i = 0; i < ma.rows(); i++) {
			
			Point p = ma.getMaxIndex();
			
			li.add(p);
			
			int row = p.y;
			
			for (int j = 0; j < ma.cols(); j++) {
				ma.set(row, j, -Double.MAX_VALUE);
			}
		}
		
		return li;
	}
	
	public static List<Point> getIndicesUniqueMaxColumnWise(Matrix maIn){
		List<Point> li = new ArrayList<Point>();
		
		Matrix ma = new Matrix(maIn);
		
		for (int col = 0; col < ma.cols(); col++) {
			
			Point p = ma.getMaxIndex();
			
			li.add(p);
			
			int colMax = p.x;
			
			for (int row = 0; row < ma.rows(); row++) {
				ma.set(row, colMax, -Double.MAX_VALUE);
			}
		}
		
		return li;
	}

    public static Matrix getScaledByFactor(Matrix ma, double fac) {
    	int cols = (int)(ma.cols()*fac+0.5);
    	int rows = (int)(ma.rows()*fac+0.5);
    	
    	Matrix maSc = new Matrix(rows, cols);
    	
    	for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				double v = ma.get((int)(i/fac), (int)(j/fac));
				maSc.set(i,j,v);
			}
		}
    	
    	return maSc;
    	
    }

    /**
     * Centers and divides by standard deviation.
     * @param ma
     * @return
     */
    public static Matrix getScaled(Matrix ma) {

	    Matrix maCentered = ma.getCenteredMatrix();

	    Matrix maSD = ma.getStandardDeviationCols();

	    int rows = ma.rows();
	    int cols = ma.cols();

        Matrix maScale = new Matrix(rows, cols);

        for (int i = 0; i < cols; i++) {

            double sdv = maSD.get(0,i);

            for (int j = 0; j < rows; j++) {
                maScale.set(j,i, maCentered.get(j,i) / sdv);
            }
        }

        return maScale;
    }



    /**
     * 03.10.04 MvK
     * @param ma matrix with objects in rows
     * @param k number of desired cluster
     * @return Matrix of cluster centers, k rows and cols equal ma.cols.
     */
    public static Matrix getKMeanClusters(Matrix ma, int k) {

        Matrix maCenters = new Matrix(k,ma.getColDim());

        int iMaxIterations = 100;
        // Array with indices
        // The index specified the corresponding mean cluster
        int[] arrIndex = new int[ma.getRowDim()];

        // Generate the first mean clusters by random
        Random rnd = new Random();
        for (int ii = 0; ii < k; ii++) {
            int rndIndex = rnd.nextInt(ma.getRowDim());
            maCenters.assignRow(ii, ma.getRow(rndIndex));
        }
        // For test
        // maCenters.assignRow(0, ma.getRow(0));
        // maCenters.assignRow(1, ma.getRow(6));
        // maCenters.assignRow(2, ma.getRow(14));

        maCenters = maCenters.getSorted();
        Matrix maCenters2 = new Matrix(maCenters);
        int counter = 0;
        do{
            maCenters = new Matrix(maCenters2);

            for (int ii = 0; ii < arrIndex.length; ii++)
                arrIndex[ii] = -1;

            // Find the next mean cluster for each object in the matrix.
            // Each col represents a mean cluster, each row represents an
            // object.
            Matrix maDist = getDistTanimotoInvReduced(ma, maCenters);
            for (int ii = 0; ii < maDist.getRowDim(); ii++) {
                int index = maDist.getMinRowIndexRND(ii);
                arrIndex[ii] = index;
            }

            // System.out.println("maDist\n" + maDist);

            // Calculate the new mean clusters.
            double[] arrNumObjects = new double[k];
            maCenters2.set(0.0);
            for (int ii = 0; ii < ma.getRowDim(); ii++) {
                int index = arrIndex[ii];
                maCenters2.add2Row(index, ma, ii);
                arrNumObjects[index]++;
            }
            // boolean bEmptyCenter = false;
            for (int ii = 0; ii < maCenters2.getRowDim(); ii++) {
                if(arrNumObjects[ii] > 0)
                    maCenters2.devideRow(ii, arrNumObjects[ii]);
                else {
                    // maCenters2.setRow(ii, -1);
                    maCenters2.assignRow(ii, maCenters.getRow(ii));
                    // bEmptyCenter = true;
                }
            }
/*
            if(bEmptyCenter) {
                String str = "Break because of empty center.\n";
                System.err.print(str);
                break;
            }
*/
            maCenters2 = maCenters2.getSorted();

            // System.out.println("maCenters2\n" + maCenters2);
            counter++;
            if(counter > iMaxIterations) {
                System.err.print("Max num iterations reached.\n");
                // (new Exception()).printStackTrace();
                break;
            }
        } while(!maCenters.equal(maCenters2));

        System.out.println("Number iterations: " + counter);

        return maCenters;
    }
    
    public static final double getCorrPearson(Matrix A, Matrix B) {
        
        final double [] a = A.toArray();
        
        final double [] aCent = ArrayUtilsCalc.getCentered(a);
        
        final double [] aCentNorm = ArrayUtilsCalc.getNormalized(aCent);
        
        final double [] b = B.toArray();
        
        final double [] bCent = ArrayUtilsCalc.getCentered(b);
        
        final double [] bCentNorm = ArrayUtilsCalc.getNormalized(bCent);
        
        final double val = ArrayUtilsCalc.getCorrPearsonStandardized(aCentNorm,bCentNorm);

        return val;
    }

    public static final double getCorrSpearman(Matrix A, Matrix B) {

        final double [] a = A.toArray();
        final double [] b = B.toArray();

        DoubleArray daA = new DoubleArray(a);
        DoubleArray daB = new DoubleArray(b);

        CorrelationCalculator correlationCalculator = new CorrelationCalculator();

        double c = correlationCalculator.calculateCorrelation(daA, daB, CorrelationCalculator.TYPE_SPEARMAN);

        return c;
    }

    public static double getCorrPearson(Matrix A, int col1, int col2) {
        double val = 0;
        Matrix Anorm = A.getCol(col1);
        Anorm = Anorm.getCenteredMatrix();
        Anorm = Anorm.getNormalizedMatrix();

        Matrix Bnorm = A.getCol(col2);
        Bnorm = Bnorm.getCenteredMatrix();
        Bnorm = Bnorm.getNormalizedMatrix();
        val = MatrixFunctions.getCorrPearsonStandardized(Anorm,Bnorm);
        return val;
    }


    private static double getCorrPearsonStandardized(Matrix A, Matrix B) {
        double val = 0;

        double covXY = getCovarianceCentered(A,B);
        double varA = A.getVarianceCentered();
        double varB = B.getVarianceCentered();
        val = covXY / (varA * varB);
        return val;
    }

    public static double getCovariance(Matrix A, Matrix B) {
        double covXY = 0;
        double Amean = A.getMean();
        double Bmean = B.getMean();
        double sum = 0;
        for (int ii = 0; ii < A.getRowDim(); ii++) {
            for (int jj = 0; jj < A.getColDim(); jj++) {
                sum += (A.get(ii,jj) - Amean) * (B.get(ii,jj) - Bmean);
            }
        }
        covXY = sum / (A.getNumElements() - 1);

        return covXY;
    }

    public static double getCovarianceCentered(Matrix A, Matrix B) {
    	
        double covXY = 0;
        
        double sum = 0;
        
    	final int cols = A.cols();
    	
    	final int rows = A.rows();
        
        for (int i = 0; i < rows; i++) {
        	
        	final double [] a = A.getRow(i);
        	
        	final double [] b = B.getRow(i);
        	
            for (int j = 0; j < cols; j++) {
            	
                sum += (a[j] * b[j]);
            }
        }
        
        covXY = sum / (A.getNumElements() - 1);

        return covXY;
    }
    /**
     * generates a matrix with double values between 0 (inclusive) and 1
     * (exclusive).
     * @param rows rows
     * @param cols columns
     * @return matrix
     */
    public static Matrix getRandomMatrix(int rows, int cols) {
        Matrix ma = new Matrix(rows, cols);

        Random rnd = new Random();

        for (int i = 0; i < ma.getRowDim(); i++) {
            for (int j = 0; j < ma.getColDim(); j++) {
                ma.set(i,j, rnd.nextDouble());
            }
        }

        return ma;
    }


    /**
     *
     * @param n
     * @param mean
     * @param sd
     * @return random one column matrix with given mean and standard deviation.
     */
    public static Matrix rnorm(int n, double mean, double sd) {

        Matrix x = new Matrix(n, 1);

        Random random = new Random();

        for (int i = 0; i < n; i++) {
            double v = random.nextGaussian();

            v = v*sd+mean;

            x.set(i,0,v);
        }
        return x;
    }

    /**
     *
     * @param X
     * @return bit wise.
     */
    public static IntVec getNonZeroCols(Matrix X) {

        int s = X.cols() / Integer.SIZE;

        if(X.cols() % Integer.SIZE>0){
            s++;
        }

        IntVec ivMask = new IntVec(s);

        int r = X.rows();
        int c = X.cols();

        for (int i = 0; i < c; i++) {

            boolean zero=true;
            for (int j = 0; j < r; j++) {

                if(Math.abs(X.get(j,i)) > Matrix.TINY16){
                    zero=false;
                    break;
                }
            }
            if(!zero) {
                ivMask.setBit(i);
            }
        }

        return ivMask;
    }

    /**
     * Converts a vector of vectors into doubles, each vector results in a row in
     * the matrix. All vectors have to be of equal length or a runtime  exception
     * is thrown.
     * @param vecvec vector on vectors, has to be converted into doubles
     * @param ma resulting matrix
     */
    public static void vecvec2Matrix(Vector<Vector<Double>> vecvec, Matrix ma) {

      Iterator<Vector<Double>> it = vecvec.iterator();
      int iLenCol0 = ((Vector<Double>) it.next()).size();

      for( ; it.hasNext(); ) {
        int iLen = ((Vector<Double>) it.next()).size();
        if(iLen != iLenCol0) {
          throw new RuntimeException("All vectors must have the same length.");
        }
      }
      int iRows = vecvec.size();

      ma.resize(iRows, iLenCol0);

      it = vecvec.iterator();
      int iRow = 0;
      for( ; it.hasNext(); ) {
        Vector<Double> vec = new Vector<Double>(((Vector<Double>) it.next()));
        for(int ii = 0; ii < vec.size(); ii++) {
          ma.set(iRow, ii, ((Double) vec.get(ii)).doubleValue());
        }
        iRow++;
      }
    }

    public static Matrix readCSV(File fi) throws IOException {

        List<String[]> li = new ArrayList<String[]>();

        BufferedReader br = new BufferedReader(new FileReader(fi));

        int cols = -1;

        String l = null;
        while ((l=br.readLine())!=null){

            String [] arr = l.split(",");

            if(cols == -1) {
                cols=arr.length;
            } else {
                if(arr.length!=cols){
                    throw new RuntimeException("Number of columns differ!");
                }
            }

            li.add(arr);
        }

        Matrix ma = new Matrix(li.size(), cols);

        for (int i = 0; i < li.size(); i++) {

            String [] arr = li.get(i);

            for (int j = 0; j < arr.length; j++) {

                double v = Double.parseDouble(arr[j]);

                ma.set(i,j,v);

            }
        }

        return ma;
    }

    public static Matrix read(File fiMatrix) throws IOException {
        return read(fiMatrix, false);
    }

    public static Matrix read(File fiMatrix, boolean skipFirstLine) throws IOException {

        FileInputStream fis = new FileInputStream(fiMatrix);

        Matrix A = read(fis, skipFirstLine);

        fis.close();

        return A;
    }

    public static Matrix read(String s) throws IOException {

        InputStream is = new ByteArrayInputStream(s.getBytes());

        Matrix A = read(is, false);

        is.close();

        return A;
    }

    public static Matrix readAsLineBase64Encoded(String s) throws IOException {

        Decoder decoder = Base64.getDecoder();

        byte [] b = decoder.decode(s);

        InputStream is = new ByteArrayInputStream(b);

        Matrix A = read(is, false);

        is.close();

        return A;
    }


    public static Matrix read(InputStream is) {
        return read(is, false);
    }

    public static Matrix read(InputStream is, boolean skipFirstLine) {

        List<DoubleArray> li = new ArrayList<DoubleArray>();

        Scanner scannerLine = new Scanner (is);

        if(skipFirstLine){
            scannerLine.nextLine();
        }

        int rows = 0;
        int cols = 0;
        while(scannerLine.hasNextLine()) {
            ++rows;
            Scanner scannerValue = new Scanner(scannerLine.nextLine());

            DoubleArray arr = null;

            if(cols==0) {
                arr = new DoubleArray();
                while(scannerValue.hasNextDouble()) {
                    double d = scannerValue.nextDouble();
                    arr.add(d);
                    cols++;
                }
            } else {
                arr = new DoubleArray(cols);
                while(scannerValue.hasNextDouble()) {
                    double d = scannerValue.nextDouble();
                    arr.add(d);
                }
            }

            li.add(arr);
        }

        scannerLine.close();

        double[][] a = new double[rows][cols];

        for (int i = 0; i < li.size(); i++) {
            DoubleArray arr = li.get(i);
            for (int j = 0; j < cols; j++) {
                a[i][j] = arr.get(j);
            }
        }

        return new Matrix(a);
    }


    public static void writeQuadraticSymmetricMatrixPairwise(Matrix ma, File fiTxt) throws IOException {

        if(ma.rows() != ma.cols()){
            throw new RuntimeException("Not a quadratic matrix");
        }

        BufferedWriter bw = new BufferedWriter(new FileWriter(fiTxt));

        int n = ma.rows();

        int rows2Write = ((n * n) - n)/2;

        int cc=0;
        for (int i = 0; i < ma.rows(); i++) {

            for (int j = i+1; j < ma.cols(); j++) {

                double v0 = ma.get(i,j);
                double v1 = ma.get(j,i);

                if(v0!=v1){
                    throw new RuntimeException("Matrix not symmetric");
                }

                StringBuilder sb = new StringBuilder();

                sb.append(i);
                sb.append("\t");
                sb.append(j);
                sb.append("\t");
                sb.append(v0);

                if(cc < rows2Write-1){
                    sb.append("\n");
                }
                cc++;

                bw.write(sb.toString());
            }
        }
        bw.close();
    }


    public static void columnIntoDoubleArray(Matrix A, int col, DoubleArray da) {

        da.clear();

        int r = A.rows();

        for (int i = 0; i < r; i++) {
            da.add(A.get(i, col));
        }
    }


    public static List<IdentifiedObject<double []>> createIdentifiedObject(Matrix A) {

        int rows = A.rows();

        List<IdentifiedObject<double []>> liIdentifiedObject = new ArrayList<>(rows);

        for (int i = 0; i < rows; i++) {

            double [] a = A.getRow(i);

            IdentifiedObject<double []> io = new IdentifiedObject<>(a, i);

            liIdentifiedObject.add(io);
        }

        return liIdentifiedObject;
    }
}