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

package com.actelion.research.calc.regression.linear.pls;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.MatrixFunctions;
import com.actelion.research.calc.MatrixTests;


/**
 * <p>Title: SimPLS </p>
 * <p>Description: SimPLS: Simple partial least squares algorithm according:
 * de Jong, S. SIMPLS: An alternative approach to partial least squares
 * regression Chemometrics and Intelligent Laboratory Systems (1993) 18
 * pp 251-263.  </p>
 * @author Modest von Korff
 * @version 1.0
 * 21.11.2003 MvK: Start implementation
 */

public class SimPLS {

    protected Matrix R, T, P, Q, U, V, B, EX, EY, ES;

    public SimPLS() {
    	
        // Weight matrix T = X R
        R = new Matrix();
        
        // X block factor scores
        T = new Matrix();
        
        // Loadings (not orthogonal)
        P = new Matrix();
        
        // Weight vectors
        Q = new Matrix();
        
        // Y block factor scores
        U = new Matrix();
        V = new Matrix();
        B = new Matrix();
        // Explained variance
        EX = new Matrix();
        EY = new Matrix();
        ES = new Matrix();
    }

    /**
     * 
     * @param B
     * @param Xtrain uncentered matrix. Is used to center Xtest.
     * @param Xtest uncentered matrix. Will be centered in method.
     * @param Ytrain uncentered matrix. Is used to de-center Yhat.
     * @return
     */
    public static Matrix invLinReg_Yhat(Matrix B, Matrix Xtrain, Matrix Xtest, Matrix Ytrain) {
    	
		Matrix XTestCentered = new Matrix(Xtest);

		Matrix YdachC = new Matrix();
		
		Matrix Ydach = new Matrix(Xtest.getRowDim(), Ytrain.getColDim());

		Matrix maMeanColXTrain = Xtrain.getMeanCols();

		double [] arrMeanColXTrain = maMeanColXTrain.toArray();

		if (B.getMaximumValue() < Double.MAX_VALUE) {

			final int rowsTest = Xtest.rows();
			
			final int cols = Xtrain.cols();
			
			// Zentrieren des Testdatensatzes
			for (int i = 0; i < rowsTest; i++) {
				
				final double [] arrRow = XTestCentered.getRow(i);
				
				for (int j = 0; j < cols; j++) {
					arrRow[j] -= arrMeanColXTrain[j];
				}
			}

			// Y Dach berechnen
			YdachC = XTestCentered.multiply(B);

			// dezentrieren von Ydach

			for (int i = 0; i < Ytrain.getColDim(); i++) {
				double meanYtrain = Ytrain.getMeanCol(i);
				for (int j = 0; j < Ydach.getRowDim(); j++) {
					double valDecentered = YdachC.get(j, i) + meanYtrain;
					Ydach.set(j, i, valDecentered);
				}
			}
		} else {
			Ydach.resize(Xtest.getRowDim(), Ytrain.getColDim());
			// Die Zahl nicht zu hoch waehlen sonst gibt es Probleme mit
			// overflow
			Ydach.set(Integer.MAX_VALUE);
		}

		return Ydach;
    }

    /**
     * 
     * @param BFromCentered must be calculated from centered X and Y data.
     * @param Xtest
     * @return
     */
    public static Matrix invLinReg_Yhat(Matrix BFromCentered, Matrix Xtest) {
    	
		Matrix Ydach = new Matrix();

		if (BFromCentered.getMaximumValue() < Double.MAX_VALUE) {
			
			// Calculate Y hat.
			Ydach = Xtest.multiply(BFromCentered);
			
		} else {
			Ydach.resize(Xtest.getRowDim(), BFromCentered.cols());
			// If set to large there is maybe a problem with overflow later on.
			Ydach.set(Integer.MAX_VALUE);
		}

		return Ydach;
    }

    public Matrix getExplainedVarS() {
        return ES;
    }
    public Matrix getExplainedVarX() {
        return EX;
    }
    public Matrix getExplainedVarY() {
        return EY;
    }
    public Matrix getP() {
        return P;
    }
    public Matrix getQ() {
        return Q;
    }
    public Matrix getR() {
        return R;
    }
    public Matrix getT() {
        return T;
    }
    public Matrix getU() {
        return U;
    }
    public Matrix getV() {
        return V;
    }

    /**
     * No explained variance is calculated.
     * B = null;
     * EX = null;
     * EY = null;
     * ES = null;
     * 
     * @param Xtrain
     * @param Ytrain
     * @param facmax
     */
    public void simPlsSave(Matrix Xtrain, Matrix Ytrain, int facmax) {
        Matrix Y0 = new Matrix(Ytrain);
        Matrix Qd;

        Matrix S = Xtrain.multiply(true, false, Y0);

        Matrix d = new Matrix();
        Matrix e = new Matrix();
        Matrix q = new Matrix();
        Matrix r, t, p, u, v, vp, vSP;
        Matrix sq;
        Matrix TP, VP;
        
        B = null;
        EX = null;
        EY = null;
        ES = null;

        // Determine number of factors
        int fac, facposs = 0;
        fac = facmax;
        // Determine the maximum number of factors
        if (Xtrain.getRowDim() <= Xtrain.getColDim())
            facposs = Xtrain.getRowDim();
        else if (Xtrain.getRowDim() > Xtrain.getColDim())
            facposs = Xtrain.getColDim();

        if (fac > facposs)
            fac = facposs;

        for (int a = 1; a < fac + 1; a++) {

            // System.out.println("Factor: " + a);

            Qd = S.multiply(true, false, S);
            Matrix.getEigenvector(Qd, Qd.getRowDim(), d, e);
            // Find the largest Eigen value (dominant Eigenvector)
            int index = 0;
            double dMax = d.get(0, 0);
            for (int i = 1; i < d.getRowDim(); i++) {
                if (dMax < d.get(i, 0)) {
                    dMax = d.get(i, 0);
                    index = i;
                }
            }

            q.resize(Qd.getRowDim(), 1);
            for (int i = 0; i < d.getRowDim(); i++)
                q.set(i, 0, Qd.get(i, index));

            /* If there is only one response column (Y), the part above the
             statement cab be replaced by the following statement. Do not
             forget to comment the following statement out if the term above
             is commented out. */
            // q = S.multiply(true, false, S);

            r = S.multiply(false, false, q);
            t = Xtrain.multiply(false, false, r);
            t = t.getCenteredMatrix();
            sq = t.multiply(true, false, t);

            if (Math.sqrt(sq.get(0, 0)) > 10e-50) {
                double normt = Math.sqrt(sq.get(0, 0));
                t = t.devide(normt);
                r = r.devide(normt);
            } else {
                System.err.println(
                    "Division by ZERO error in SimPls(...), normt. Factor " + a + ".");
                // t = t;
                // r = r;
                break;
            }

            p = Xtrain.multiply(true, false, t);
            q = Y0.multiply(true, false, t);
            u = Y0.multiply(false, false, q);
            v = p;
            if (a > 1) {
                VP = V.multiply(true, false, p);
                v = v.subtract(V.multiply(false,false, VP));
                TP = T.multiply(true, false, u);
                u = u.subtract(T.multiply(false,false, TP));
            }

            vp = v.multiply(true, false, v);

            if (vp.get(0, 0) > 10e-50) {
                double normv = Math.sqrt(vp.get(0, 0));
                v = v.devide(normv);
            } else {
                System.err.println(
                    "Division by ZERO error in SimPlsSave(...), normv. Factor " + a + ".");
                // v = v;
                break;
            }
            vSP = v.multiply(true, false, S);
            S = S.subtract(v.multiply(false,false,vSP));

            R.resize(r.getRowDim(), a);
            R.assignCol(a - 1, r);
            T.resize(t.getRowDim(), a);
            T.assignCol(a - 1, t);
            P.resize(p.getRowDim(), a);
            P.assignCol(a - 1, p);
            Q.resize(q.getRowDim(), a);
            Q.assignCol(a - 1, q);
            U.resize(u.getRowDim(), a);
            U.assignCol(a - 1, u);
            V.resize(v.getRowDim(), a);
            V.assignCol(a - 1, v);

        }
    }
    
    /**
     * The explained variance is calculated. This is a time consuming process, 
     * the calculation of the correlation coefficient is takes longer than the PLS calculation itself. 
     * @param Xtrain
     * @param Ytrain
     * @param facmax
     */
    public void simPlsSaveExplainedVariance(Matrix Xtrain, Matrix Ytrain, int facmax) {
        Matrix Y0 = new Matrix(Ytrain);
        Matrix Qd;

        Matrix S = Xtrain.multiply(true, false, Y0);
        Matrix Sorig = new Matrix(S);

        Matrix d = new Matrix();
        Matrix e = new Matrix();
        Matrix q = new Matrix();
        Matrix r, t, p, u, v, vp, vSP;
        Matrix sq;
        Matrix TP, VP;
        Matrix h, varX, varY;
        Matrix qt, tmp;

        // System.out.println("Y0\n" + Y0.toString(3));

        // System.out.println("SimPLS\nS: " + S.getRowDim() + " " + S.getColDim());

        // Determine number of factors
        int fac, facposs = 0;
        fac = facmax;
        // Determine the maximum number of factors
        if (Xtrain.getRowDim() <= Xtrain.getColDim())
            facposs = Xtrain.getRowDim();
        else if (Xtrain.getRowDim() > Xtrain.getColDim())
            facposs = Xtrain.getColDim();

        if (fac > facposs)
            fac = facposs;

        for (int a = 1; a < fac + 1; a++) {

            // System.out.println("Factor: " + a);

            Qd = S.multiply(true, false, S);
            Matrix.getEigenvector(Qd, Qd.getRowDim(), d, e);
            // Find the largest eigen value (dominant eigenvector)
            int index = 0;
            double dMax = d.get(0, 0);
            for (int i = 1; i < d.getRowDim(); i++) {
                if (dMax < d.get(i, 0)) {
                    dMax = d.get(i, 0);
                    index = i;
                }
            }

            q.resize(Qd.getRowDim(), 1);
            for (int i = 0; i < d.getRowDim(); i++)
                q.set(i, 0, Qd.get(i, index));

            /* If there is only one response column (Y), the part above the
             statement cab be replaced by the following statement. Do not
             forget to comment the following statement out if the term above
             is commented out. */
            // q = S.multiply(true, false, S);

            r = S.multiply(false, false, q);
            t = Xtrain.multiply(false, false, r);
            t = t.getCenteredMatrix();
            sq = t.multiply(true, false, t);

            if (Math.sqrt(sq.get(0, 0)) > 10e-50) {
                double normt = Math.sqrt(sq.get(0, 0));
                t = t.devide(normt);
                r = r.devide(normt);
            }
            else {
                System.err.println(
                    "Division by ZERO error in SimPlsSave(...), normt. Factor " + a + ".");
                // t = t;
                // r = r;
                break;
            }

            p = Xtrain.multiply(true, false, t);
            q = Y0.multiply(true, false, t);
            u = Y0.multiply(false, false, q);
            v = p;
            if (a > 1) {
                VP = V.multiply(true, false, p);
                v = v.subtract(V.multiply(false,false, VP));
                TP = T.multiply(true, false, u);
                u = u.subtract(T.multiply(false,false, TP));
            }


            vp = v.multiply(true, false, v);

            if (vp.get(0, 0) > 10e-50) {
                double normv = Math.sqrt(vp.get(0, 0));
                v = v.devide(normv);
            }
            else {
                System.err.println(
                    "Division by ZERO error in SimPlsSave(...), normv. Factor " + a + ".");
                // v = v;
                break;
            }
            vSP = v.multiply(true, false, S);
            S = S.subtract(v.multiply(false,false,vSP));

            R.resize(r.getRowDim(), a);
            R.assignCol(a - 1, r);
            T.resize(t.getRowDim(), a);
            T.assignCol(a - 1, t);
            P.resize(p.getRowDim(), a);
            P.assignCol(a - 1, p);
            Q.resize(q.getRowDim(), a);
            Q.assignCol(a - 1, q);
            U.resize(u.getRowDim(), a);
            U.assignCol(a - 1, u);
            V.resize(v.getRowDim(), a);
            V.assignCol(a - 1, v);

            // Explained variance in X
            Matrix Xhat = T.multiply(false,true, P);
            double corrEX = MatrixFunctions.getCorrPearson(Xhat, Xtrain);
            EX.resize(1, a);
            EX.set(0, a-1, corrEX);

            tmp = r.multiply(false, true, q);
            B.resize(tmp.getRowDim(), tmp.getColDim());
            B = B.plus(tmp);
            Matrix Yhat = Xtrain.multiply(B);
            EY.resize(1, a);
            double corrEY = MatrixFunctions.getCorrPearson(Yhat, Y0);
            EY.set(0, a-1, corrEY);

            Matrix Shat = Xhat.multiply(true, false, Yhat);
            double corrES = MatrixFunctions.getCorrPearson(Shat, Sorig);
            ES.resize(1, a);
            ES.set(0, a-1, corrES);

        }
    }

    /**
     * Explained variance in X
     * @param A
     * @param Ahat
     * @return
     */
    private static double explainedVariance(Matrix A, Matrix Ahat) {

        double mean = 0;
        Matrix percent = null;

        // System.out.println("Xhat\n" + Xhat.toString(3));

        // percent = A.subtract(Ahat);
        percent = A.devideDivisorBigger(Ahat);
        System.out.println("percent\n" + percent.toString(3));

        // percent = percent.getAbs().subtract(100);
        // Matrix meanCols = percent.getMeanCols();
        mean = Math.abs(percent.getMean());

        return mean;
    }

    public void simPlsSave(Matrix Xtrain, Matrix Ytrain, Matrix Xtest, Matrix Ytest, int facmax) {

        Matrix Y0 = new Matrix(Ytrain);
        Matrix Qd;
        Matrix S = Xtrain.multiply(true, false, Y0);

        Matrix d = new Matrix();
        Matrix e = new Matrix();
        Matrix q = new Matrix();

        
        Matrix r, t, p, u, v, vp, vSP;
        Matrix sq;
        Matrix TP, VP;

        Matrix h, varX, varY;
        Matrix qt, tmp;

        // Determine number of factors
        int fac, facposs = 0;
        fac = facmax;
        // Determine the maximum number of factors
        if (Xtrain.getRowDim() <= Xtrain.getColDim())
            facposs = Xtrain.getRowDim();
        else if (Xtrain.getRowDim() > Xtrain.getColDim())
            facposs = Xtrain.getColDim();

        if (fac > facposs)
            fac = facposs;

        for (int a = 1; a < fac + 1; a++) {

            Qd = S.multiply(true, false, S);
            Matrix.getEigenvector(Qd, Qd.getRowDim(), d, e);
            // Find the largest eigen value (dominant eigenvector)
            int index = 0;
            double dMax = d.get(0, 0);
            for (int i = 1; i < d.getRowDim(); i++) {
                if (dMax < d.get(i, 0)) {
                    dMax = d.get(i, 0);
                    index = i;
                }
            }
            q.resize(Qd.getRowDim(), 1);
            for (int i = 0; i < d.getRowDim(); i++)
                q.set(i, 0, Qd.get(i, index));

            /* If there is only one response column (Y), the part above the
             statement cab be replaced by the following statement. Do not
             forget to comment the following statement out if the term above
             is commented out. */
            // q = S.multiply(true, false, S);

            r = S.multiply(false, false, q);
            t = Xtrain.multiply(false, false, r);
            t = t.getCenteredMatrix();
            sq = t.multiply(true, false, t);

            if (Math.sqrt(sq.get(0, 0)) > 10e-50) {
                double normt = Math.sqrt(sq.get(0, 0));
                t = t.devide(normt);
                r = r.devide(normt);
            }
            else {
                System.err.println(
                    "Division by ZERO error in SimPlsSave(...), normt.");
                // t = t;
                // r = r;
                B.set(Float.MAX_VALUE);
                break;
            }

            p = Xtrain.multiply(true, false, t);
            q = Y0.multiply(true, false, t);
            u = Y0.multiply(false, false, q);
            v = p;
            if (a > 1) {
                VP = V.multiply(true, false, p);
                v = v.subtract(V.multiply(false,false, VP));
                TP = T.multiply(true, false, u);
                u = u.subtract(T.multiply(false,false, TP));
            }
            vp = v.multiply(true, false, v);

            if (vp.get(0, 0) > 10e-50) {
                double normv = Math.sqrt(vp.get(0, 0));
                v = v.devide(normv);
            }
            else {
                System.err.println(
                    "Division by ZERO error in SimPlsSave(...), normv.");
                B.set(Float.MAX_VALUE);
                // v = v;
                break;
            }
            vSP = v.multiply(true, false, S);
            S = S.subtract(v.multiply(false,false,vSP));

            R.resize(r.getRowDim(), a);
            R.assignCol(a - 1, r);
            T.resize(t.getRowDim(), a);
            T.assignCol(a - 1, t);
            P.resize(p.getRowDim(), a);
            P.assignCol(a - 1, p);
            Q.resize(q.getRowDim(), a);
            Q.assignCol(a - 1, q);
            U.resize(u.getRowDim(), a);
            U.assignCol(a - 1, u);
            V.resize(v.getRowDim(), a);
            V.assignCol(a - 1, v);

            tmp = r.multiply(false, true, q);
            B.resize(tmp.getRowDim(), tmp.getColDim());
            B = B.plus(tmp);

            // Matrix Yhat = invLinReg_Yhat(B, Xtrain, Xtest, Ytrain);
            // System.out.println("B:\n" + B + "\n");
            // System.out.println("Yhat:\n" + Yhat.toString(5) + "\n");
        }
    }

    /**
     *
     * @param XPreprocessed has to be preprovcessed with XTrain. I.e. cebntered by the mean values of XTrain.
     * @return
     */
    public Matrix getT(Matrix XPreprocessed){

        Matrix Tcalc = XPreprocessed.multiply(getR());

        return Tcalc;

    }


    public static void main(String[] args) {
       // Matrix X = MatrixFunctions.testMatrix_XWine();
       // Matrix Y = MatrixFunctions.testMatrix_YWine();
       // Matrix X = MatrixFunctions.testLonglyX();
       // Matrix Y = MatrixFunctions.testLonglyY();
       // Matrix X = MatrixFunctions.test08();
       // Matrix Y = MatrixFunctions.test06();

       Matrix X = MatrixTests.testDescriptor01X();
       Matrix Y = MatrixTests.testDescriptor01Y();

       // X = X.log();

       Matrix Xc = X.getCenteredMatrix();
       Matrix Yc = Y.getCenteredMatrix();

       SimPLS simPLS = new SimPLS();
       int iFactors = 2;
       // simPLS.simPlsSave(Xs, Ys, Xs, Ys, iFactors);
       simPLS.simPlsSave(Xc, Yc, iFactors);

       System.out.println(simPLS.toString(4));

       Matrix P = simPLS.getP();
       Matrix R = simPLS.getR();
       Matrix U = simPLS.getU();
       Matrix V = simPLS.getV();
       Matrix Q = simPLS.getQ();
       Matrix T = simPLS.getT();

       Matrix That = Xc.multiply(R);
       System.out.println("That\n" + That.toString(4));

       Matrix G = R.multiply(false,true,Q);
       System.out.println("U\n" + U.toString(4));
       System.out.println("G\n" + G.toString(4));

       Matrix Xhat = T.multiply(false,true, P);
       System.out.println(Xhat.toString(4));

       Matrix XG = X.multiply(false,false,G);
       System.out.println("XG\n" + XG.toString(4));

    }
    public String toString(int iDigits) {

        String sOut = "";
        sOut += "E [%]\n" + EX.toString(iDigits) + "\n";
        sOut += "P:\n" + P.toString(iDigits) + "\n";
        sOut += "Q:\n" + Q.toString(iDigits) + "\n";
        sOut += "R:\n" + R.toString(iDigits) + "\n";
        sOut += "T:\n" + T.toString(iDigits) + "\n";
        sOut += "U:\n" + U.toString(iDigits) + "\n";
        sOut += "V:\n" + V.toString(iDigits) + "\n";

        return sOut;
    }
}