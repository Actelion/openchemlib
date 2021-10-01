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
package com.actelion.research.chem.optimization;

import java.util.Arrays;




/**
 *
 * taken from DD_chem3d, small changes necessary because a different output structure is needed
 * it returns not only a value (the objective function), but the transformation array for achieving the best alignment
 */
public class OptimizerLBFGS   {

	int maxIterations;	
	double minRMS;
	

	
	public OptimizerLBFGS(int maxIterations, double minRMS){
		this.maxIterations = maxIterations;
		this.minRMS = minRMS;
		
	}

	

	/**
	 * Optimization routine using the limited Broyden-Fletcher-Goldfarb-Shanno
	 * algorithm
	 */

	public synchronized double[] optimize(Evaluable eval) {
	//public synchronized double optimize(AbstractEvaluable eval) {
		//System.out.println("new opti");
		double[] initial = eval.getState();
		int N = initial.length;

		final int MSAV = N; // Math.max(1, Math.min(12, N));
		double[] alpha = new double[MSAV];
		double[] rho = new double[MSAV];
		double gamma = 1;
		int m = 0;
		int nErrors = 0;

		// evaluate the function and get the initial gradient
		double[] grad = new double[N];
		double f = eval.getFGValue(grad);
		//System.out.println(f);
		double fOld = f;
		double f0 = f;
		double gNorm;
		double fMove = 0;
		
		if (N == 0) {
			return eval.getState();
		}

			//eturn f;
		double[] oldX = initial;
		double[] oldGradient = new double[N];
		double[][] s = new double[MSAV][N];
		double[][] y = new double[MSAV][N];
		double[] h0 = new double[N];
		double[] q = new double[N];
		double[] r = new double[N];

		boolean restart = true;
		int mUse = 0;
		int iteration ;
		for (iteration = 1; iteration <= maxIterations; iteration++) {
			if (restart) {

				mUse = 0;
				f = eval.getFGValue(grad);
				//System.out.println(f);
				gamma = 1;
				fMove = .25 * getNorm(grad);
				restart = false;
			}

			gNorm = getNorm(grad);
			double RMS = gNorm / Math.sqrt(N);

			if (RMS < minRMS) {
				break;

			} else if (nErrors > 2) {
				break;
			}
			// Estimate Hessian diagonal
			m = (m + 1) % MSAV;
			Arrays.fill(h0, gamma);
			System.arraycopy(grad, 0, q, 0, N);

			int k = m;
			for (int j = 0; j < mUse; j++) {
				k = k == 0 ? MSAV - 1 : k - 1;
				alpha[k] = 0;
				for (int i = 0; i < N; i++)
					alpha[k] += s[k][i] * q[i];
				alpha[k] *= rho[k];
				for (int i = 0; i < N; i++)
					q[i] -= y[k][i] * alpha[k];
			}
			for (int i = 0; i < N; i++)
				r[i] = h0[i] * q[i];
			for (int j = 0; j < mUse; j++) {
				double beta = 0;
				for (int i = 0; i < N; i++)
					beta += y[k][i] * r[i];
				beta *= rho[k];
				for (int i = 0; i < N; i++)
					r[i] += s[k][i] * (alpha[k] - beta);
				k = (k + 1) % MSAV;
			}

			// set search direction and store current point and gradient
			for (int i = 0; i < N; i++) {
				r[i] = -r[i];
			}

			// Memorize position
			oldX = eval.getState();
			System.arraycopy(grad, 0, oldGradient, 0, N);

			// perform line search along the new conjugate direction
			Object[] res = Lnsrch.minimizeEnergyAroundDirection(eval, f, grad, r, fMove);
			f = (Double) res[0];
			grad = (double[]) res[1];
			if (res[2] == Boolean.FALSE) {
				nErrors++;
				restart = true;
			}

			// Update variables
			double ys = 0, yy = 0;
			double[] newState = eval.getState();
			for (int i = 0; i < N; i++) {
				s[m][i] = newState[i] - oldX[i];
				y[m][i] = grad[i] - oldGradient[i];
				ys += y[m][i] * s[m][i];
				yy += y[m][i] * y[m][i];
			}
			gamma = Math.abs(ys / yy);
			if (ys == 0) {
				restart = true;
				continue;
			}

			rho[m] = 1.0 / ys;
			fMove = fOld - f;
			fOld = f;
			mUse = Math.min(mUse + 1, MSAV);

		}

		if (f > f0) {
			eval.setState(initial);
			f = f0;
		}
		return eval.getState();
	}

	public final static double getRMS(double[] vector) {return Math.sqrt(getNormSq(vector) / vector.length);}	
	
	public static double getNorm(double[] vector) {
		 return Math.sqrt(getNormSq(vector));

}
	
	public final static double getNormSq(double[] vector) {
		double res = 0;
		for (int i = 0; i < vector.length; i++) res += vector[i] * vector[i];
		return res;
	}


}
