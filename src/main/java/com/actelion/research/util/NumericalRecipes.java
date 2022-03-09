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

/*
 * @(#)NumericalRecipes.java
 *
 * contains routines from 'Numerical Recipes in C' ported from C to Java
 *
 * @original authors W.H.Press, B.P.Fannery, S.A.Teukolsky, W.T.Vettering
 */

package com.actelion.research.util;

public class NumericalRecipes {
	private double _fit_siga,_fit_sigb,_fit_a,_fit_b,_fit_chi2,_fit_q;	// by reference params of fit()
	private double _gser_gamser, _gser_gln;								// by reference params of gser()
	private double _gcf_gammcf, _gcf_gln;								// by reference params of gcf()

	private double ochisq;			// used by mrqmin
	private double[] atry,beta,da;
	private double[][] oneda;
	private int mfit;

	public static void svbksb(double[][] u, double[] w, double[][] v, int m, int n, double[] b, double[] x) {
//	Solves A·X = B for a vector X, where A is specified by the arrays u[1..m][1..n], w[1..n],
//	v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for
//	square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
//	No input quantities are destroyed, so the routine may be called sequentially with different b's.
		double[] tmp = new double[n];
		for (int j=0; j<n; j++) {	// Calculate UTB.
			double s=0.0;
			if (w[j] != 0.0) {	// Nonzero result only if wj is nonzero.
				for (int i=0; i<m; i++)
					s += u[i][j]*b[i];
				s /= w[j];	// This is the divide by wj .
				}
			tmp[j]=s;
			}
		for (int j=0; j<n; j++) {	// Matrix multiply by V to get answer.
			double s=0.0;
			for (int jj=0; jj<n; jj++)
				s += v[j][jj]*tmp[jj];
			x[j]=s;
			}
		}


	public static void svdcmp(double[][] a, int m, int n, double[] w, double[][] v) throws Exception {
//	Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A =
//	U·W·V T. Thematrix U replaces a on output. The diagonal matrix of singular values W is output
//	as a vector w[1..n]. Thematrix V (not the transpose V T ) is output as v[1..n][1..n].
		double anorm,c,f,g,h,s,scale,x,y,z;

		double[] rv1 = new double[n];
		g=scale=0.0;	// Householder reduction to bidiagonal form.
		anorm=0.0;
		int l=0;
		for (int i=0; i<n; i++) {
			l=i+1;
			rv1[i]=scale*g;
			g=s=scale=0.0;
			if (i < m) {
				for (int k=i; k<m; k++)
					scale += Math.abs(a[k][i]);
				if (scale != 0.0) {
					for (int k=i; k<m; k++) {
						a[k][i] /= scale;
						s += a[k][i]*a[k][i];
						}
					f=a[i][i];
					g = -SIGN(Math.sqrt(s),f);
					h=f*g-s;
					a[i][i]=f-g;
					for (int j=l; j<n; j++) {
						s=0.0;
						for (int k=i; k<m; k++)
							s += a[k][i]*a[k][j];
						f=s/h;
						for (int k=i; k<m; k++)
							a[k][j] += f*a[k][i];
						}
					for (int k=i; k<m; k++)
						a[k][i] *= scale;
					}
				}
			w[i]=scale *g;
			g=s=scale=0.0;
			if (i < m && i != n-1) {
				for (int k=l; k<n; k++)
					scale += Math.abs(a[i][k]);
				if (scale != 0.0) {
					for (int k=l; k<n; k++) {
						a[i][k] /= scale;
						s += a[i][k]*a[i][k];
						}
					f=a[i][l];
					g = -SIGN(Math.sqrt(s),f);
					h=f*g-s;
					a[i][l]=f-g;
					for (int k=l; k<n; k++)
						rv1[k]=a[i][k]/h;
					for (int j=l; j<m; j++) {
						s=0.0;
						for (int k=l; k<n; k++)
							s += a[j][k]*a[i][k];
						for (int k=l; k<n; k++)
							a[j][k] += s*rv1[k];
						}
					for (int k=l; k<n; k++)
						a[i][k] *= scale;
					}
				}
			anorm=Math.max(anorm,(Math.abs(w[i])+Math.abs(rv1[i])));
			}
		for (int i=n-1; i>=0; i--) {	// Accumulation of right-hand transformations.
			if (i < n-1) {
				if (g != 0.0) {
					for (int j=l; j<n; j++)	// Double division to avoid possible underflow.
						v[j][i]=(a[i][j]/a[i][l])/g;
					for (int j=l; j<n; j++) {
						s=0.0;
						for (int k=l; k<n; k++)
							s += a[i][k]*v[k][j];
						for (int k=l; k<n; k++)
							v[k][j] += s*v[k][i];
						}
					}
				for (int j=l; j<n; j++)
					v[i][j]=v[j][i]=0.0;
				}
			v[i][i]=1.0;
			g=rv1[i];
			l=i;
			}
		for (int i=Math.min(m,n)-1; i>=0; i--) {	// Accumulation of left-hand transformations.
			l=i+1;
			g=w[i];
			for (int j=l; j<n; j++)
				a[i][j]=0.0;
			if (g != 0.0) {
				g=1.0/g;
				for (int j=l; j<n; j++) {
					s=0.0;
					for (int k=l; k<m; k++)
						s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					for (int k=i; k<m; k++)
						a[k][j] += f*a[k][i];
					}
				for (int j=i; j<m; j++)
					a[j][i] *= g;
				}
			else
				for (int j=i; j<m; j++)
					a[j][i]=0.0;
			++a[i][i];
			}
		for (int k=n-1; k>=0; k--) {	// Diagonalization of the bidiagonal form: Loop over
										// singular values, and over allowed iterations.
			for (int its=1;its<=30;its++) {
				boolean flag = true;
				int nm=0;
				for (l=k; l>=0; l--) {	// Test for splitting.
					nm=l-1;				// Note that rv1[1] is always zero.
					if ((double)(Math.abs(rv1[l])+anorm) == anorm) {
						flag = false;
						break;
						}
					if ((double)(Math.abs(w[nm])+anorm) == anorm)
						break;
					}
				if (flag) {
					c=0.0;	// Cancellation of rv1[l], if l > 1.
					s=1.0;
					for (int i=l; i<k; i++) {
						f=s*rv1[i];
						rv1[i]=c*rv1[i];
						if ((double)(Math.abs(f)+anorm) == anorm)
							break;
						g=w[i];
						h=pythag(f,g);
						w[i]=h;
						h=1.0/h;
						c=g*h;
						s = -f*h;
						for (int j=0; j<m; j++) {
							y=a[j][nm];
							z=a[j][i];
							a[j][nm]=y*c+z*s;
							a[j][i]=z*c-y*s;
							}
						}
					}
				z=w[k];
				if (l == k) {	// Convergence.
					if (z < 0.0) {	// Singular value is made nonnegative.
						w[k] = -z;
						for (int j=0; j<n; j++)
							v[j][k] = -v[j][k];
						}
					break;
					}
				if (its == 30)
					nrerror("no convergence in 30 svdcmp iterations");
				x=w[l];	// Shift from bottom 2-by-2 minor.
				nm=k-1;
				y=w[nm];
				g=rv1[nm];
				h=rv1[k];
				f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
				g=pythag(f,1.0);
				f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
				c=s=1.0;	// Next QR transformation:
				for (int j=l; j<=nm; j++) {
					int i=j+1;
					g=rv1[i];
					y=w[i];
					h=s*g;
					g=c*g;
					z=pythag(f,h);
					rv1[j]=z;
					c=f/z;
					s=h/z;
					f=x*c+g*s;
					g = g*c-x*s;
					h=y*s;
					y *= c;
					for (int jj=0; jj<n; jj++) {
						x=v[jj][j];
						z=v[jj][i];
						v[jj][j]=x*c+z*s;
						v[jj][i]=z*c-x*s;
						}
					z=pythag(f,h);
					w[j]=z;	// Rotation can be arbitrary if z = 0.
					if (z != 0.0) {
						z=1.0/z;
						c=f*z;
						s=h*z;
						}
					f=c*g+s*y;
					x=c*y-s*g;
					for (int jj=0; jj<m; jj++) {
						y=a[jj][j];
						z=a[jj][i];
						a[jj][j]=y*c+z*s;
						a[jj][i]=z*c-y*s;
						}
					}
				rv1[l]=0.0;
				rv1[k]=f;
				w[k]=x;
				}
			}
		}


	public void fit(double x[], double y[], int ndata, double sig[], int mwt) throws Exception {
//	Given a set of data points x[1..ndata],y[1..ndata] with individual standard deviations
//	sig[1..ndata], fit them to a straight line y = a + bx by minimizing ÷2. Returned are
//	a,b and their respective probable uncertainties siga and sigb, the chi-square chi2, and the
//	goodness-of-fit probability q (that the fit would have ÷2 this large or larger). If mwt=0 on
//	input, then the standard deviations are assumed to be unavailable: q is returned as 1.0 and
//	the normalization of chi2 is to unit standard deviation on all points.

		double sx = 0.0;
		double sy = 0.0;
		double st2 = 0.0;
		double ss = 0.0;

		_fit_b = 0.0;
		if (mwt != 0) {	// Accumulate sums ...
			ss=0.0;
			for (int i=0; i<ndata; i++) {	// ...with weights
				double wt = 1.0/SQR(sig[i]);
				ss += wt;
				sx += x[i]*wt;
				sy += y[i]*wt;
				}
			}
		else {
			for (int i=0; i<ndata; i++) {	// ...or without weights.
				sx += x[i];
				sy += y[i];
				}
			ss=ndata;
			}

		double sxoss=sx/ss;

		if (mwt == 0) {
			for (int i=0; i<ndata; i++) {
				double t=(x[i]-sxoss)/sig[i];
				st2 += t*t;
				_fit_b += t*y[i]/sig[i];
				}
			}
		else {
			for (int i=0; i<ndata; i++) {
				double t = x[i]-sxoss;
				st2 += t*t;
				_fit_b += t*y[i];
				}
			}

		_fit_b /= st2;	// Solve for a, b, óa, and ób.
		_fit_a = (sy-sx*(_fit_b))/ss;
		_fit_siga = Math.sqrt((1.0+sx*sx/(ss*st2))/ss);
		_fit_sigb = Math.sqrt(1.0/st2);
		_fit_chi2 = 0.0;	// Calculate ÷2.
		_fit_q = 1.0;
		if (mwt == 0) {
			for (int i=0; i<ndata; i++)
				_fit_chi2 += SQR(y[i]-(_fit_a)-(_fit_b)*x[i]);
			double sigdat = Math.sqrt((_fit_chi2)/(ndata-2));
				// For unweighted data evaluate typical
				// sig using chi2, and adjust
				// the standard deviations.
			_fit_siga *= sigdat;
			_fit_sigb *= sigdat;
			}
		else {
			for (int i=0; i<ndata; i++)
				_fit_chi2 += SQR((y[i]-(_fit_a)-(_fit_b)*x[i])/sig[i]);
			if (ndata>2)
				_fit_q = gammq(0.5*(ndata-2),0.5*(_fit_chi2));	// Equation (15.2.12).
			}
		}


	public void mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
				int ma, double[][] covar, double[][] alpha, FittingFunction funcs) throws Exception {
//	Levenberg-Marquardt method, attempting to reduce the value x^2 of a fit between a set of data
//	points x[1..ndata], y[1..ndata] with individual standard deviations sig[1..ndata],
//	and a nonlinear function dependent on ma coeficients a[1..ma]. The input array ia[1..ma]
//	indicates by nonzero entries those components of a that should be fitted for, and by zero
//	entries those components that should be held fixed at their input values. The program returns
//	current best-fit values for the parameters a[1..ma], and x^2 = chisq. The arrays
//	covar[1..ma][1..ma], alpha[1..ma][1..ma] are used as working space during most
//	iterations. Supply a routine funcs(x,a,yfit,dyda,ma) that evaluates the fitting function
//	yfit, and its derivatives dyda[1..ma] with respect to the fitting parameters a at x. On
//	the first call provide an initial guess for the parameters a, and set alamda<0 for initialization
//	(which then sets alamda=.001). If a step succeeds chisq becomes smaller and alamda decreases
//	by a factor of 10. If a step fails alamda grows by a factor of 10. You must call this
//	routine repeatedly until convergence is achieved. Then, make one final call with alamda=0, so
//	that covar[1..ma][1..ma] returns the covariance matrix, and alpha the curvature matrix.
//	(Parameters held fixed will return zero covariances.)
		if (funcs.alamda < 0.0) {	/* Initialization. */
			atry = new double[ma];
			beta = new double[ma];
			da = new double[ma];
			mfit = 0;
			for (int j=0;j<ma;j++)
				if (ia[j] != 0)
					mfit++;
			oneda = new double[mfit][1];
			funcs.alamda=0.001;
			mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,funcs);
			ochisq=(funcs.chisq);
			for (int j=0;j<ma;j++)
				atry[j]=a[j];
			}

		for (int j=0;j<mfit;j++) {	/* Alter linearized fitting matrix, by augmenting diagonal elements. */
			for (int k=0;k<mfit;k++)
				covar[j][k]=alpha[j][k];
			covar[j][j]=alpha[j][j]*(1.0+(funcs.alamda));
			oneda[j][0]=beta[j];
			}

		gaussj(covar,mfit,oneda,1);	/* Matrix solution. */
		for (int j=0;j<mfit;j++)
			da[j]=oneda[j][0];

		if (funcs.alamda == 0.0) {	/* Once converged, evaluate covariance matrix. */
			covsrt(covar,ma,ia,mfit);
			oneda = null;
			da = null;
			beta = null;
			atry = null;
			return;
			}

		int jj=-1;
		for (int l=0;l<ma;l++)	/* Did the trial succeed? */
			if (ia[l] != 0)
				atry[l]=a[l]+da[++jj];

		mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,funcs);
		if (funcs.chisq < ochisq) {	/* Success, accept the new solution. */
			funcs.alamda *= 0.1;
			ochisq=(funcs.chisq);
			for (int j=0;j<mfit;j++) {
				for (int k=0;k<mfit;k++)
					alpha[j][k]=covar[j][k];
				beta[j]=da[j];
				}

			for (int l=0;l<ma;l++)
				a[l]=atry[l];
			}
		else {	/* Failure, increase alamda and return. */
			funcs.alamda *= 10.0;
			funcs.chisq = ochisq;
			}
		}


	private void gser(double a, double x) throws Exception {
//	Returns the incomplete gamma function P(a, x) evaluated
//	by its series representation as gamser.
//	Also returns ln Ã(a) as gln.
		final int ITMAX = 100;
		final double EPS = 3.0e-7;

		_gser_gln = gammln(a);
		if (x <= 0.0) {
			if (x < 0.0)
				nrerror("x less than 0 in routine gser");
			_gser_gamser = 0.0;
			return;
			}
		else {
			double ap = a;
			double del = 1.0/a;
			double sum = del;
			for (int n=1; n<=ITMAX; n++) {
				ap += 1.0;
				del *= x/ap;
				sum += del;
				if (Math.abs(del) < Math.abs(sum)*EPS) {
					_gser_gamser = sum * Math.exp(-x+a*Math.log(x)-(_gser_gln));
					return;
					}
				}
			nrerror("a too large, ITMAX too small in routine gser");
			return;
			}
		}


	private void gcf(double a, double x) throws Exception {
//	original params: (float *gammcf, float a, float x, float *gln)
// Returns the incomplete gamma function Q(a, x) evaluated by its continued
// fraction representation as gammcf. Also returns lnÃ(a) as gln.

		final int ITMAX = 100;			// Maximum allowed number of iterations.
		final double EPS = 3.0e-7;		// Relative accuracy.
		final double FPMIN = 1.0e-30;	// Number near the smallest representable
										// floating-point number.
		_gcf_gln = gammln(a);
		double b = x+1.0-a;	//	Set up for evaluating continued fraction
							//	by modified Lentz’s method (§5.2)
							//	with b0 = 0.
		double c = 1.0 / FPMIN;
		double d = 1.0 / b;
		double h = d;
		int i;
		for (i=1; i<=ITMAX; i++) {	// Iterate to convergence.
			double an = -i*(i-a);
			b += 2.0;
			d = an*d+b;
			if (Math.abs(d) < FPMIN)
				d=FPMIN;
			c = b+an/c;
			if (Math.abs(c) < FPMIN)
				c=FPMIN;
			d = 1.0/d;
			double del = d*c;
			h *= del;
			if (Math.abs(del-1.0) < EPS)
				break;
			}
		if (i > ITMAX)
			nrerror("a too large, ITMAX too small in gcf");
		_gcf_gammcf = Math.exp(-x+a*Math.log(x)-(_gcf_gln))*h; // Put factors in front.
		}


	private double gammq(double a, double x) throws Exception {
//	Returns the incomplete gamma function Q(a, x) \uFFFDß 1 \uFFFD| P(a, x).
		if (x < 0.0 || a <= 0.0)
			nrerror("Invalid arguments in routine gammq");
		if (x < (a+1.0)) {	// Use the series representation
			gser(a, x);
			return 1.0 - _gser_gamser;	// and take its complement.
			}
		else {	// Use the continued fraction representation.
			gcf(a, x);
			return _gcf_gammcf;
			}
		}


	private double gammln(double xx) {
//	Returns the value ln[Ã(xx)] for xx > 0.

//	Internal arithmetic will be done in double precision,
//	a nicety that you can omit if .ve-.gure
//	accuracy is good enough.
		final double[] cof = {	 76.18009172947146,
								-86.50532032941677,
								 24.01409824083091,
								 -1.231739572450155,
								  0.1208650973866179e-2,
								 -0.5395239384953e-5	};
		double x = xx;
		double y = xx;
		double tmp = x + 5.5;
		tmp -= (x+0.5) * Math.log(tmp);
		double ser = 1.000000000190015;
		for (int j=0; j<=5; j++)
			ser += cof[j] / ++y;
		return -tmp + Math.log(2.5066282746310005*ser/x);
		}


	private static void covsrt(double covar[][], int ma, int ia[], int mfit) {
//	Expand in storage the covariance matrix covar, so as to take into account parameters
//	that are being held fixed. (For the latter, return zero covariances.)
		for (int i=mfit; i<ma; i++)
			for (int j=0; j<=i; j++)
				covar[i][j] = covar[j][i] = 0.0;

		int k=mfit-1;
		for (int j=ma-1;j>=0;j--) {
			if (ia[j] != 0) {
				for (int i=0;i<ma;i++) SWAP(covar[i][k],covar[i][j]);
				for (int i=0;i<ma;i++) SWAP(covar[k][i],covar[j][i]);
				k--;
				}
			}
		}


	private static void mrqcof(double x[], double y[], double sig[], int ndata, double a[], int ia[],
					int ma, double alpha[][], double beta[], FittingFunction funcs) {
//	Used by mrqmin to evaluate the linearized fitting matrix alpha, and vector beta as in (15.5.8),
//	and calculate X2
		double[] dyda = new double[ma];

		int mfit = 0;
		for (int j=0; j<ma; j++)
			if (ia[j] != 0)
				mfit++;

		for (int j=0;j<mfit;j++) {	/* Initialize (symmetric) alpha, beta.	*/
			for (int k=0;k<=j;k++)
				alpha[j][k]=0.0;
			beta[j]=0.0;
			}

		funcs.chisq=0.0;
		for (int i=0;i<ndata;i++) { /* Summation loop over all data. */
			double ymod = funcs.fittingFunction(x[i], a, dyda, ma);
			double sig2i=1.0/(sig[i]*sig[i]);
			double dy=y[i]-ymod;
			int j = -1;
			for (int l=0;l<ma;l++) {
				if (ia[l] != 0) {
					double wt=dyda[l]*sig2i;
					j++;
					int k = -1;
					for (int m=0;m<=l;m++)
						if (ia[m] != 0) alpha[j][++k] += wt*dyda[m];
					beta[j] += dy*wt;
					}
				}
			funcs.chisq += dy*dy*sig2i;
			}

		for (int j=1;j<mfit;j++)	/* Fill in the symmetric side. */
			for (int k=0;k<j;k++)
				alpha[k][j]=alpha[j][k];
		}


	private static void gaussj(double a[][], int n, double b[][], int m) throws Exception {
//	Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a[1..n][1..n]
//	is the input matrix. b[1..n][1..m] is input containing the m right-hand side vectors. On
//	output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution
//	vectors.
		int[] indxc = new int[n];	// The integer arrays ipiv, indxr, and indxc are
		int[] indxr = new int[n];	// used for bookkeeping on the pivoting.
		int[] ipiv = new int[n];
		for (int j=0; j<n; j++)
			ipiv[j]=0;  

		for (int i=0; i<n; i++) {// This is the main loop over the columns to be reduced.
			int irow = 0;
			int icol = 0;
			double big=0.0;
			for (int j=0; j<n; j++)	// This is the outer loop of the search for a pivot element.
				if (ipiv[j] != 1)
					for (int k=0; k<n; k++) {
						if (ipiv[k] == 0) {
							if (Math.abs(a[j][k]) >= big) {
								big=Math.abs(a[j][k]);
								irow=j;
								icol=k;
								}
							}
						}
			++(ipiv[icol]);
		/*	We now have the pivot element, so we interchange rows, if needed, to put the pivot
			element on the diagonal. The columns are not physically interchanged, only relabeled:
			indxc[i], the column of the ith pivot element, is the ith column that is reduced, while
			indxr[i] is the row in which that pivot element was originally located. If indxr[i] =
			indxc[i] there is an implied column interchange. With this form of bookkeeping, the
			solution b’s will end up in the correct order, and the inverse matrix will be scrambled
			by columns.	*/
			if (irow != icol) {
				for (int l=0;l<n;l++)
					SWAP(a[irow][l],a[icol][l]);
				for (int l=0;l<m;l++)
					SWAP(b[irow][l],b[icol][l]);
				}
			indxr[i]=irow;	/*	We are now ready to divide the pivot row by the
										pivot element, located at irow and icol.*/
			indxc[i]=icol;
			if (a[icol][icol] == 0.0)
				nrerror("gaussj: Singular Matrix");

			double pivinv=1.0/a[icol][icol];
			a[icol][icol]=1.0;
			for (int l=0;l<n;l++)
				a[icol][l] *= pivinv;
			for (int l=0;l<m;l++)
				b[icol][l] *= pivinv;
			for (int ll=0;ll<n;ll++)	// Next, we reduce the rows...
				if (ll != icol) {		// ...except for the pivot one, of course.
					double dum=a[ll][icol];
					a[ll][icol]=0.0;
					for (int l=0;l<n;l++)
						a[ll][l] -= a[icol][l]*dum;
					for (int l=0;l<m;l++)
						b[ll][l] -= b[icol][l]*dum;
					}
			}
		/*	This is the end of the main loop over columns of the reduction. It only remains to unscramble
			the solution in view of the column interchanges. We do this by interchanging pairs of
			columns in the reverse order that the permutation was built up.	*/
		for (int l=n-1;l>=0;l--) {
			if (indxr[l] != indxc[l])
				for (int k=0;k<n;k++)
					SWAP(a[k][indxr[l]],a[k][indxc[l]]);
			}	// And we are done.
		}


	private static double pythag(double a, double b) {
//	Computes (a2 + b2)1/2 without destructive underflow or overflow.
		double absa = Math.abs(a);
		double absb = Math.abs(b);
		return (absa > absb) ? absa*Math.sqrt(1.0+SQR(absb/absa))
							 : ((absb == 0.0) ? 0.0 : absb*Math.sqrt(1.0+SQR(absa/absb)));
		}


	private static double SQR(double a) {
		return (a == 0.0) ? 0.0 : a*a;
		}


	private static double SIGN(double a, double b) {
		return (b >= 0.0) ? Math.abs(a) : -Math.abs(a);
		}


	private static void SWAP(double a, double b) {
		double temp = a;
		a = b;
		b = temp;
		}


	private static void nrerror(String err) throws Exception {
		throw new Exception("Numerical Recipes run-time error: "+err);
		}
	}
