package com.actelion.research.util;


/**
 * Smoothing Spline Interpolator based on the algorithm described at
 * http://www.qmw.ac.uk/~ugte133/book/11_tsd/splines.pdf
 * <br>
 * <br>
 * The Smoothing Spline is used to minimize Sum(sqr((Si-yi)/sigmai)) + lambda*Sum(sqr(S''i))
 * <li>If lambda=0 (default), this is equivalent to the cubic spline interpolation
 * <li>If lambda=Infinity, this is equivalent to the least square fitting
 * 
 * <pre>FastSpline spline = new SmoothingSplineInterpolator().interpolate(X, Y);</pre>
 * 
 * @author freyssj
 */
public class SmoothingSplineInterpolator  {

	private double lambda = 0;
	private double[] sigma = null;
	
	private double residuals;
	private double smoothing;
	
	/**
	 * 
	 */
	private static void quincunx(int n, double[] u, double[] v, double[] w, double[] q) {			
		//Factorisation	
		u[0] = 0;
		v[0] = 0;
		w[0] = 0;
		
		for (int j = 1; j <= n-1; j++) {
			u[j] = u[j] - (j-2>=0? u[j-2]*sqr(w[j-2]): 0) - u[j-1]*sqr(v[j-1]);
			v[j] = (v[j] - u[j-1] * v[j-1] * w[j-1]) / u[j];
			w[j] = w[j] / u[j];			
		}
		
		//Forward Substitution
		q[0] = 0;
		for (int j = 1; j <= n-1; j++) {
			q[j] = q[j] - v[j-1] * q[j-1] - (j-2>=0? w[j-2]*q[j-2] :0);			
		}
		for (int j = 1; j <= n-1; j++) {
			q[j] = q[j] / u[j];
		}
				
		//Back Substitution
		q[n-1] = q[n-1];
		q[n-2] = q[n-2] - v[n-2] * q[n-2+1];
		for (int j = n-3; j >=1; j--) {
			q[j] = q[j] - v[j] * q[j+1] - w[j] * q[j+2];
		}
	}
	
	/**
	 * @see org.apache.commons.math.analysis.UnivariateRealInterpolator#interpolate(double[], double[])
	 */
	public FastSpline interpolate(double[] x, double[] y) {

		if (x.length != y.length) throw new IllegalArgumentException("Dataset arrays must have same length.");
		if (x.length < 3) throw new IllegalArgumentException("At least 3 datapoints are required to compute a spline interpolant");
		if (sigma!=null && sigma.length != x.length) throw new IllegalArgumentException("Sigma and dataset arrays must have same length");

		for (int i = 0; i < x.length-1; i++) {
			if(x[i+1]<=x[i]) throw new IllegalArgumentException("the X must be strictly increasing");
		}
		

		int n = x.length;	
		double[] h = new double[n-1];
		double[] r = new double[n];
		double[] f = new double[n];
		double[] p = new double[n];
		double[] q = new double[n];
		double[] u = new double[n];
		double[] v = new double[n];
		double[] w = new double[n];
		
		for (int i = 0; i < n-1; i++) {		
			h[i] = x[i + 1] - x[i];
			r[i] = 3d / h[i];
		}
		for (int i = 1; i < n-1; i++) {		
			f[i] = -(r[i-1] + r[i]);
			p[i] = 2d * (x[i + 1] - x[i - 1]);
			q[i] = 3d * (y[i + 1] - y[i]) / h[i]
				 - 3d * (y[i] - y[i - 1]) / h[i - 1];
		}
		
		r[n-1] = 0;
		f[n-1] = 0;
		for (int i = 1; i <= n-1; i++) {		
			u[i] = sqr(r[i-1])*sigma(i-1) + sqr(f[i])*sigma(i) + sqr(r[i])*sigma(i + 1);
			u[i] = lambda * u[i] + p[i];
			if(i<n-1) v[i] = f[i] * r[i] * sigma(i) + r[i] * f[i + 1] * sigma(i + 1);
			if(i<n-1) v[i] = lambda * v[i] + h[i];
			if(i<n-1) w[i] = lambda * r[i] * r[i + 1] * sigma(i + 1);
		}

		//Solve the system
		//lambda=0 => u=p, v=h, w=0
		quincunx(n, u, v, w, q);
		
		//Spline Parameters
		residuals = smoothing = 0;
		double a[] = new double[n];
		double b[] = new double[n];
		double c[] = new double[n];
		double d[] = new double[n];
		d[0] = y[0] - lambda * r[0] * q[1] * sigma(0);
		d[1] = y[1] - lambda * (f[1] * q[1] + r[1] * q[2]) * sigma(1);
		a[0] = q[1] / (3d * h[0]);
		b[0] = 0;
		c[0] = (d[1] - d[0])/h[0] - q[1] * h[0]/3d;		
		r[0] = 0;
		
		for (int j = 1; j < n-1; j++) {			
			a[j] = (q[j + 1]-q[j]) / (3d * h[j]);
			b[j] = q[j];
			c[j] = (q[j] + q[j-1]) * h[j-1] + c[j-1];
			d[j] = r[j-1] * q[j-1] + f[j] * q[j] + r[j] * q[j+1];
			d[j] = y[j] - lambda * d[j] * sigma(j);
		}
		for (int j = 0; j < n-1; j++) {
			if(sigma(j)>0) residuals += (d[j]-y[j]) * (d[j]-y[j]) / (sigma(j) * sigma(j));
			smoothing += a[j]*a[j];
		}

		FastSpline.Polynome polynomials[] = new FastSpline.Polynome[n-1];
		for (int i = 0; i < polynomials.length; i++) {
			polynomials[i] = new FastSpline.Polynome(new double[] {d[i], c[i], b[i], a[i]});
		}
		return new FastSpline(x, polynomials);
	}
	
	public double getResiduals() {
		return residuals;
	}
	
	public double getSmoothing() {
		return smoothing;
	}

	public double getLambda() {
		return lambda;
	}

	public double[] getSigma() {
		return sigma;
	}

	public void setLambda(double d) {
		lambda = d;
	}

	public void setSigma(double[] ds) {
		sigma = ds;
	}	
	
	private static final double sqr(double v) {
		return v*v;
	}
	private final double sigma(int index) {
		if(sigma==null || index>=sigma.length) return 1;
		return sigma[index];
	}
	
}
