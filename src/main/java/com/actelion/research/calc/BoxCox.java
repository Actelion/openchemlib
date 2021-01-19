package com.actelion.research.calc;


import com.actelion.research.util.ArrayUtils;

import java.text.DecimalFormat;

/**
 * 
 * BoxCox
 * @author Modest von Korff
 * 2005 MvK Start implementation
 * http://www.itl.nist.gov/div898/handbook/eda/section3/boxcox.htm
 * https://en.wikipedia.org/wiki/Power_transform
 */
public class BoxCox {

	private static final double TINY = 0.0000001;

	private double lambda;

	private boolean zero;

	public BoxCox(double lambda){
		setLambda(lambda);
	}

	public void setLambda(double lambda) {
		this.lambda = lambda;
		if(lambda<TINY){
			zero=true;
		}else{
			zero=false;
		}
	}

	public double get(double val){
		double r = 0;

		boolean negative=false;
		if(val<0){
			val*=-1;
			negative=true;
		}

		if(zero){
			r = Math.log(val);
		} else {
			r = (Math.pow(val, lambda)-1.0)/lambda;
		}

		if(negative){
			r *= -1;
		}

		return r;
	}

	public double inverse(double v) {

		boolean negative=false;
		if(v<0){
			v*=-1;
			negative=true;
		}

		double r = 0;

		if(zero){
			r = Math.exp(r);
		} else {
			r = Math.pow(((v * lambda) + 1), 1/lambda);

			if(negative){
				r *= -1;
			}
		}




		return r;
	}

	public double getLambda() {
		return lambda;
	}

	public static Matrix transform(Matrix A, BoxCox boxCox) {

		int rows = A.rows();
		int cols = A.cols();

		Matrix B = new Matrix(rows, cols);

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				double v = A.get(i,j);

				double bc = boxCox.get(v);

				B.set(i,j, bc);
			}
		}
		return B;
	}

	public static Matrix reTransform(Matrix A, BoxCox boxCox) {

		int rows = A.rows();
		int cols = A.cols();

		Matrix B = new Matrix(rows, cols);

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				double v = A.get(i,j);

				double bc = boxCox.inverse(v);

				B.set(i,j,bc);
			}
		}
		return B;
	}



	public static void main(String[] args) {
		test01();
	}

	public static void test01() {

		double [] arr = {-64,-32,-16,-8,-4,-2,-1,0,1,2,4,8,16,32,64};

		double lambdaStart = -3;
		double lambdaEnd = 3;

		double step = 0.1;

		DecimalFormat df = new DecimalFormat("0.000");
		DecimalFormat dfLambda = new DecimalFormat("0.00");

		for (double lambda = lambdaStart; lambda < lambdaEnd+step; lambda+=step) {

			BoxCox boxCox = new BoxCox(lambda);

			double [] arrTrans = new double[arr.length];
			double [] arrReTrans = new double[arr.length];

			for (int i = 0; i < arr.length; i++) {
				arrTrans[i] = boxCox.get(arr[i]);

				arrReTrans[i]=boxCox.inverse(arrTrans[i]);
			}

			System.out.println(dfLambda.format(lambda) + "\t" + ArrayUtils.toString(arrTrans, df));
			System.out.println(dfLambda.format(lambda) + "\t" + ArrayUtils.toString(arrReTrans, df));


		}
	}


	
}
