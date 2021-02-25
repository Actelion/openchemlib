
package com.actelion.research.util.datamodel;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.MatrixFunctions;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


/**
 * ModelDataXY
 * X and Y are assumed to have the same row dim.
 * @author Modest von Korff
 * Aug 4, 2011 MvK: Start implementation
 */
public class ModelXY implements IModelCloneable<ModelXY> {
		
	public Matrix X;
	
	public Matrix Y;
		
	public ModelXY() {
		
	}

	public ModelXY(int rows, int colsX, int colsY) {

		X = new Matrix(rows, colsX);
		Y = new Matrix(rows, colsY);
	}
	
	/**
	 * @param x
	 * @param y
	 */
	public ModelXY(Matrix x, Matrix y) {
		super();
		X = x;
		Y = y;
	}

	public ModelXY(List<XY> li) {

		double [][] X = new double[li.size()][];
		double [][] Y = new double[li.size()][];

		for (int i = 0; i < li.size(); i++) {

			XY xy = li.get(i);

			X[i]=xy.x;
			Y[i]=xy.y;
		}

		this.X = new Matrix(X);
		this.Y = new Matrix(Y);
	}

	/**
	 * Deep copy constructor
	 * @param dataXY
	 */
	public ModelXY(ModelXY dataXY) {
		deepCopy(dataXY);
	}

	public void deepCopy(ModelXY source) {

		X = new Matrix(source.X.rows(), source.X.cols());

		Y = new Matrix(source.Y.rows(), source.Y.cols());

		X.copy(source.X);

		Y.copy(source.Y);
	}


	
	/**
	 * 
	 * @return rows in X.
	 */
	public int size(){
		return X.rows();
	}


	@Override
	public ModelXY getDeepClone() {
		return new ModelXY(this);
	}

	/**
	 *
	 * @return Deep copy.
	 */
	public List<XY> getAsList(){
		int rows = X.rows();
		List<XY> li = new ArrayList<>(rows);
		for (int i = 0; i < rows; i++) {
			li.add(new XY(X.getRowCopy(i), Y.getRowCopy(i)));
		}
		return li;
	}

	public List<int []> toIntegerListX(){
		int rows = X.rows();
		List<int []> li = new ArrayList<>(rows);
		for (int i = 0; i < rows; i++) {

			double [] a = X.getRow(i);

			int [] b = new int[a.length];

			for (int j = 0; j < a.length; j++) {
				b[j]=(int)a[j];
			}
			li.add(b);
		}
		return li;
	}

	public ModelXY getSortedByY(int colY){

		List<XY>  li = getAsList();

		Collections.sort(li, XY.getComparatorY(colY));

		double [][] arrX = new double[li.size()][];
		double [][] arrY = new double[li.size()][];

		for (int i = 0; i < li.size(); i++) {

			arrX[i]=li.get(i).x;
			arrY[i]=li.get(i).y;

		}

		ModelXY modelSorted = new ModelXY();

		modelSorted.X = new Matrix(arrX);
		modelSorted.Y = new Matrix(arrY);

		return modelSorted;
	}

	public static ModelXY appendRows(ModelXY a, ModelXY b){
			Matrix X = MatrixFunctions.appendRows(a.X, b.X);
			Matrix Y = MatrixFunctions.appendRows(a.Y, b.Y);
			return new ModelXY(X,Y);
	}

}