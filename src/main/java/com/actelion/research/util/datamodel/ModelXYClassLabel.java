
package com.actelion.research.util.datamodel;

import java.util.ArrayList;
import java.util.List;


/**
 * 
 * 
 * ModelXYClassLabel
 * @author Modest von Korff
 * 02.01.2019 MvK: Start implementation
 */
public class ModelXYClassLabel extends ModelXYIndex {


	// Class label. To enable a controlled leave out cross validation.
	public List<Integer> liClassLabel;

	public ModelXYClassLabel() {

	}

	/**
	 * Deep copy constructor
	 * @param dataXY
	 */
	public ModelXYClassLabel(ModelXYClassLabel dataXY) {
		super(dataXY);

		if(dataXY.liClassLabel!=null) {
			liClassLabel = new ArrayList<Integer>(dataXY.X.rows());

			for (int index : dataXY.liClassLabel) {
				liClassLabel.add(index);
			}
		}
	}

	public ModelXYClassLabel sub(List<Integer> liIndexRow){

		ModelXYClassLabel modelDataXY = new ModelXYClassLabel();
		
		modelDataXY.X = X.getSubMatrix(liIndexRow);
		
		modelDataXY.Y = Y.getSubMatrix(liIndexRow);
		
		modelDataXY.liIndex = new ArrayList<>(liIndexRow);

		return modelDataXY;
	}

	public void sub(List<Integer> liIndexRow, ModelXYClassLabel modelXYIndex){


		for (int i = 0; i < liIndexRow.size(); i++) {

			int index = liIndexRow.get(i);

			double [] a = X.getRow(index);

			modelXYIndex.X.setRow(i, a);

			double [] b = Y.getRow(index);

			modelXYIndex.Y.setRow(i, b);
		}
	}

	@Override
	public ModelXYClassLabel getDeepClone() {
		return new ModelXYClassLabel(this);
	}

}