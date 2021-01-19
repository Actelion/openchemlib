package com.actelion.research.calc.regression.median;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.regression.ARegressionMethod;
import com.actelion.research.util.datamodel.ModelXYIndex;

/**
 * MedianRegression
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 13.02.19.
 *
 * The zero model.
 */
public class MedianRegression extends ARegressionMethod<ParameterMedian> {

    private Matrix maMedian;

    public MedianRegression() {
        setParameterRegressionMethod(new ParameterMedian());
    }

    public Matrix createModel(ModelXYIndex modelXYIndexTrain) {

        maMedian = modelXYIndexTrain.Y.getMedianCols();

        Matrix maYHat = calculateYHat(modelXYIndexTrain.X);

        return maYHat;
    }

    public Matrix calculateYHat(Matrix X){

        int rows = X.rows();
        int cols = maMedian.cols();

        Matrix maYHat = new Matrix(rows, cols);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                maYHat.set(i,j,maMedian.get(0,j));
            }
        }

        return maYHat;
    }

    @Override
    public double calculateYHat(double[] arrRow) {
        return maMedian.get(0,0);
    }
}
