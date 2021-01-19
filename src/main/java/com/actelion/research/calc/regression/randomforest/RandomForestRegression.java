package com.actelion.research.calc.regression.randomforest;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.regression.ARegressionMethod;
import com.actelion.research.util.datamodel.ModelXYIndex;
import smile.regression.RandomForest;

/**
 * RandomForestRegression
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 14.01.19.
 */
public class RandomForestRegression extends ARegressionMethod<ParameterRandomForest> implements Comparable<RandomForestRegression> {

    public static final int MIN_NUM_VAR_SPLIT = 3;

    private RandomForest forest;


    public RandomForestRegression() {
        setParameterRegressionMethod(new ParameterRandomForest());
        // To prevent multi-core execution on Random Forest level
        // On the grid permission is denied.
        try {
            System.setProperty("smile.threads", "1");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public RandomForestRegression(ParameterRandomForest parameterRandomForest) {
        setParameterRegressionMethod(parameterRandomForest);
    }

    @Override
    public Matrix createModel(ModelXYIndex modelXYIndexTrain) {

        Matrix YHat = null;
        try {
            ParameterRandomForest parameterRandomForest = getParameter();

            int mTry = (int)(modelXYIndexTrain.X.cols() * parameterRandomForest.getFractionMTry()+0.5);

            if(mTry<MIN_NUM_VAR_SPLIT){
                mTry=MIN_NUM_VAR_SPLIT;
            }

            forest = new RandomForest(modelXYIndexTrain.X.getArray(), modelXYIndexTrain.Y.getColAsDouble(0),
                    parameterRandomForest.getNumberOfTrees(),
                    parameterRandomForest.getMaxNodes(),
                    parameterRandomForest.getNodeSize(),
                    mTry);

            YHat = calculateYHat(modelXYIndexTrain.X);
        } catch (Exception e) {
            e.printStackTrace();
        }

        return YHat;
    }

    @Override
    public Matrix calculateYHat(Matrix X) {

        double [] arrY = new double[X.rows()];

        for (int i = 0; i < X.rows(); i++) {

            double [] arrRow = X.getRow(i);

            double y = forest.predict(arrRow);

            arrY[i]=y;
        }
        return new Matrix(false, arrY);
    }

    @Override
    public double calculateYHat(double[] arrRow) {

        double y;

        synchronized (this) {
            y = forest.predict(arrRow);
        }
        return y;
    }

    @Override
    public int compareTo(RandomForestRegression o) {
        return getParameter().compareTo(o.getParameter());
    }


}
