package com.actelion.research.calc.regression.gaussianprocess;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.regression.ARegressionMethod;
import com.actelion.research.util.datamodel.ModelXYIndex;
import smile.clustering.KMeans;
import smile.math.Math;
import smile.math.kernel.GaussianKernel;
import smile.math.kernel.MercerKernel;

/**
 * GaussianProcessRegression
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 01.04.19.
 */
public class GaussianProcessRegression extends ARegressionMethod<ParameterGaussianProcess> implements Comparable<GaussianProcessRegression> {




    private static final int MIN_K = 3;
    private static final int K_DEVISOR = 10;
    // private static final int K = 11;

    private smile.regression.GaussianProcessRegression<double[]> gaussianProcessRegression;


    public GaussianProcessRegression() {
        setParameterRegressionMethod(new ParameterGaussianProcess());

        // To prevent multi-core execution on Random Forest level
        // On the grid the permissions are denied.
        try {
            System.setProperty("smile.threads", "1");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public GaussianProcessRegression(ParameterGaussianProcess parameterGaussianProcess) {
        setParameterRegressionMethod(parameterGaussianProcess);
    }

    public void setLambda(double lambda){
        getParameter().setLambda(lambda);
    }



    @Override
    public Matrix createModel(ModelXYIndex modelXYIndexTrain) {

        Matrix YHat = null;
        try {

            ParameterGaussianProcess parameterGaussianProcess = getParameter();

            int rows = modelXYIndexTrain.X.rows();

            if(modelXYIndexTrain.Y.cols()!=1){
                throw new RuntimeException("Only one column for y is allowed!");
            } else if(rows < MIN_K){
                throw new RuntimeException("Unsufficient number of objects for regression.");
            }

            double [][] X = modelXYIndexTrain.X.getArray();

            double [] y = modelXYIndexTrain.Y.getColAsDouble(0);

            int k = rows / K_DEVISOR;

            if(rows>1000){
                k = rows / 100;
            } else if(rows>10000){
                k = rows / 1000;
            }

            if(k < MIN_K){
                k=MIN_K;
            }

            // System.out.println("GaussianProcessRegression rows train " + rows + ", k " + k + ".");

            KMeans kmeans = new KMeans(X, k, 10);

            double[][] centers = kmeans.centroids();
            double r0 = 0.0;
            for (int l = 0; l < centers.length; l++) {
                for (int j = 0; j < l; j++) {
                    r0 += Math.distance(centers[l], centers[j]);
                }
            }
            r0 /= (2 * centers.length);

            MercerKernel mercerKernel = new GaussianKernel(r0);

            gaussianProcessRegression = new smile.regression.GaussianProcessRegression(X, y, mercerKernel, parameterGaussianProcess.getLambda());

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

            double y = gaussianProcessRegression.predict(arrRow);

            arrY[i]=y;
        }
        return new Matrix(false, arrY);
    }

    @Override
    public double calculateYHat(double[] arrRow) {

        double yHat;

        synchronized (this) {
            yHat = gaussianProcessRegression.predict(arrRow);
        }

        return yHat;
    }

    @Override
    public int compareTo(GaussianProcessRegression o) {
        return getParameter().compareTo(o.getParameter());
    }


}
