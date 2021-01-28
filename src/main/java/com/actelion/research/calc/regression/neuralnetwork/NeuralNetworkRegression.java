package com.actelion.research.calc.regression.neuralnetwork;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.regression.ARegressionMethod;
import com.actelion.research.util.datamodel.ModelXYIndex;
import smile.regression.NeuralNetwork;

/**
 * NeuralNetworkRegression
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 02.04.19.
 */
public class NeuralNetworkRegression extends ARegressionMethod<ParameterNeuralNetwork> implements Comparable<NeuralNetworkRegression> {

    private NeuralNetwork neuralNetwork;

    public NeuralNetworkRegression() {
        setParameterRegressionMethod(new ParameterNeuralNetwork());
        // To prevent multi-core execution on Random Forest level
        System.setProperty("smile.threads", "1");
    }

    public NeuralNetworkRegression(ParameterNeuralNetwork parameterNeuralNetwork) {
        setParameterRegressionMethod(parameterNeuralNetwork);
    }

    @Override
    public Matrix createModel(ModelXYIndex modelXYIndexTrain) {

        Matrix YHat = null;
        try {

            ParameterNeuralNetwork parameterNeuralNetwork = getParameter();

            if(modelXYIndexTrain.Y.cols()!=1){
                throw new RuntimeException("Only one column for y is allowed!");
            }

            double [][] X = modelXYIndexTrain.X.getArray();

            double [] y = modelXYIndexTrain.Y.getColAsDouble(0);

            int [] arrNetworkArchitecture = new int[parameterNeuralNetwork.getArrInnerLayerArchitecture().length + 2];

            System.arraycopy(parameterNeuralNetwork.getArrInnerLayerArchitecture(), 0, arrNetworkArchitecture, 1, parameterNeuralNetwork.getArrInnerLayerArchitecture().length);

            arrNetworkArchitecture[0]=modelXYIndexTrain.X.cols();

            arrNetworkArchitecture[arrNetworkArchitecture.length-1]=1;

            neuralNetwork = new NeuralNetwork(parameterNeuralNetwork.getActivationFunction(), arrNetworkArchitecture);

            neuralNetwork.learn(X,y);

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

            double y = neuralNetwork.predict(arrRow);

            arrY[i]=y;
        }
        return new Matrix(false, arrY);
    }

    @Override
    public double calculateYHat(double[] arrRow) {

        double y;

        synchronized (this) {
            y = neuralNetwork.predict(arrRow);
        }

        return y;
    }

    @Override
    public int compareTo(NeuralNetworkRegression o) {
        return getParameter().compareTo(o.getParameter());
    }


}
