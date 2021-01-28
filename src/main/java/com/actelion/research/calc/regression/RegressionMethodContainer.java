package com.actelion.research.calc.regression;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.regression.gaussianprocess.GaussianProcessRegression;
import com.actelion.research.calc.regression.knn.KNNRegression;
import com.actelion.research.calc.regression.linear.pls.PLSRegressionModelCalculator;
import com.actelion.research.calc.regression.linear.pls.boxcox.PLSBoxCoxY;
import com.actelion.research.calc.regression.median.MedianRegression;
import com.actelion.research.calc.regression.neuralnetwork.NeuralNetworkRegression;
import com.actelion.research.calc.regression.randomforest.RandomForestRegression;
import com.actelion.research.calc.regression.svm.SVMRegression;
import com.actelion.research.util.datamodel.ModelXYIndex;

/**
 * RegressionMethodContainer
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 03.12.18.
 */
public class RegressionMethodContainer {


    private ARegressionMethod regressionMethod;

    private MedianRegression calculatorMedianRegression;

    private PLSRegressionModelCalculator calculatorPLS;

    private PLSBoxCoxY calculatorPLSBoxCox;

    private KNNRegression calculatorKNNRegression;

    private SVMRegression calculatorSVMRegression;

    private RandomForestRegression randomForestRegression;

    private GaussianProcessRegression gaussianProcessRegression;

    private NeuralNetworkRegression neuralNetworkRegression;


    public RegressionMethodContainer() {

        calculatorMedianRegression = new MedianRegression();

        calculatorPLS = new PLSRegressionModelCalculator();

        calculatorPLSBoxCox = new PLSBoxCoxY();

        calculatorKNNRegression = new KNNRegression();

        calculatorSVMRegression = new SVMRegression();

        randomForestRegression = new RandomForestRegression();

        gaussianProcessRegression = new GaussianProcessRegression();

        neuralNetworkRegression = new NeuralNetworkRegression();

        regressionMethod = calculatorPLS;

    }

    public MedianRegression getMedian(){
        return calculatorMedianRegression;
    }

    public PLSRegressionModelCalculator getPLS(){
        return calculatorPLS;
    }

    public PLSBoxCoxY getPLSBoxCox(){
        return calculatorPLSBoxCox;
    }

    public KNNRegression getKNN(){
        return calculatorKNNRegression;
    }

    public SVMRegression getSVM(){
        return calculatorSVMRegression;
    }

    public RandomForestRegression getRandomForestRegression(){
        return randomForestRegression;
    }

    public GaussianProcessRegression getGaussianProcessRegression() {
        return gaussianProcessRegression;
    }

    public NeuralNetworkRegression getNeuralNetworkRegression() {
        return neuralNetworkRegression;
    }

    public Matrix createModel(ModelXYIndex modelXYIndexTrain) {
        return regressionMethod.createModel(modelXYIndexTrain);
    }

    public Matrix calculateYHat(Matrix X) {
        return regressionMethod.calculateYHat(X);
    }

    public void setMethod(String txt){

        if(calculatorMedianRegression.getName().equals(txt.trim())) {
            regressionMethod = calculatorMedianRegression;
        } else if(calculatorKNNRegression.getName().equals(txt.trim())) {
            regressionMethod = calculatorKNNRegression;
        } else if(calculatorPLSBoxCox.getName().equals(txt.trim())) {
            regressionMethod = calculatorPLSBoxCox;
        } else if(calculatorPLS.getName().equals(txt.trim())) {
            regressionMethod = calculatorPLS;
        } else if(calculatorSVMRegression.getName().equals(txt.trim())) {
            regressionMethod = calculatorSVMRegression;
        } else if(randomForestRegression.getName().equals(txt.trim())) {
            regressionMethod = randomForestRegression;
        } else if(gaussianProcessRegression.getName().equals(txt.trim())) {
            regressionMethod = gaussianProcessRegression;
        } else if(neuralNetworkRegression.getName().equals(txt.trim())) {
            regressionMethod = neuralNetworkRegression;
        } else {
            regressionMethod = calculatorPLS;
            System.err.println("RegressionMethodFactory setMethod unknown model " + txt + " set to " + calculatorPLS.getName());
        }
    }


    public ARegressionMethod getRegressionMethod(){
        return regressionMethod;
    }

    public ARegressionMethod getRegressionMethod(String txt){
        setMethod(txt);
        return regressionMethod;
    }

    public static ARegressionMethod createRegressionMethod(String txt){
        ARegressionMethod regressionMethod = null;

        if(ConstantsRegressionMethods.MODEL_MEDIAN.equals(txt.trim())) {
            regressionMethod = new MedianRegression();
        } else if(ConstantsRegressionMethods.MODEL_KNN.equals(txt.trim())) {
            regressionMethod = new KNNRegression();
        } else if(ConstantsRegressionMethods.MODEL_PLS.equals(txt.trim())) {
            regressionMethod = new PLSRegressionModelCalculator();
        } else if(ConstantsRegressionMethods.MODEL_PLS_POWER.equals(txt.trim())) {
            regressionMethod = new PLSBoxCoxY();
        } else if(ConstantsRegressionMethods.MODEL_SVM.equals(txt.trim())) {
            regressionMethod = new SVMRegression();
        } else if(ConstantsRegressionMethods.MODEL_RND_FOREST.equals(txt.trim())) {
            regressionMethod = new RandomForestRegression();
        } else if(ConstantsRegressionMethods.MODEL_GAUSSIAN_PROCESS.equals(txt.trim())) {
            regressionMethod = new GaussianProcessRegression();
        } else if(ConstantsRegressionMethods.MODEL_NEURAL_NETWORK.equals(txt.trim())) {
            regressionMethod = new NeuralNetworkRegression();
        } else {
            throw new RuntimeException("Unknown regression type '" + txt + "'.");
        }

        return regressionMethod;
    }

    /**
     * Sets the parameter
     * @param parameterRegressionMethod
     * @return
     */
    public static ARegressionMethod createRegressionMethod(ParameterRegressionMethod parameterRegressionMethod){

        ARegressionMethod regressionMethod = createRegressionMethod(parameterRegressionMethod.getName());

        regressionMethod.setParameterRegressionMethod(parameterRegressionMethod);

        return regressionMethod;
    }


}
