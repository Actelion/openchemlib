package com.actelion.research.calc.regression.neuralnetwork;

import com.actelion.research.calc.ArrayUtilsCalc;
import com.actelion.research.calc.regression.ConstantsRegressionMethods;
import com.actelion.research.calc.regression.ParameterRegressionMethod;
import smile.regression.NeuralNetwork;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * ParameterNeuralNetwork
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 02.04.19.
 */
public class ParameterNeuralNetwork extends ParameterRegressionMethod {

    // int [] arrNetworkArchitecture = new int[]{modelXYIndexTrain.X.cols(), 100, 30, 10,1};

    public static final String TAG_ACTIVATION_FCT="ActivationFunction";
    public static final String TAG_INNERLAYER_ARCITECTURE="InnerLayerArchitecture";



    private int [] arrInnerLayerArchitecture;

    private NeuralNetwork.ActivationFunction activationFunction;


    public ParameterNeuralNetwork() {
        super(ConstantsRegressionMethods.MODEL_NEURAL_NETWORK);

        setActivationFunction(NeuralNetwork.ActivationFunction.LOGISTIC_SIGMOID);

        setArrInnerLayerArchitecture(new int[0]);
    }


    @Override
    public int compareTo(ParameterRegressionMethod o) {

        int cmp = 0;

        ParameterNeuralNetwork p = (ParameterNeuralNetwork)o;


        return cmp;
    }

    public void setArrInnerLayerArchitecture(int[] arrInnerLayerArchitecture) {
        this.arrInnerLayerArchitecture = arrInnerLayerArchitecture;
        properties.put(TAG_INNERLAYER_ARCITECTURE, ArrayUtilsCalc.toString(arrInnerLayerArchitecture));
    }

    public void setActivationFunction(NeuralNetwork.ActivationFunction activationFunction) {

        properties.put(TAG_ACTIVATION_FCT, NeuralNetworkParameterHelper.getActivationFunctionName(activationFunction));

        this.activationFunction = activationFunction;
    }

    public NeuralNetwork.ActivationFunction getActivationFunction() {
        return activationFunction;
    }

    public int [] getArrInnerLayerArchitecture() {
        return arrInnerLayerArchitecture;
    }

    public static List<String> getHeader(){

        List<String> li = ParameterRegressionMethod.getHeader();

        li.add(TAG_ACTIVATION_FCT);
        li.add(TAG_INNERLAYER_ARCITECTURE);

        return li;
    }

    @Override
    public void decodeProperties2Parameter() {

        String sActivationFctType = properties.getProperty(TAG_ACTIVATION_FCT);
        activationFunction = NeuralNetworkParameterHelper.getActivationFunction(sActivationFctType);

        String sArrInnerLayerArchitecture = properties.getProperty(TAG_INNERLAYER_ARCITECTURE);
        arrInnerLayerArchitecture = ArrayUtilsCalc.readIntArray(sArrInnerLayerArchitecture);
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("ParameterNeuralNetwork{");
        sb.append("arrInnerLayerArchitecture=").append(Arrays.toString(arrInnerLayerArchitecture));
        sb.append(", activationFunction=").append(NeuralNetworkParameterHelper.getActivationFunctionName(activationFunction));
        sb.append('}');
        return sb.toString();
    }

    public static void main(String[] args) throws IOException {

        File dir = new File("/home/korffmo1/tmp/tmp00");

        File fiProp = new File(dir, "neuralNetwork.properties");

        ParameterNeuralNetwork parameter = new ParameterNeuralNetwork();

        parameter.setActivationFunction(NeuralNetwork.ActivationFunction.LOGISTIC_SIGMOID);
        parameter.setArrInnerLayerArchitecture(new int[]{512,64,32,8,64});

        parameter.write(fiProp);

        ParameterNeuralNetwork parameterIn = new ParameterNeuralNetwork();

        parameterIn.read(fiProp);

        System.out.println(parameterIn.toString());

    }

}
