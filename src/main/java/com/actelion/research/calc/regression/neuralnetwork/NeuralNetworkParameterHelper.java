package com.actelion.research.calc.regression.neuralnetwork;

import smile.regression.NeuralNetwork;


/**
 * NeuralNetworkParameterHelper
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 02.04.19.
 */
public class NeuralNetworkParameterHelper {


    public static final String ACTIVATION_FCT_LOG = "LogisticSigmoid";

    public static final String ACTIVATION_FCT_TAN = "Tangens";


    public static String getActivationFunctionName(NeuralNetwork.ActivationFunction activationFunction) {
        String type = "";

        switch (activationFunction){
            case LOGISTIC_SIGMOID:
                type = ACTIVATION_FCT_LOG;
                break;
            case TANH:
                type = ACTIVATION_FCT_TAN;
                break;
        }

        return type;
    }

    public static NeuralNetwork.ActivationFunction getActivationFunction(String type) {
        NeuralNetwork.ActivationFunction activationFunction = null;

        if(ACTIVATION_FCT_LOG.equals(type)) {
            activationFunction = NeuralNetwork.ActivationFunction.LOGISTIC_SIGMOID;

        }else if(ACTIVATION_FCT_TAN.equals(type)) {
            activationFunction = NeuralNetwork.ActivationFunction.TANH;

        } else {
            throw new RuntimeException("Unknown activation function type " + type + ".");
        }

        return activationFunction;
    }


}
