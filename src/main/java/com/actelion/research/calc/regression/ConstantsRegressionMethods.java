package com.actelion.research.calc.regression;

/**
 * ConstantsRegressionMethods
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 28.11.18.
 */
public class ConstantsRegressionMethods {

    public static final String  MODEL_PLS = "PLS";
    public static final String  MODEL_PLS_POWER = "PLS Power";
    public static final String  MODEL_MEDIAN = "Median";
    public static final String  MODEL_KNN = "KNN regression";
    public static final String  MODEL_SVM = "SVM regression";
    public static final String  MODEL_RND_FOREST = "Random Forest regression";
    public static final String  MODEL_GAUSSIAN_PROCESS = "Gaussian process regression";
    public static final String  MODEL_NEURAL_NETWORK = "Neural network regression";

    // method code that is guaranteed to stay unchanged for all time
    public static final String [] ARR_MODEL_CODE = {"pls", "plsPower", "knn", "svm"};
    public static final String [] ARR_MODEL = {MODEL_PLS, MODEL_PLS_POWER, MODEL_KNN, MODEL_SVM};
    public static final int DEFAULT_MODEL = 3;

    public static final String  RESULT_EVALUATOR_BALANCED_R2_TRAIN_TEST = "BalancedR2TrainTest";

    public static final String  RESULT_EVALUATOR_MAXIMUM_R2_TEST = "MaximumR2Test";

    public static final String [] ARR_RESULT_EVALUATOR = {RESULT_EVALUATOR_BALANCED_R2_TRAIN_TEST, RESULT_EVALUATOR_MAXIMUM_R2_TEST};

}
