package com.actelion.research.calc.regression.svm;

import org.machinelearning.svm.libsvm.svm_parameter;


/**
 * SVMParameterHelper
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 27.11.18.
 */
public class SVMParameterHelper {

    public static final String SVM_TYPE_C = "C_SVC";
    public static final String SVM_TYPE_NU = "NU_SVC";
    public static final String SVM_TYPE_ONE_CLASS  = "OneClass";
    public static final String SVM_TYPE_EPSILON_SVR = "EpsSVR";
    public static final String SVM_TYPE_NU_SVR = "NU_SVR";



    public static final String KERNEL_TYPE_LINEAR = "Linear";
    public static final String KERNEL_TYPE_POLY = "Poly";
    public static final String KERNEL_TYPE_RBF = "RBF";
    public static final String KERNEL_TYPE_SIGMOID = "Sigmoid";
    public static final String KERNEL_TYPE_PRECOMPUTED = "Precomputed";

    public static final int DEGREE_ANALYTICALLY_PARAMETER_CALC = -1;



    public static svm_parameter standard() {

        svm_parameter param = new svm_parameter();
        // default values
        param.svm_type = svm_parameter.C_SVC;
        param.kernel_type = svm_parameter.RBF;
        param.degree = 3;
        param.gamma = 0;	// 1/num_features
        param.coef0 = 0;
        param.nu = 0.5;
        param.cache_size = 100;
        param.C = 1;
        param.eps = 1e-3;
        param.p = 0.1;
        param.shrinking = 1;
        param.probability = 0;
        param.nr_weight = 0;
        param.weight_label = new int[0];
        param.weight = new double[0];

        return param;
    }

    /**
     * Regression epsilon-SVR
     * Degree set to trigger analytical parameter determination in createModel(....)
     * @return
     */
    public static svm_parameter regressionEpsilonSVR() {

        svm_parameter param = new svm_parameter();
        // default values
        param.svm_type = svm_parameter.EPSILON_SVR;
        param.kernel_type = svm_parameter.RBF;
        param.degree = DEGREE_ANALYTICALLY_PARAMETER_CALC;
        param.gamma = 0;	// 1/num_features
        param.coef0 = 0;
        param.nu = 0.5;
        param.cache_size = 100;
        param.C = 130;
        param.eps = 5;
        param.p = 0.1;
        param.shrinking = 1;
        param.probability = 0;
        param.nr_weight = 0;
        param.weight_label = new int[0];
        param.weight = new double[0];

        return param;
    }



//  Until 01.11.2019
//    public static svm_parameter regressionEpsilonSVR() {
//
//        svm_parameter param = new svm_parameter();
//        // default values
//        param.svm_type = svm_parameter.EPSILON_SVR;
//        param.kernel_type = svm_parameter.RBF;
//        param.degree = 3;
//        param.gamma = 0;	// 1/num_features
//        param.coef0 = 0;
//        param.nu = 0.5;
//        param.cache_size = 100;
//        param.C = 1;
//        param.eps = 1e-3;
//        param.p = 0.1;
//        param.shrinking = 1;
//        param.probability = 0;
//        param.nr_weight = 0;
//        param.weight_label = new int[0];
//        param.weight = new double[0];
//
//        return param;
//    }

    /**
     * Regression nu-SVR
     *
     * @param nu, default was 0.5
     * @return
     */
    public static svm_parameter regressionNuSVR(double nu) {

        svm_parameter param = new svm_parameter();
        // default values
        param.svm_type = svm_parameter.NU_SVR;
        param.kernel_type = svm_parameter.RBF;
        param.degree = 3;
        param.gamma = 1;	// 1/num_features
        param.coef0 = 0;
        param.nu = nu;
        param.cache_size = 100;
        param.C = 1;
        param.eps = 1e-3;
        param.p = 0.1;
        param.shrinking = 1;
        param.probability = 0;
        param.nr_weight = 0;
        param.weight_label = new int[0];
        param.weight = new double[0];

        return param;
    }


    public static String getSVMType(int svmType) {
        String strType = "";

        switch (svmType){
            case svm_parameter.C_SVC:
                strType = SVM_TYPE_C;
                break;
            case svm_parameter.NU_SVC:
                strType = SVM_TYPE_NU;
                break;
            case svm_parameter.ONE_CLASS:
                strType = SVM_TYPE_ONE_CLASS;
                break;
            case svm_parameter.NU_SVR:
                strType = SVM_TYPE_NU_SVR;
                break;
            case svm_parameter.EPSILON_SVR:
                strType = SVM_TYPE_EPSILON_SVR;
                break;
        }

        return strType;
    }

    public static int getSVMType(String strSVMType) {
        int svmType = -1;

        if(SVM_TYPE_C.equals(strSVMType)) {
            svmType = svm_parameter.C_SVC;

        }else if(SVM_TYPE_NU.equals(strSVMType)) {
            svmType = svm_parameter.NU_SVC;

        }else if(SVM_TYPE_ONE_CLASS.equals(strSVMType)) {
            svmType = svm_parameter.ONE_CLASS;

        }else if(SVM_TYPE_NU_SVR.equals(strSVMType)) {
            svmType = svm_parameter.NU_SVR;

        }else if(SVM_TYPE_EPSILON_SVR.equals(strSVMType)) {
            svmType = svm_parameter.EPSILON_SVR;
        } else {
            throw new RuntimeException("Unknown svm type " + strSVMType + ".");
        }

        return svmType;
    }

    public static String getKernelType(int svmType) {
        String strType = "";

        switch (svmType){
            case svm_parameter.LINEAR:
                strType = KERNEL_TYPE_LINEAR;
                break;
            case svm_parameter.POLY:
                strType = KERNEL_TYPE_POLY;
                break;
            case svm_parameter.RBF:
                strType = KERNEL_TYPE_RBF;
                break;
            case svm_parameter.SIGMOID:
                strType = KERNEL_TYPE_SIGMOID;
                break;
            case svm_parameter.PRECOMPUTED:
                strType = KERNEL_TYPE_PRECOMPUTED;
                break;
        }

        return strType;
    }

    public static int getKernelType(String strKernelType) {
        int svmType = -1;

        if(KERNEL_TYPE_LINEAR.equals(strKernelType)) {
            svmType = svm_parameter.LINEAR;

        }else if(KERNEL_TYPE_POLY.equals(strKernelType)) {
            svmType = svm_parameter.POLY;

        }else if(KERNEL_TYPE_RBF.equals(strKernelType)) {
            svmType = svm_parameter.RBF;

        }else if(KERNEL_TYPE_SIGMOID.equals(strKernelType)) {
            svmType = svm_parameter.SIGMOID;

        }else if(KERNEL_TYPE_PRECOMPUTED.equals(strKernelType)) {
            svmType = svm_parameter.PRECOMPUTED;
        } else {
            throw new RuntimeException("Unknown kernel type " + strKernelType + ".");
        }

        return svmType;
    }

}
