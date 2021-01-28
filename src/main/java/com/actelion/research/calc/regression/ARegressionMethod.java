package com.actelion.research.calc.regression;

import com.actelion.research.calc.ProgressController;

import java.util.Properties;

/**
 * ARegressionMethod
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 28.11.18.
 */
abstract public class ARegressionMethod<T extends ParameterRegressionMethod> implements ICalculateModel {


    protected T parameterRegressionMethod;

    private ProgressController progressController;

    public void setParameterRegressionMethod(T parameterRegressionMethod) {
        this.parameterRegressionMethod = parameterRegressionMethod;
    }

    public T getParameter() {
        return parameterRegressionMethod;
    }

    public String getName(){
        return parameterRegressionMethod.getName();
    }

    public Properties getProperties() {return parameterRegressionMethod.getProperties();}

    public void decodeProperties2Parameter() {
        parameterRegressionMethod.decodeProperties2Parameter();
    }


    public ProgressController getProgressController() {
        return progressController;
    }

    public void setProgressController(ProgressController progressController) {
        this.progressController = progressController;
    }

}
