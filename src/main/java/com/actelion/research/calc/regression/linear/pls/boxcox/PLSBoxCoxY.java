package com.actelion.research.calc.regression.linear.pls.boxcox;

import com.actelion.research.calc.BoxCox;
import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.regression.linear.pls.PLSRegressionModelCalculator;
import com.actelion.research.util.datamodel.ModelXYColTags;
import com.actelion.research.util.datamodel.ModelXYIndex;

/**
 * PLSBoxCoxY
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 06.12.18.
 */
public class PLSBoxCoxY extends PLSRegressionModelCalculator {


    public static final int FACTORS = 15;

    public static final double LAMBDA = 0.6;

    private BoxCox boxCox;

    public PLSBoxCoxY() {
        setParameterRegressionMethod(new ParameterPLSBoxCox(FACTORS));

        ((ParameterPLSBoxCox)getParameter()).setLambda(LAMBDA);

        boxCox = new BoxCox(((ParameterPLSBoxCox)getParameter()).getLambda());
    }

    public PLSBoxCoxY(ParameterPLSBoxCox parameterPLSBoxCox) {
        setParameterRegressionMethod(parameterPLSBoxCox);

        ((ParameterPLSBoxCox)getParameter()).setLambda(LAMBDA);

        boxCox = new BoxCox(((ParameterPLSBoxCox)getParameter()).getLambda());
    }

    public double getLambda(){
        return ((ParameterPLSBoxCox)getParameter()).getLambda();
    }

    public void setLambda(double lambda){
        ((ParameterPLSBoxCox)getParameter()).setLambda(lambda);
    }

    @Override
    public Matrix createModel(ModelXYIndex modelXYIndexTrain) {

        super.setCenterData(true);

        ModelXYColTags modelXYColTagsPowerTrans = new ModelXYColTags();

        modelXYColTagsPowerTrans.X = modelXYIndexTrain.X;

        boxCox.setLambda(((ParameterPLSBoxCox)getParameter()).getLambda());

        modelXYColTagsPowerTrans.Y = BoxCox.transform(modelXYIndexTrain.Y, boxCox);

        Matrix yHat = super.createModel(modelXYColTagsPowerTrans);

        return yHat;
    }

    public Matrix calculateYHat(Matrix Xtest){

        Matrix YHatTest = super.calculateYHat(Xtest);

        Matrix YHatTestRe = BoxCox.reTransform(YHatTest, boxCox);

        return YHatTestRe;
    }

    @Override
    public double calculateYHat(double[] arrRow) {

        double yHatUntrans = super.calculateYHat(arrRow);

        double yHat = boxCox.inverse(yHatUntrans);

        return  yHat;

    }

    @Override
    public ParameterPLSBoxCox getParameter() {
        return (ParameterPLSBoxCox)super.getParameter();
    }
}
