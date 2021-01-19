package com.actelion.research.calc.regression.linear.pls.boxcox;

import com.actelion.research.calc.MatrixFunctions;
import com.actelion.research.calc.regression.ConstantsRegressionMethods;
import com.actelion.research.calc.regression.ParameterRegressionMethod;
import com.actelion.research.calc.regression.linear.pls.ParameterPLS;

import java.text.DecimalFormat;
import java.util.List;

/**
 * ParameterPLSBoxCox
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 06.12.18.
 */
public class ParameterPLSBoxCox extends ParameterPLS {

    public static final double LAMBDA= 0.4;

    public static final String TAG_LAMBDA="Lambda";

    private double lambda;


    public ParameterPLSBoxCox() {
        super(ConstantsRegressionMethods.MODEL_PLS_POWER, FACTORS);
        setLambda(LAMBDA);
    }

    public ParameterPLSBoxCox(int factors) {
        super(ConstantsRegressionMethods.MODEL_PLS_POWER, factors);
    }

    public ParameterPLSBoxCox(ParameterPLSBoxCox orig) {
        super(orig);
        copy(orig);
    }

    public void copy(ParameterPLSBoxCox orig){
        super.copy(orig);
        setLambda(orig.getLambda());
    }

    @Override
    public boolean equals(Object obj) {

        if(!(obj instanceof ParameterPLSBoxCox)){
            return false;
        }

        boolean eq = super.equals(obj);

        ParameterPLSBoxCox p = (ParameterPLSBoxCox)obj;

        if(!MatrixFunctions.equals(getLambda(), p.getLambda())){
            eq = false;
        }

        return eq;
    }


    public double getLambda() {
        return lambda;
    }

    public void setLambda(double lambda) {
        this.lambda = lambda;
        properties.put(TAG_LAMBDA, Double.toString(lambda));
    }

    @Override
    public int compareTo(ParameterRegressionMethod o) {

        int cmp = 0;

        ParameterPLSBoxCox parameterPLSBoxCox = (ParameterPLSBoxCox)o;

        if(getFactors()>parameterPLSBoxCox.getFactors()) {
            cmp=1;
        } else if(getFactors()<parameterPLSBoxCox.getFactors()) {
            cmp=-1;
        }

        if(cmp==0){

            if(lambda>parameterPLSBoxCox.lambda){
                cmp=1;
            }else if(lambda<parameterPLSBoxCox.lambda){
                cmp=-1;
            }

        }

        return cmp;
    }

    @Override
    protected void decodeProperties2Parameter() {
        factors = Integer.parseInt(properties.getProperty(TAG_FACTORS));
        lambda = Double.parseDouble(properties.getProperty(TAG_LAMBDA));
    }

    @Override
    public String toString() {

        DecimalFormat df = new DecimalFormat("0.0##");

        final StringBuilder sb = new StringBuilder("ParameterPLSBoxCox{");
        sb.append("name=").append(getName());
        sb.append(", factors=").append(getFactors());
        sb.append(", lambda=").append(df.format(lambda));
        sb.append('}');
        return sb.toString();
    }

    public static List<String> getHeader(){

        List<String> li = ParameterRegressionMethod.getHeader();

        li.add(TAG_FACTORS);
        li.add(TAG_LAMBDA);

        return li;
    }

}
