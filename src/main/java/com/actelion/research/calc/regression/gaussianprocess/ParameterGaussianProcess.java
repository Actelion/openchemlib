package com.actelion.research.calc.regression.gaussianprocess;

import com.actelion.research.util.Formatter;
import com.actelion.research.calc.regression.ConstantsRegressionMethods;
import com.actelion.research.calc.regression.ParameterRegressionMethod;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * ParameterGaussianProcess
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 01.04.19.
 */
public class ParameterGaussianProcess extends ParameterRegressionMethod {

    public static final String TAG_LAMBDA="Lambda";

    public static final double LAMBDA = 0.005;

    private double lambda;


    public ParameterGaussianProcess() {
        super(ConstantsRegressionMethods.MODEL_GAUSSIAN_PROCESS);
        setLambda(LAMBDA);
    }


    @Override
    public int compareTo(ParameterRegressionMethod o) {

        int cmp = 0;

        ParameterGaussianProcess p = (ParameterGaussianProcess)o;

        return cmp;
    }


    public double getLambda() {
        return lambda;
    }

    public void setLambda(double lambda) {
        this.lambda = lambda;
        properties.put(TAG_LAMBDA, Double.toString(lambda));
    }

    @Override
    protected void decodeProperties2Parameter() {
        lambda = Double.parseDouble(properties.getProperty(TAG_LAMBDA));
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("ParameterGaussianProcess{");
        sb.append("lambda=").append(Formatter.format3(lambda));
        sb.append('}');
        return sb.toString();
    }

    public static List<String> getHeader(){

        List<String> li = ParameterRegressionMethod.getHeader();

        // li.add(TAG_KERNEL);
        li.add(TAG_LAMBDA);

        return li;
    }

    public static void main(String[] args) throws IOException {

        File dir = new File("/home/korffmo1/tmp/tmp00");

        File fiProp = new File(dir, "gaussianProcess.properties");

        ParameterGaussianProcess parameter = new ParameterGaussianProcess();

        parameter.lambda = 1.123456;

        parameter.write(fiProp);

        ParameterGaussianProcess parameterIn = new ParameterGaussianProcess();

        parameterIn.read(fiProp);

        System.out.println(parameterIn.toString());

    }

}
