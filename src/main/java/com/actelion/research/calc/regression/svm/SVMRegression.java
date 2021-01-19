package com.actelion.research.calc.regression.svm;

import com.actelion.research.calc.CorrelationCalculator;
import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.ProgressController;
import com.actelion.research.util.Formatter;
import com.actelion.research.util.datamodel.DoubleArray;
import com.actelion.research.util.datamodel.ModelXYIndex;
import com.actelion.research.calc.regression.ARegressionMethod;
import org.machinelearning.svm.libsvm.svm_model;
import org.machinelearning.svm.libsvm.svm_problem;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * SVMRegression
 *
 * Calculates the SVM regression for one Y column.
 *
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 27.11.18.
 *
 * http://www.svms.org/parameters/
 */
public class SVMRegression extends ARegressionMethod<ParameterSVM> implements Comparable<SVMRegression> {


    public static final int LIMIT_ROWS_ANALYTICAL = 1000;

    public static boolean VERBOSE = false;

    public static final double NU = 0.1;

    private svm_model modelSVM;

    public SVMRegression() {
        setParameterRegressionMethod(new ParameterSVM(SVMParameterHelper.regressionEpsilonSVR()));
    }

    public SVMRegression(ParameterSVM parameterSVM) {

        setParameterRegressionMethod(parameterSVM);
    }

    public Matrix createModel(ModelXYIndex modelXYIndex) {

        if(parameterRegressionMethod.getSvmParameter().degree==SVMParameterHelper.DEGREE_ANALYTICALLY_PARAMETER_CALC) {

            ModelXYIndex modelXYIndexAnalytical = null;

            int rows = modelXYIndex.X.rows();

            if(rows>LIMIT_ROWS_ANALYTICAL) {

                List<Integer> liIndex = new ArrayList<>(rows);

                for (int i = 0; i < rows; i++) {
                    liIndex.add(i);
                }

                Collections.shuffle(liIndex);

                List<Integer> liIndexSub = new ArrayList<>(liIndex.subList(0, LIMIT_ROWS_ANALYTICAL));

                modelXYIndexAnalytical = modelXYIndex.sub(liIndexSub);
            } else{
                modelXYIndexAnalytical = modelXYIndex;
            }

            ParameterSVM parameterSVM = AnalyticalParameterCalculatorSVM.calculate(modelXYIndexAnalytical);

            setParameterRegressionMethod(parameterSVM);

        }

        int rows = modelXYIndex.X.rows();
        int colsY = modelXYIndex.Y.cols();

        if(colsY!=1){
            throw new RuntimeException("Only one y column allowed!");
        }

        Matrix YHat = new Matrix(rows, colsY);

        svm_problem svm_problem = new svm_problem();

        svm_problem.l = rows;

        svm_problem.x = Matrix2SVMNodeConverter.convert(modelXYIndex.X);

        if(getParameter().getGamma()<=0) {
            getParameter().setGamma(1.0 / modelXYIndex.X.cols());
        }

        boolean failed = false;

        try {
            double [] arrY = modelXYIndex.Y.getColAsDouble(0);

            DoubleArray daY = new DoubleArray(arrY);

            svm_problem.y = arrY;

            String error_msg = svm.svm_check_parameter(svm_problem, getParameter().getSvmParameter());

            if(error_msg != null) {
                System.err.print("SVMRegressionMultiY svm_check_parameter error: "+error_msg+"\n");
                failed = true;
            }

            ProgressController pg = getProgressController();

            modelSVM = svm.svm_train(svm_problem, getParameter().getSvmParameter(), pg);

            // System.out.println("SVM model created");

            if(pg!=null) {
                pg.startProgress("Calculate train data fit", 0, rows);
            }

            DoubleArray daYHat = new DoubleArray(rows);
            for (int j = 0; j < rows; j++) {
                double yHat = svm.svm_predict(modelSVM, svm_problem.x[j]);
                daYHat.add(yHat);
                YHat.set(j,0,yHat);


                if(pg!=null) {
                    pg.updateProgress(j);
                    if (pg.threadMustDie()) {
                        failed = true;
                        break;
                    }
                }
            }

            if(VERBOSE) {
                double r = new CorrelationCalculator().calculateCorrelation(daY, daYHat, CorrelationCalculator.TYPE_BRAVAIS_PEARSON);
                double r2 = r * r;
                System.out.println("SVMRegressionMultiY model r2 " + Formatter.format4(r2) + ".");
            }
        } catch (Exception e) {
            e.printStackTrace();
            if(VERBOSE)
                System.err.println("SVMRegressionMultiY break.");
            failed = true;
        }



        if(failed) {
            YHat=null;
        }

        return YHat;
    }

    @Override
    public Matrix calculateYHat(Matrix X) {

        int rows = X.rows();


        Matrix YHat = new Matrix(rows, 1);

        svm_problem svm_problem = new svm_problem();

        svm_problem.l = rows;

        svm_problem.x = Matrix2SVMNodeConverter.convert(X);

        for (int j = 0; j < rows; j++) {
            double yHat = svm.svm_predict(modelSVM, svm_problem.x[j]);
            YHat.set(j,0,yHat);
        }

        return YHat;
    }

    @Override
    public double calculateYHat(double[] arrRow) {

        double yHat;

        synchronized (this) {

            svm_problem svm_problem = new svm_problem();

            svm_problem.l = 1;

            svm_problem.x = Matrix2SVMNodeConverter.convertSingleRow(arrRow);

            yHat = svm.svm_predict(modelSVM, svm_problem.x[0]);
        }

        return yHat;
    }

    public void setNu(double nu){
        getParameter().setNu(nu);
    }
    public void setGamma(double gamma){
        getParameter().setGamma(gamma);
    }

    @Override
    public int compareTo(SVMRegression o) {
        return getParameter().compareTo(o.getParameter());
    }
}
