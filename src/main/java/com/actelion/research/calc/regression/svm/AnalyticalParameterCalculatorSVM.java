/*
 * Copyright (c) 1997 - 2016
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the the copyright holder nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Modest v. Korff
 */

package com.actelion.research.calc.regression.svm;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.statistics.ModelStatisticsOverview;
import com.actelion.research.calc.statistics.StatisticsOverview;
import com.actelion.research.util.Formatter;
import com.actelion.research.util.datamodel.DoubleArray;
import com.actelion.research.util.datamodel.ModelXYIndex;
import com.actelion.research.calc.regression.ModelError;
import com.actelion.research.calc.regression.knn.KNNRegression;
import org.machinelearning.svm.libsvm.svm_parameter;

/**
 * AnalyticalParameterCalculatorSVM
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 01.11.19.
 */
public class AnalyticalParameterCalculatorSVM {

    public static ParameterSVM calculate(ModelXYIndex modelXYTrain){


        System.out.println("AnalyticalParameterCalculatorSVM");

        //
        // Parameter C
        //

        if(modelXYTrain.Y.cols() != 1) {
            throw new RuntimeException("Only single col allowed for Y!");
        }


        double rows = modelXYTrain.X.rows();

        System.out.println("Rows X train " + (int)rows);

        DoubleArray daY = new DoubleArray(modelXYTrain.Y.getColAsDouble(0));

        StatisticsOverview statisticsOverview = new StatisticsOverview(daY);

        ModelStatisticsOverview modelStatisticsOverviewY = statisticsOverview.evaluate();

        double sigmaY = modelStatisticsOverviewY.sdv;

        System.out.println("Sigma in y " + Formatter.format3(sigmaY));

        double avr = modelStatisticsOverviewY.avr;

        System.out.println("Average in y " + Formatter.format3(avr));

        //
        // Cherkassky, Vladimir, and Yunqian Ma.
        // "Practical selection of SVM parameters and noise estimation for SVM regression."
        // Neural networks 17.1 (2004): 113-126.
        // p5 equation 11
        double C = Math.max(Math.abs(avr-sigmaY*3), Math.abs(avr+sigmaY*3));

        System.out.println("C " + Formatter.format3(C));

        KNNRegression knnRegression = new KNNRegression();

        final int k = 3;

        knnRegression.setNeighbours(k);

        Matrix yHat = knnRegression.createModel(modelXYTrain);

        ModelError modelErrorKNN = ModelError.calculateError(modelXYTrain.Y, yHat);

        // Cherkassky, Vladimir, and Yunqian Ma.
        // "Practical selection of SVM parameters and noise estimation for SVM regression."
        // Neural networks 17.1 (2004): 113-126.
        // p18 equation 23
        double sigmaSquaredYHat = k/(k-1) * 1.0/rows * modelErrorKNN.errSumSquared;

        System.out.println("Sigma squared y hat " + Formatter.format3(sigmaSquaredYHat));

        // Cherkassky, Vladimir, and Yunqian Ma.
        // "Practical selection of SVM parameters and noise estimation for SVM regression."
        // Neural networks 17.1 (2004): 113-126.
        // p6 equation 14
        final double tau = 3;

        double epsilon = tau * Math.sqrt(sigmaSquaredYHat) * Math.sqrt(Math.log(rows)/rows);

        System.out.println("Epsilon " + Formatter.format3(epsilon));

        double gamma = 1.0 / modelXYTrain.X.cols();

        System.out.println("Gamma " + Formatter.format4(gamma));

        int kernelType = svm_parameter.RBF;

        svm_parameter svmParameter = SVMParameterHelper.regressionEpsilonSVR();

        svmParameter.kernel_type = kernelType;

        svmParameter.eps = epsilon;

        svmParameter.C = C;

        svmParameter.gamma = gamma;

        svmParameter.degree = 0;

        ParameterSVM parameterSVM = new ParameterSVM(svmParameter);

        return parameterSVM;

    }

}
