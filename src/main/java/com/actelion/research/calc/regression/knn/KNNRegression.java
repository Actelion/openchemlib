package com.actelion.research.calc.regression.knn;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.MatrixFunctions;
import com.actelion.research.calc.SimilarityMulticore;
import com.actelion.research.calc.regression.ARegressionMethod;
import com.actelion.research.chem.descriptor.SimilarityCalculatorDoubleArray;
import com.actelion.research.util.datamodel.IdentifiedObject;
import com.actelion.research.util.datamodel.ModelXYIndex;

import java.util.ArrayList;
import java.util.List;

/**
 * KNNRegression
 *
 * kNN regression seems to be very inappropriate for regression problems. However, as base-line method it has a value.
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 24.11.18.
 * 20.09.2019 Index bug fixed.
 */
public class KNNRegression extends ARegressionMethod<ParameterKNN> {



    public static final int NEIGHBOURS = 3;

    private static final double TINY = 10e-6;

    private SimilarityCalculatorDoubleArray similarityCalculatorDoubleArray;
    private SimilarityMulticore<double []> similarityMulticore;

    private List<IdentifiedObject<double []>> liXTrain;

    private Matrix XTrain;
    private Matrix YTrain;

    public KNNRegression() {
        setParameterRegressionMethod(new ParameterKNN(NEIGHBOURS));
        similarityCalculatorDoubleArray = new SimilarityCalculatorDoubleArray();
        similarityMulticore = new SimilarityMulticore<>(similarityCalculatorDoubleArray);
    }

    public int getNeighbours() {
        return getParameter().getNeighbours();
    }

    public void setNeighbours(int neighbours) {
        getParameter().setNeighbours(neighbours);
    }

    public Matrix createModel(ModelXYIndex modelXYIndexTrain) {

        XTrain = modelXYIndexTrain.X;
        YTrain = modelXYIndexTrain.Y;

        liXTrain = MatrixFunctions.createIdentifiedObject(modelXYIndexTrain.X);
        similarityMulticore.run(liXTrain, liXTrain);
        Matrix maSimilarity = similarityMulticore.getSimilarityMatrix();
        Matrix maYHat = calculateYHat(maSimilarity, YTrain, getParameter().getNeighbours());
        return maYHat;
    }

    /**
     * Not thread save. Should not be called from several threads.
     * @param X
     * @return
     */
    @Override
    public Matrix calculateYHat(Matrix X){

        List<IdentifiedObject<double []>> liX = MatrixFunctions.createIdentifiedObject(X);

        similarityMulticore.run(liXTrain, liX);

        Matrix maSimMatrixTrainTest = similarityMulticore.getSimilarityMatrix();
// For debug
//        Matrix maSimMatrixTrainTest = MatrixFunctions.calculateSimilarityMatrixRowWise(XTrain, X);

        return calculateYHat(maSimMatrixTrainTest, YTrain, getParameter().getNeighbours());
    }

    /**
     * Thread save method. Can be called from different threads.
     * @param arrRow
     * @return
     */
    @Override
    public double calculateYHat(double[] arrRow) {

        int k = getParameter().getNeighbours();

        double [] arrSim = new double[k];
        int [] arrIndex = new int[k];

        for (int i = 0; i < liXTrain.size(); i++) {

            IdentifiedObject<double[]> idObj = liXTrain.get(i);

            double [] arrXTrain = idObj.getData();

            double sim = similarityCalculatorDoubleArray.getSimilarity(arrXTrain, arrRow);

            for (int j = k-1; j >= 0; j--) {
                if(sim>arrSim[j]){
                    arrSim[j] = sim;
                    arrIndex[j]=i;
                    break;
                }
            }
        }

        double [] arrNYTrain = new double[k];


        for (int i = 0; i < k; i++) {
            int indexTrain = arrIndex[i];
            arrNYTrain[i]= YTrain.get(indexTrain, 0);
        }

        double yHat = calculateYHat(arrSim, arrNYTrain);


        return yHat;
    }

    private static Matrix calculateYHat(Matrix maSimMatrixTrainTest, Matrix YTrain, int neighbours){

        int rowsTrain = maSimMatrixTrainTest.rows();
        int rowsTest = maSimMatrixTrainTest.cols();

        int colsY = YTrain.cols();

        // One entry in list for each test object
        // The arrays have the dimension k
        // The double array contains the k highest similarity values. The int array contains the corresponding object
        // indices of the train data.
        List<double[]> liSimilarityMostNSimilar = new ArrayList<>(rowsTest);
        List<int[]> liIndexMostNSimilar = new ArrayList<>(rowsTest);
        for (int i = 0; i < rowsTest; i++) {

            double [] arrSimilarityNSimilar= new double[neighbours];
            int [] arrNSimilarIndex= new int[neighbours];

            double minSim = 0;

            int indexMinSim = 0;

            for (int j = 0; j < rowsTrain; j++) {

                double v = maSimMatrixTrainTest.get(j,i);

                if(v>minSim){
                    arrSimilarityNSimilar[indexMinSim] = v;
                    arrNSimilarIndex[indexMinSim] = j;

                    minSim = Double.MAX_VALUE;
                    indexMinSim = -1;
                    for (int k = 0; k < arrSimilarityNSimilar.length; k++) {
                        if(arrSimilarityNSimilar[k] < minSim) {
                            minSim = arrSimilarityNSimilar[k];
                            indexMinSim = k;
                        }
                    }
                }
            }
            liSimilarityMostNSimilar.add(arrSimilarityNSimilar);
            liIndexMostNSimilar.add(arrNSimilarIndex);
        }

        double [][] arrYHat = new double[rowsTest][colsY];

        for (int h = 0; h < colsY; h++) {

            for (int i = 0; i < rowsTest; i++) {

                int[] arrNSimilarIndex = liIndexMostNSimilar.get(i);

                // Array with the y values from the train data which were most similar in X.
                double[] arrYN = new double[neighbours];

                for (int j = 0; j < arrNSimilarIndex.length; j++) {
                    arrYN[j] = YTrain.get(arrNSimilarIndex[j], h);
                }

                double yHat = calculateYHat(liSimilarityMostNSimilar.get(i), arrYN);

                arrYHat[i][h] = yHat;

            }
        }

        return new Matrix(arrYHat);
    }

    /**
     * Weights y linear with the similarity of the train data.
     * @param arrSimilarityNMostSimilar
     * @param arrNYTrain
     * @return
     */
    private static double calculateYHat(double [] arrSimilarityNMostSimilar, double [] arrNYTrain) {

        double sumSimilarity = 0;

        for (double v : arrSimilarityNMostSimilar) {
            sumSimilarity+=v;
        }

        double yHat = 0;

        if(sumSimilarity < TINY) {
            for (int i = 0; i < arrNYTrain.length; i++) {
                yHat += arrNYTrain[i];
            }

            yHat /= arrNYTrain.length;
        } else {
            for (int i = 0; i < arrNYTrain.length; i++) {

                double v = arrSimilarityNMostSimilar[i];
                yHat += arrNYTrain[i] * v;
            }

            yHat /= sumSimilarity;
        }

        return yHat;

    }

}
