package com.actelion.research.calc.histogram;

import com.actelion.research.calc.Matrix;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;

/**
 * MatrixBasedHistogram
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 21.12.18.
 */
public class MatrixBasedHistogram {

    public static boolean VERBOSE = false;

    public static Histogram getFromQuadraticMatrix(Matrix ma, double min, double max, int bins){

        int r = ma.rows();

        int n = ((r * r) - r)/2;

        double [] a = new double[n];

        int c = 0;
        for (int i = 0; i < r; i++) {

            for (int j = i+1; j < r; j++) {
                a[c++] = ma.get(i,j);
            }
        }

        return new Histogram(a, min, max, bins);

    }

    /**
     *
     * @param ma matrix values are written into the histogram
     * @param maHistogramBins matrix with two rows row 1: lower bins; row 2: upper
     * bins.
     * @return matrix with 3 rows. Row 0 lower bins, row 1: upper bins, row 2
     * frequency.
     */
    public static Matrix getHistogram(Matrix ma, Matrix maHistogramBins) {

        Matrix maHistogram = new Matrix(maHistogramBins);

        maHistogram.resize(3, maHistogramBins.getColDim());

        for (int i = 0; i < ma.getRowDim(); i++) {
            for (int j = 0; j < ma.getColDim(); j++) {
                double v = ma.get(i, j);
                placeValueInHistogramBin(v, maHistogram);
            }
        }

        return maHistogram;
    }

    public static Matrix getHistogram(double [] arrValues, Matrix maHistogramBins) {
        Matrix maHistogram = new Matrix(maHistogramBins);

        maHistogram.resize(3, maHistogramBins.getColDim());

        for (int i = 0; i < arrValues.length; i++) {
            placeValueInHistogramBin(arrValues[i], maHistogram);
        }

        return maHistogram;
    }

    private static void placeValueInHistogramBin(double v, Matrix maHist){

        if(Double.isNaN(v)) {
            if(VERBOSE)
                System.err.println("Warning NaN in placeValueInHistogramBin(...).");
            return;
        }

        if (v < maHist.get(0, 0) || v > maHist.get(1, maHist.cols()-1))
            return;

        int maxloops = maHist.cols();

        boolean bEnd = false;

        int pos = maHist.cols() / 2;

        int posLow=0;

        int posUp = maHist.cols()-1;

        int cc=0;
        while(!bEnd){

            if (v >= maHist.get(0, pos) && v < maHist.get(1, pos)) {
                maHist.increase(2, pos, 1);
                bEnd=true;
            } else if (v < maHist.get(0, pos)) {
                posUp = pos;
                pos = posLow + (pos-posLow)/2;
            } else if (v >= maHist.get(1, pos)) {
                posLow = pos;
                pos = pos + (int)(((double)posUp-pos)/2+0.5);
            }

            if(cc==maxloops){

                String msg = "Fitting bin for value " + v + " not found!\n" +
                        MatrixBasedHistogram.histogram2String(maHist, 2, 8);

                throw new RuntimeException("Fitting bin for value " + v + " not found");
            }
            cc++;

        }

    }

    public static Matrix getHistogram(Matrix maHistogramBins) {
        Matrix maHistogram = new Matrix(maHistogramBins);
        maHistogram.resize(3, maHistogramBins.getColDim());
        return maHistogram;
    }

    /**
     * Gets the upper and lower limit of the most occupied bin.
     * @param maHistogram
     * @param radius so many bins are taken from the left and the right.
     * @return arr[0]: lower limit, arr[1]: upper limit
     */
    public static double [] getBordersMostFreqOccBin(Matrix maHistogram, int radius) {
        double [] arr = new double[2];

        double maxFreq = 0;
        int index=-1;
        for (int i = radius; i < maHistogram.cols() - radius; i++) {

            int sumFreq = 0;
            for (int j = -radius; j < radius+1; j++) {
                sumFreq += maHistogram.get(2, i+j);
            }

            if(sumFreq>maxFreq){
                maxFreq=sumFreq;
                index=i;
            }
        }

        int indexLower = Math.max(0, index-radius);

        int indexUpper = Math.min(maHistogram.cols()-1, index+radius);

        arr[0]=maHistogram.get(0, indexLower);

        arr[1]=maHistogram.get(1, indexUpper);

        return arr;
    }

    public static Matrix getHistogram(float [] ma, Matrix maHistogramBins) {
        Matrix maHistogram = new Matrix(maHistogramBins);
        maHistogram.resize(3, maHistogramBins.getColDim());
        for (int ii = 0; ii < ma.length; ii++) {
            for (int kk = 0; kk < maHistogramBins.getColDim(); kk++) {
                if (ma[ii] >= maHistogram.get(0, kk) &&
                        ma[ii] < maHistogram.get(1, kk)) {
                    int iCounter = (int) maHistogram.get(2, kk);
                    iCounter++;
                    maHistogram.set(2, kk, iCounter);
                }
            }
        }

        return maHistogram;
    }

    public static Matrix getHistogram(float [][] ma, int col, Matrix maHistogramBins) {
        Matrix maHistogram = new Matrix(maHistogramBins);
        maHistogram.resize(3, maHistogramBins.getColDim());
        for (int ii = 0; ii < ma.length; ii++) {
            for (int kk = 0; kk < maHistogramBins.getColDim(); kk++) {
                if (ma[ii][col] >= maHistogram.get(0, kk) &&
                        ma[ii][col] < maHistogram.get(1, kk)) {
                    int iCounter = (int) maHistogram.get(2, kk);
                    iCounter++;
                    maHistogram.set(2, kk, iCounter);
                }
            }
        }

        return maHistogram;
    }

    public static Matrix getHistogram(Matrix ma, int numBins) {
        double min = ma.getMin() - ma.getMin() * ConstantsHistogram.TINY_FACTOR;
        double max = ma.getMax() + ma.getMax() * ConstantsHistogram.TINY_FACTOR;

        Matrix maBins = getHistogramBins(min,max, numBins);

        Matrix maHist = getHistogram(ma, maBins);

        return maHist;
    }

    /**
     *
     * @param maHist histogram
     * @return the lower limit of the first occupied bin in the histogram.
     */
    public static double getMinOccBin(Matrix maHist) {
        double min = 0;
        for (int i = 0; i < maHist.getColDim(); i++) {
            if(maHist.get(2,i) > 0) {
                min = maHist.get(0,i);
                break;
            }
        }
        return min;
    }

    /**
     *
     * @param maHist histogram
     * @return the higher limit of the last occupied bin in the histogram.
     */
    public static double getMaxOccBin(Matrix maHist) {
        double max = 0;
        for (int i = maHist.getColDim() - 1; i >= 0; i--) {
            if(maHist.get(2,i) > 0) {
                max = maHist.get(1,i);
                break;
            }
        }
        return max;
    }

    /**
     *
     * @param dMin smallest value to put into the histogram
     * @param dMax maximum value to be considered.
     * @param iNumBins number of bins, between min and max.
     * @return matrix with two rows, the lower and upper bins.
     * The last bin is a little bit bigger than the other bins. So the highest
     * value fits into it.
     */
    public static Matrix getHistogramBins(double dMin, double dMax, int iNumBins) {
        Matrix maHistogramBins = new Matrix(2, iNumBins);

        double dDelta = dMax - dMin;
        double dBinWidth = dDelta / iNumBins;

        double dIncrementLast = dBinWidth * 0.0001;

        double dLow = dMin;

        int iCols = maHistogramBins.getColDim();
        for (int ii = 0; ii < iCols; ii++) {
            maHistogramBins.set(0, ii, dLow);
            double dUp = dLow + dBinWidth;

            if (ii == (iCols - 1))
                dUp += dIncrementLast;

            maHistogramBins.set(1, ii, dUp);
            dLow = dUp;
        }

        return maHistogramBins;
    }

    public static void writeHistogram(String sFile, Matrix hist, boolean bApppend, int digits, int totalWidth) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
                sFile), bApppend));
        DecimalFormat dfBins = Matrix.format(digits);
        String sVal = "";
        for (int ii = 0; ii < 2; ii++) {
            sVal = "";
            for (int jj = 0; jj < hist.getColDim(); jj++) {
                sVal += Matrix.format(hist.get(ii,jj), dfBins, totalWidth) + Matrix.OUT_SEPARATOR_COL;

            }
            sVal += Matrix.OUT_SEPARATOR_ROW;
            writer.write(sVal);
        }

        DecimalFormat dfFreq = new DecimalFormat();
        sVal = "";
        for (int jj = 0; jj < hist.getColDim(); jj++) {
            sVal += Matrix.format(hist.get(2, jj), dfFreq, totalWidth) + Matrix.OUT_SEPARATOR_COL;

        }
        writer.write(sVal + "\n");

        writer.flush();
        writer.close();

    }

    /**
     *
     * @param hist
     * @param digits
     * @param totalWidth for one number.
     * @return
     */
    public static String histogram2String(Matrix hist, int digits, int totalWidth) {

        StringBuilder sb = new StringBuilder();


        DecimalFormat dfBins = Matrix.format(digits);
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < hist.getColDim(); j++) {
                sb.append(Matrix.format(hist.get(i,j), dfBins, totalWidth) + Matrix.OUT_SEPARATOR_COL);

            }
            sb.append(Matrix.OUT_SEPARATOR_ROW);
        }

        DecimalFormat dfFreq = new DecimalFormat();
        for (int i = 0; i < hist.getColDim(); i++) {
            sb.append(Matrix.format(hist.get(2, i), dfFreq, totalWidth) + Matrix.OUT_SEPARATOR_COL);

        }
        sb.append("\n");

        return sb.toString();
    }
}
