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
*/

package com.actelion.research.calc;


public class CorrelationCalculator {
    public static final String[] TYPE_LONG_NAME = { "Bravais-Pearson (linear correlation)", "Spearman (correlation of ranks)" };
    public static final String[] TYPE_NAME = { "Bravais-Pearson", "Spearman" };
    public static final String[] TYPE_CODE = { "bravais-pearson", "spearman" };
    public static final int TYPE_NONE = -1;
    public static final int TYPE_BRAVAIS_PEARSON = 0;
    public static final int TYPE_SPEARMAN = 1;

    private int mValueCount;

    private int[][] mValueCountMatrix;

    /**
     * Calculates the correlation coefficient between two columns of data.
     * Use the TYPE_BRAVAIS_PEARSON for normal distributed data and the
     * more robust TYPE_SPEARMAN if the data is not normal distributed.
     * If type==TYPE_BRAVAIS_PEARSON and one of a row's values is Double.NaN,
     * then the row is skipped.
     * If less than two valid rows are found or if both columns have a
     * different number of values, than Double.NaN is returned.
     * @param column1
     * @param column2
     * @param correlationType
     * @return
     */
    public double calculateCorrelation(INumericalDataColumn column1,
                                              INumericalDataColumn column2,
                                              int correlationType) {
        int valueCount = column1.getValueCount();
        if (valueCount != column2.getValueCount())
            return Double.NaN;

        double r = Double.NaN;

        if (correlationType == TYPE_BRAVAIS_PEARSON) {
            // http://de.wikibooks.org/wiki/Mathematik:_Statistik:_Korrelationsanalyse
            mValueCount = 0;
            double xMean = 0;
            double yMean = 0;
            for (int i=0; i<valueCount; i++) {
                double x = column1.getValueAt(i);
                double y = column2.getValueAt(i);
                if (!Double.isNaN(x) && !Double.isNaN(y)) {
                    xMean += x;
                    yMean += y;
                    mValueCount++;
                    }
                }

            if (mValueCount < 2)
                return Double.NaN;

            xMean /= mValueCount;
            yMean /= mValueCount;

            double sumdxdx = 0;
            double sumdxdy = 0;
            double sumdydy = 0;
            for (int i=0; i<valueCount; i++) {
                double x = column1.getValueAt(i);
                double y = column2.getValueAt(i);
                if (!Double.isNaN(x) && !Double.isNaN(y)) {
                    double dx = x - xMean;
                    double dy = y - yMean;
                    sumdxdx += dx*dx;
                    sumdxdy += dx*dy;
                    sumdydy += dy*dy;
                    }
                }
            r = sumdxdy / Math.sqrt(sumdxdx * sumdydy);
            }
        else if (correlationType == TYPE_SPEARMAN) {
            if (valueCount < 2)
                return Double.NaN;

            mValueCount = 0;
            double[] xValue = new double[valueCount];
            double[] yValue = new double[valueCount];
            for (int i=0; i<valueCount; i++) {
                xValue[mValueCount] = column1.getValueAt(i);
                yValue[mValueCount] = column2.getValueAt(i);
                if (!Double.isNaN(xValue[mValueCount]) && !Double.isNaN(yValue[mValueCount]))
                    mValueCount++;
                }
            for (int i=mValueCount; i<valueCount; i++) {
                xValue[i] = Double.NaN;
                yValue[i] = Double.NaN;
                }
            java.util.Arrays.sort(xValue);
            java.util.Arrays.sort(yValue);

            double sumdxdx = 0;
            double sumdxdy = 0;
            double sumdydy = 0;
            double mean = (double)(mValueCount + 1) / 2;
            for (int i=0; i<valueCount; i++) {
            	if (!Double.isNaN(column1.getValueAt(i)) && !Double.isNaN(column2.getValueAt(i))) {
	                double xPosition = getPosition(xValue, column1.getValueAt(i));
	                double yPosition = getPosition(yValue, column2.getValueAt(i));
	                double dx = xPosition - mean;
	                double dy = yPosition - mean;
	                sumdxdx += dx*dx;
	                sumdxdy += dx*dy;
	                sumdydy += dy*dy;
            		}
                }
            r = sumdxdy / Math.sqrt(sumdxdx * sumdydy);
            }

        return r;
        }

    /**
     * @return number of non-null values used for calculation
     */
    public int getValueCount() {
        return mValueCount;
        }

    public int[][] getValueCountMatrix() {
        return mValueCountMatrix;
    }

    /**
     * Calculates a half correlation matrix of all passed numerical columns
     * @param numericalColumn
     * @param type
     * @return half matrix with matrix.length=numericalColumn.length and matrix[i].length=i
     */
    public double[][] calculateMatrix(final INumericalDataColumn[] numericalColumn, int type) {
        double[][] matrix = new double[numericalColumn.length][];
        mValueCountMatrix = new int[numericalColumn.length][];
        for (int i=1; i<numericalColumn.length; i++) {
            matrix[i] = new double[i];
            mValueCountMatrix[i] = new int[i];
            for (int j=0; j<i; j++) {
                matrix[i][j] = calculateCorrelation(numericalColumn[i], numericalColumn[j], type);
                mValueCountMatrix[i][j] = mValueCount;
                }
            }
        return matrix;
        }

    private double getPosition(double[] array, double value) {
        int position = java.util.Arrays.binarySearch(array, value);
        int position1 = position;
        while (position1 > 0 && array[position1-1] == value)
            position1--;
        int position2 = position;
        while (position2 < array.length-1 && array[position2+1] == value)
            position2++;
        return ((double)(position1+position2))/2 + 1;
        }
    }
