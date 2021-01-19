package com.actelion.research.calc.regression.svm;

import com.actelion.research.calc.Matrix;
import org.machinelearning.svm.libsvm.svm_node;

/**
 * Matrix2SVMNodeConverter
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 27.11.18.
 */
public class Matrix2SVMNodeConverter {


    public static svm_node[][] convert(Matrix A) {

        svm_node[][] arrSVMNode = new svm_node[A.rows()][A.cols()];

        int rows = A.rows();
        int cols = A.cols();

        for (int i = 0; i < rows; i++) {

            for (int j = 0; j < cols; j++) {

                svm_node node = new svm_node();

                node.index = j;
                node.value = A.get(i,j);

                arrSVMNode[i][j]=node;
            }
        }

        return arrSVMNode;
    }

    public static svm_node[][] convertSingleRow(double [] a) {

        int cols = a.length;

        svm_node[] arrSVMNode = new svm_node[cols];

        for (int i = 0; i < cols; i++) {

            svm_node node = new svm_node();

            node.index = i;
            node.value = a[i];
            arrSVMNode[i]=node;
        }

        svm_node[][] arrSVMNodeRow = new svm_node[1][cols];

        arrSVMNodeRow[0]=arrSVMNode;

        return arrSVMNodeRow;
    }

}
