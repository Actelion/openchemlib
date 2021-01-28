package com.actelion.research.calc.regression;

import com.actelion.research.calc.Matrix;

/**
 * ICalculateYHat
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 27.11.18.
 */
public interface ICalculateYHat {

    Matrix calculateYHat(Matrix X);

    double calculateYHat(double [] arrRow);

}
