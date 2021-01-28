package com.actelion.research.calc.regression;

import com.actelion.research.calc.Matrix;
import com.actelion.research.util.datamodel.ModelXYIndex;

/**
 * ICalculateModel
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 28.11.18.
 */
public interface ICalculateModel extends ICalculateYHat {

    Matrix createModel(ModelXYIndex modelXYIndexTrain);

    ParameterRegressionMethod getParameter();

}
