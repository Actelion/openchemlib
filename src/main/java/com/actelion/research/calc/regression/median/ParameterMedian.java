package com.actelion.research.calc.regression.median;

import com.actelion.research.calc.regression.ConstantsRegressionMethods;
import com.actelion.research.calc.regression.ParameterRegressionMethod;

import java.util.List;

/**
 * ParameterMedian
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 13.02.19.
 */
public class ParameterMedian extends ParameterRegressionMethod {

    public ParameterMedian() {
        super(ConstantsRegressionMethods.MODEL_MEDIAN);
    }


    @Override
    public int compareTo(ParameterRegressionMethod o) {

        int cmp = 0;

        return cmp;
    }

    @Override
    protected void decodeProperties2Parameter() {
        // Nothing to do.
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("ParameterMedian{");
        sb.append("name=").append(getName());
        sb.append('}');
        return sb.toString();
    }

    public static List<String> getHeader(){

        List<String> li = ParameterRegressionMethod.getHeader();

        return li;
    }

}
