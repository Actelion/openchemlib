package com.actelion.research.util.datamodel;

/**
 * Modest v. Korff
 * Idorsia Pharmaceuticals Ltd.
 * 08.02.2021 Start implementation
 **/
public class XYIndex extends XY {

    int index;

    public XYIndex(double[] x, double[] y) {
        super(x, y);
    }

    public XYIndex(double[] x, double[] y, int index) {
        super(x, y);
        this.index = index;
    }
}
