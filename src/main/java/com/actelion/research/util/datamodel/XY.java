/*
 * Copyright (c) 2019.
 * Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 *  This file is part of DataWarrior.
 *
 *  DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.
 *
 *  DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License along with DataWarrior.
 *  If not, see http://www.gnu.org/licenses/.
 *
 *  @author Modest v. Korff
 *
 */

package com.actelion.research.util.datamodel;

import java.util.Comparator;

/**
 * XY
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 05.11.19.
 */
public class XY {

    double [] x;
    double [] y;

    public XY(double[] x, double[] y) {
        this.x = x;
        this.y = y;
    }


    public static Comparator<XY> getComparatorY(int colY){
        return new ComparatorColY(colY);
    }

    private static class ComparatorColY implements Comparator<XY> {

        private int colY;

        public ComparatorColY(int colY) {
            this.colY = colY;
        }

        @Override
        public int compare(XY o1, XY o2) {

            int cmp=0;
            if(o1.y[colY]> o2.y[colY]){
                cmp = 1;
            }else if(o1.y[colY] < o2.y[colY]){
                cmp = -1;
            }

            return cmp;
        }
    }
}
