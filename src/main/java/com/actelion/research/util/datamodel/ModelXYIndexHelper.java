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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * ModelXYIndexHelper
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 08.11.19.
 */
public class ModelXYIndexHelper {

    /**
     * Sorts the data for the first column in Y and splits them.
     * @param m
     * @param fractionFirst
     * @return
     */
    public static ModelXYIndex [] splitLowHighY(ModelXYIndex m, double fractionFirst){

        List<XY> li = m.getAsList();

        Collections.sort(li, XY.getComparatorY(0));

        int nFirst = (int)(li.size() * fractionFirst);

        List<XY> liFirst = new ArrayList<>();
        for (int i = 0; i < nFirst; i++) {
            liFirst.add(li.get(i));
        }

        List<XY> liSec = new ArrayList<>();
        for (int i = nFirst; i < li.size(); i++) {
            liSec.add(li.get(i));
        }

        ModelXYIndex XYFirst = new ModelXYIndex(new ModelXY(liFirst));
        ModelXYIndex XYSecond = new ModelXYIndex(new ModelXY(liSec));

        ModelXYIndex [] arr = {XYFirst, XYSecond};

        return arr;
    }

    public static ModelXYIndex [] splitShuffled(ModelXYIndex m, double fractionFirst){

        List<XY> li = m.getAsList();

        Collections.shuffle(li);

        int nFirst = (int)(li.size() * fractionFirst);

        List<XY> liFirst = new ArrayList<>();
        for (int i = 0; i < nFirst; i++) {
            liFirst.add(li.get(i));
        }

        List<XY> liSec = new ArrayList<>();
        for (int i = nFirst; i < li.size(); i++) {
            liSec.add(li.get(i));
        }

        ModelXYIndex XYFirst = new ModelXYIndex(new ModelXY(liFirst));
        ModelXYIndex XYSecond = new ModelXYIndex(new ModelXY(liSec));

        ModelXYIndex [] arr = {XYFirst, XYSecond};

        return arr;
    }

}
