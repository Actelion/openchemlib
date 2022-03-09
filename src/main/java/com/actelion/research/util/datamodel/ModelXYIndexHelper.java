/*
 * Copyright (c) 1997 - 2022
 * Idorsia Pharmaceuticals Ltd.
 * Hegenheimermattweg 91
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
 * 3. Neither the name of the copyright holder nor the
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
 * @author Modest v. Korff
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
