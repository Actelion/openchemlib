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

package com.actelion.research.chem.properties.fractaldimension;

import com.actelion.research.chem.ExtendedMoleculeFunctions;
import com.actelion.research.util.Formatter;

import java.util.ArrayList;
import java.util.List;

/**
 * ResultFracDimCalc
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 28.08.18.
 */
public class ResultFracDimCalc extends InputObjectFracDimCalc {



    int idMolecule;

    double fractalDimension;

    int atomCount;
    int bondCount;
    int bondsAtMaxFrag;

    int maxNumUniqueFrags;

    int sumUniqueFrags;

    String message;

    public ResultFracDimCalc(InputObjectFracDimCalc inputObjectFracDimCalc) {
        super(inputObjectFracDimCalc);

        idMolecule = -1;

        fractalDimension = Double.NaN;

        atomCount = ExtendedMoleculeFunctions.getNumNonHydrogenAtoms(inputObjectFracDimCalc.getData());
        bondCount = ExtendedMoleculeFunctions.getNumBondsNoHydrogen(inputObjectFracDimCalc.getData());

        bondsAtMaxFrag = -1;

        maxNumUniqueFrags = -1;

        sumUniqueFrags = -1;

        message = "";
    }

    public double getFractalDimension() {
        return fractalDimension;
    }

    public int getBondsAtMaxFrag() {
        return bondsAtMaxFrag;
    }

    public int getBondCount() {
        return bondCount;
    }

    public int getAtomCount() {
        return atomCount;
    }

    public int getMaxNumUniqueFrags() {
        return maxNumUniqueFrags;
    }

    public int getSumUniqueFrags() {
        return sumUniqueFrags;
    }

    public String getMessage() {
        return message;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();

        sb.append(getSmiles());
        sb.append(ResultFracDimCalcHeaderTags.SEP);
        sb.append(getId());
        sb.append(ResultFracDimCalcHeaderTags.SEP);
        sb.append(getSumUniqueFrags());
        sb.append(ResultFracDimCalcHeaderTags.SEP);
        sb.append(getAtomCount());
        sb.append(ResultFracDimCalcHeaderTags.SEP);
        sb.append(getBondCount());
        sb.append(ResultFracDimCalcHeaderTags.SEP);
        sb.append(getBondsAtMaxFrag());
        sb.append(ResultFracDimCalcHeaderTags.SEP);
        sb.append(getMaxNumUniqueFrags());
        sb.append(ResultFracDimCalcHeaderTags.SEP);
        sb.append(Formatter.format3(getFractalDimension()));
        sb.append(ResultFracDimCalcHeaderTags.SEP);
        sb.append(getMessage());

        return sb.toString();
    }

}
