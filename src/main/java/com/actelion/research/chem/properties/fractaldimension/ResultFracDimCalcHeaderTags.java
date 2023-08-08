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

import com.actelion.research.util.ConstantsDWAR;

import java.util.ArrayList;
import java.util.List;

/**
 * ResultFracDimCalcHeaderTags
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 21.07.23.
 */
public class ResultFracDimCalcHeaderTags {



    public static final String TAG_OUTCOME = "SuccessCalcFracDim";

    public static final String ATTR_SUCCESS = "true";
    public static final String ATTR_FAILURE = "false";

    private static final String TAG_SMILES = "SMILES";

    private static final String TAG_ID = "MoleculeId";

    private static final String TAG_SUM_UNIQUE_FRAGMENTS_CALC = "SumUniqueFragmentsCalculated";

    private static final String TAG_ATOM_COUNT = "AtomCountNonH";
    private static final String TAG_BOND_COUNT = "BondCountNonH";
    private static final String TAG_BONDS_AT_MAX_FRAGS_CALC = "BondNumberAtMaxNumFragCalculated";

    private static final String TAG_MAX_NUM_UNIQUE_FRAGMENTS_CALC = "MaxNumUniqueFragmentsCalculated";

    private static final String TAG_FRACTAL_DIM = "FractalDimension";
    private static final String TAG_MESSAGE = "Message";

    private static final String [] ARR_TAGS = {
            TAG_SMILES,
            TAG_ID,
            TAG_SUM_UNIQUE_FRAGMENTS_CALC,
            TAG_BONDS_AT_MAX_FRAGS_CALC,
            TAG_ATOM_COUNT,
            TAG_BOND_COUNT,
            TAG_MAX_NUM_UNIQUE_FRAGMENTS_CALC,
            TAG_FRACTAL_DIM,
            TAG_MESSAGE};


    public static final String SEP = "\t";

    private String tagColumnIdCode;

    public ResultFracDimCalcHeaderTags() {
        this(ConstantsDWAR.TAG_IDCODE2);
    }
    public ResultFracDimCalcHeaderTags(String tagColumnIdCode) {
        this.tagColumnIdCode = tagColumnIdCode;
    }

    public String getTagSmiles(){
        return TAG_SMILES+"_"+tagColumnIdCode;
    }

    public String getTagId(){
        return TAG_ID+"_"+tagColumnIdCode;
    }

    public String getTagSumUniqueFragmentsCalc(){
        return TAG_SUM_UNIQUE_FRAGMENTS_CALC+"_"+tagColumnIdCode;
    }

    public String getTagAtomCount(){
        return TAG_ATOM_COUNT+"_"+tagColumnIdCode;
    }

    public String getTagBondCount(){
        return TAG_ATOM_COUNT+"_"+tagColumnIdCode;
    }

    public String getTagBondsAtMaxFragsCalc(){
        return TAG_BONDS_AT_MAX_FRAGS_CALC+"_"+tagColumnIdCode;
    }

    public String getTagMaxNumUniqueFragmentsCalc(){
        return TAG_MAX_NUM_UNIQUE_FRAGMENTS_CALC+"_"+tagColumnIdCode;
    }

    public String getTagFractalDimension(){
        return TAG_FRACTAL_DIM+"_"+tagColumnIdCode;
    }

    public String getTagMessage(){
        return TAG_MESSAGE+"_"+tagColumnIdCode;
    }
    public String getTagOutcome(){
        return TAG_OUTCOME+"_"+tagColumnIdCode;
    }

    public String toStringHeader() {

        StringBuilder sb = new StringBuilder();

        sb.append(getTagSmiles());
        sb.append(SEP);
        sb.append(getTagId());
        sb.append(SEP);
        sb.append(getTagSumUniqueFragmentsCalc());
        sb.append(SEP);
        sb.append(getTagAtomCount());
        sb.append(SEP);
        sb.append(getTagBondCount());
        sb.append(SEP);
        sb.append(getTagBondsAtMaxFragsCalc());
        sb.append(SEP);
        sb.append(getTagMaxNumUniqueFragmentsCalc());
        sb.append(SEP);
        sb.append(getTagFractalDimension());
        sb.append(SEP);
        sb.append(getTagMessage());

        return sb.toString();
    }

    public List<String> getHeaderTags() {

        List<String> li = new ArrayList<>();

        li.add(getTagSumUniqueFragmentsCalc());
        li.add(getTagAtomCount());
        li.add(getTagBondCount());
        li.add(getTagBondsAtMaxFragsCalc());
        li.add(getTagMaxNumUniqueFragmentsCalc());
        li.add(getTagFractalDimension());
        li.add(getTagMessage());
        li.add(getTagOutcome());

        return li;
    }
    public List<String> getHeaderTags4SMILES() {

        List<String> li = new ArrayList<>(getHeaderTags());

        li.add(getTagSmiles());
        li.add(getTagId());

        return li;
    }
}
