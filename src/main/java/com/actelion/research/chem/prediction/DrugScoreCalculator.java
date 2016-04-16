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

package com.actelion.research.chem.prediction;

/**
 * Created by rufenec on 26/10/15.
 */
public class DrugScoreCalculator
{
    public static double calculate(double mCLogP,double mSolubility,double mMolweight,double mDruglikeness,int[] toxRisks)
    {

        double cLogPScore = 1/(1+Math.exp(mCLogP-5));
        double solubilityScore = 1-1/(1+Math.exp(mSolubility+5));
        double molweightScore = 1/(1+Math.exp(0.012*mMolweight-6));
//        mTPSAScore = 1/(1+Math.exp(0.05*mTPSA-5));
        double drugLikenessScore = 1-1/(1+Math.exp(mDruglikeness));

        double drugScore = (0.5+cLogPScore/2)
                * (0.5+solubilityScore/2)
                * (0.5+molweightScore/2)
                * (0.5+drugLikenessScore/2);

        for (int i=0; toxRisks != null && i<toxRisks.length; i++) {
            if (toxRisks[i] == ToxicityPredictor.cLowRisk)
                drugScore *= 0.80;
            else if (toxRisks[i] == ToxicityPredictor.cHighRisk)
                drugScore *= 0.60;
        }
        return drugScore;


    }
}
