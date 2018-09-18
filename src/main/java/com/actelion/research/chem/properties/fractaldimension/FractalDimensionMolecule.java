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

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.properties.complexity.BitArray128;
import com.actelion.research.chem.properties.complexity.ExhaustiveFragmentsStatistics;
import com.actelion.research.chem.properties.complexity.ModelExhaustiveStatistics;
import com.actelion.research.chem.properties.complexity.ResultFragmentsStatistic;
import com.actelion.research.util.PointUtils;

import java.awt.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * FractalDimensionMolecule
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 28.08.18.
 */
public class FractalDimensionMolecule {

    private static final int MAX_THREADS_BOND_VECTOR_TO_IDCODE = 3;

    private ExhaustiveFragmentsStatistics exhaustiveFragmentsStatistics;

    public FractalDimensionMolecule(int totalCapacity, boolean elusive) {

        ExhaustiveFragmentsStatistics.setELUSIVE(elusive);

        int threadsBondVector2IdCode = Runtime.getRuntime().availableProcessors()-1;

        threadsBondVector2IdCode = Math.min(MAX_THREADS_BOND_VECTOR_TO_IDCODE, threadsBondVector2IdCode);

        if(threadsBondVector2IdCode==0){
            threadsBondVector2IdCode=1;
        }

        exhaustiveFragmentsStatistics = new ExhaustiveFragmentsStatistics(BitArray128.MAX_NUM_BITS, threadsBondVector2IdCode, totalCapacity);
    }

    public ResultFracDimCalc process(InputObjectFracDimCalc inputObjectFracDimCalc){

        ResultFracDimCalc resultFracDimCalc = new ResultFracDimCalc(inputObjectFracDimCalc);

        StereoMolecule mol = inputObjectFracDimCalc.getData();

        int bonds = inputObjectFracDimCalc.getData().getBonds();

        if(bonds<3){
            resultFracDimCalc.message = "Num bonds in molecule below limit of 3.";
            return resultFracDimCalc;
        } else if(bonds > exhaustiveFragmentsStatistics.getMaximumNumberBondsInMolecule()){
            resultFracDimCalc.message = "Num bonds in molecule above limit of " + exhaustiveFragmentsStatistics.getMaximumNumberBondsInMolecule() + ".";
            return resultFracDimCalc;
        }

        System.out.println("Process molecule " + inputObjectFracDimCalc.getId() + " with " + bonds + " bonds." );

        ResultFragmentsStatistic resultFragmentsStatistic = exhaustiveFragmentsStatistics.create(mol, bonds-1);

        List<ModelExhaustiveStatistics> liModelExhaustiveStatistics = resultFragmentsStatistic.getExhaustiveStatistics();


        List<Point> liFragBnds_NumUniqueFrags = new ArrayList<Point>();

        for (ModelExhaustiveStatistics modelExhaustiveStatistic : liModelExhaustiveStatistics) {
            liFragBnds_NumUniqueFrags.add(new Point(modelExhaustiveStatistic.getNumBondsInFragment(), modelExhaustiveStatistic.getUnique()));
        }

        Collections.sort(liFragBnds_NumUniqueFrags, PointUtils.getComparatorX());

        Point pBnds_MaxNumUniqueFrags = getMaxNumUniqueFrags(liFragBnds_NumUniqueFrags);

        int nBondsAtMaxNumFrags = pBnds_MaxNumUniqueFrags.x;
        int nMaxFrags = pBnds_MaxNumUniqueFrags.y;

        resultFracDimCalc.fractalDimension = Math.log10(nMaxFrags) / Math.log10(nBondsAtMaxNumFrags);

        resultFracDimCalc.bondsAtMaxFrag = nBondsAtMaxNumFrags;

        resultFracDimCalc.maxNumUniqueFrags = nMaxFrags;

        resultFracDimCalc.sumUniqueFrags = getSumUniqueFrags(liFragBnds_NumUniqueFrags);

        return resultFracDimCalc;
    }


    protected void finalizeThreads() throws Throwable {
        if(exhaustiveFragmentsStatistics!=null) {
            exhaustiveFragmentsStatistics.finalize();
        }
    }

    public static int getSumUniqueFrags(List<Point> liFragBnds_NumUniqueFrags) {

        int sum = 0;

        int indexMaxNumUniqueFrags = getIndexMaxNumUniqueFrags(liFragBnds_NumUniqueFrags);

        int end = indexMaxNumUniqueFrags+1;

        for (int i = 0; i < end; i++) {
            int nFragsUnique = liFragBnds_NumUniqueFrags.get(i).y;
            sum += nFragsUnique;
        }

        return sum;
    }

    public static Point getMaxNumUniqueFrags(List<Point> liFragBnds_NumUniqueFrags) {

        int indexMaxNumUniqueFrags = getIndexMaxNumUniqueFrags(liFragBnds_NumUniqueFrags);

        return liFragBnds_NumUniqueFrags.get(indexMaxNumUniqueFrags);
    }

    public static int getIndexMaxNumUniqueFrags(List<Point> liFragBnds_NumUniqueFrags){

        int indexMaxNumUniqueFrags = -1;

        int maxNumUniqueFrags = 0;

        for (int i = 0; i < liFragBnds_NumUniqueFrags.size(); i++) {
            Point pFragBnds_NumUniqueFrags = liFragBnds_NumUniqueFrags.get(i);

            int nUniqueFrags = pFragBnds_NumUniqueFrags.y;

            if(nUniqueFrags>maxNumUniqueFrags){
                maxNumUniqueFrags=nUniqueFrags;
                indexMaxNumUniqueFrags=i;
            }
        }

        return indexMaxNumUniqueFrags;
    }

}
