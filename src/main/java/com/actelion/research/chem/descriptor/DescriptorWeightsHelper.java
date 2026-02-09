package com.actelion.research.chem.descriptor;

import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.descriptor.flexophore.generator.CreatorMolDistHistViz;
import com.actelion.research.chem.descriptor.flexophore.redgraph.SubGraphIndices;

import java.util.Arrays;
import java.util.List;

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
 */
public class DescriptorWeightsHelper {

    public static final String TAG_WEIGHTS_ATOMS_FLEXOPHORE = "Atom Weights Flexophore";
    public static final String TAG_WEIGHTS_ATOMS_PHESA = "Atom Weights PheSA";


    public static final int LABEL_WEIGHT_LOW = 0;
    public static final int LABEL_WEIGHT_NORMAL = 1;
    public static final int LABEL_WEIGHT_MANDATORY = 2;


    private CreatorMolDistHistViz creatorMolDistHistViz;



    public DescriptorWeightsHelper() {
        this.creatorMolDistHistViz = new CreatorMolDistHistViz();

    }

    /**
     * The labels are used to calculate the weight values for Flexophore and PheSA.
     * The labels can be given by the users in the
     * @param molecule3D
     * @return
     */
    public int [] calcCategory4WeightsFlexophore(Molecule3D molecule3D){
        List<SubGraphIndices> liSubGraphIndices = creatorMolDistHistViz.getSubGraphIndices(molecule3D);
        int [] category = calcCategory4Weights(liSubGraphIndices, molecule3D);
        return category;
    }

    /**
     * - End standing pp points are set mandatory.
     * - Charged pp points are set mandatory.
     * @param liSubGraphIndices
     * @param molecule3D
     * @return array with dimension molecule3D.getAtoms().
     */
    public static int [] calcCategory4Weights(List<SubGraphIndices> liSubGraphIndices, Molecule3D molecule3D){

        int [] weights = getBasisWeightCategories(molecule3D.getAtoms());

        for (int i = 0; i < liSubGraphIndices.size(); i++) {
            SubGraphIndices sgi = liSubGraphIndices.get(i);
            int weight = DescriptorWeightsHelper.LABEL_WEIGHT_NORMAL;
            if (SubGraphIndices.isLinker(molecule3D, liSubGraphIndices, i)) {
                if(!SubGraphIndices.isOnlyCarbon(molecule3D, sgi)) {
                    weight = DescriptorWeightsHelper.LABEL_WEIGHT_MANDATORY;
                }
                if(SubGraphIndices.isCharged(molecule3D, sgi)) {
                    weight = DescriptorWeightsHelper.LABEL_WEIGHT_MANDATORY;
                }
            } else {
                weight = DescriptorWeightsHelper.LABEL_WEIGHT_MANDATORY;
            }
            for (int atomIndex : sgi.getAtomIndices()) {
                weights[atomIndex] = weight;
            }
        }

        return weights;
    }

    public static int [] getBasisWeightCategories(int atoms){
        int [] weights = new int[atoms];
        Arrays.fill(weights, DescriptorWeightsHelper.LABEL_WEIGHT_NORMAL);
        return weights;
    }
    public static void setWeightCategories(int [] basisArray, int category, int [] atsInd){
        for (int atInd : atsInd) {
            basisArray[atInd]=category;
        }
    }

    public static String toStringCategoryWeights(int [] weightLabels){
        StringBuilder weightBuilder = new StringBuilder();
        for (int weightLabel : weightLabels) {
            weightBuilder.append((char)('0'+weightLabel));
        }
        return weightBuilder.toString();
    }

    /**
     *
     * @param s string with single digits '0123456789'.
     * @return
     */
    public static int [] parseSingleDigitString(String s) {
        if (s == null)
            return null;
        int[] arr = new int[s.length()];
        for (int i = 0; i < s.length(); i++) {
            int c = Integer.parseInt(Character.toString(s.charAt(i)));
            arr[i] = c;
        }
        return arr;
    }
}
