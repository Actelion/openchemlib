package com.actelion.research.chem.descriptor;

import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.descriptor.flexophore.generator.CreatorMolDistHistViz;
import com.actelion.research.chem.descriptor.flexophore.redgraph.SubGraphIndices;

import java.util.Arrays;
import java.util.List;

/*
 
 Copyright (c) 2024 Alipheron AG. All rights reserved.
 
 This file is part of the Alipheron AG software suite.
 
 Licensed under the Alipheron AG Software License Agreement (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at the company's official website or upon request.
 
 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License. 
 
 Created by Modest von Korff 
 26/07/2024
 
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
