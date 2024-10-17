package com.actelion.research.chem.descriptor;

import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.descriptor.flexophore.ConstantsFlexophore;
import com.actelion.research.chem.descriptor.flexophore.generator.CreatorMolDistHistViz;
import com.actelion.research.chem.descriptor.flexophore.redgraph.SubGraphIndices;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.PheSAMolecule;
import com.actelion.research.chem.phesa.pharmacophore.pp.IPharmacophorePoint;
import com.actelion.research.chem.phesa.pharmacophore.pp.PPGaussian;
import com.actelion.research.util.datamodel.IntArray;

import java.util.ArrayList;
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
    public static final int LABEL_WEIGHT_HIGH_USER = 3;

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
    public int [] calcWeightLabelsFlexophore(Molecule3D molecule3D){
        List<SubGraphIndices> liSubGraphIndices = creatorMolDistHistViz.getSubGraphIndices(molecule3D);
        int [] weights = calcWeightLabels(liSubGraphIndices, molecule3D);
        return weights;
    }

    /**
     * - End standing pp points are set mandatory.
     * - Charged pp points are set mandatory.
     * @param liSubGraphIndices
     * @param molecule3D
     * @return array with dimension molecule3D.getAtoms().
     */
    public static int [] calcWeightLabels(List<SubGraphIndices> liSubGraphIndices, Molecule3D molecule3D){

        int [] weights = getBasisWeightLabels(molecule3D);

        for (int i = 0; i < liSubGraphIndices.size(); i++) {
            SubGraphIndices sgi = liSubGraphIndices.get(i);
            int indexWeight = DescriptorWeightsHelper.LABEL_WEIGHT_NORMAL;
            if (SubGraphIndices.isLinker(molecule3D, liSubGraphIndices, i)) {
                if(!SubGraphIndices.isOnlyCarbon(molecule3D, sgi)) {
                    indexWeight = DescriptorWeightsHelper.LABEL_WEIGHT_MANDATORY;
                }
                if(SubGraphIndices.isCharged(molecule3D, sgi)) {
                    indexWeight = DescriptorWeightsHelper.LABEL_WEIGHT_MANDATORY;
                }
            } else {
                indexWeight = DescriptorWeightsHelper.LABEL_WEIGHT_MANDATORY;
            }
            for (int atomIndex : sgi.getAtomIndices()) {
                weights[atomIndex] = indexWeight;
            }
        }

        return weights;
    }

    public static int [] getBasisWeightLabels(Molecule3D molecule3D){
        int [] weights = new int[molecule3D.getAtoms()];
        Arrays.fill(weights, DescriptorWeightsHelper.LABEL_WEIGHT_NORMAL);
        return weights;
    }

    public static int [] mergeWeightLabels(int[] arrWeightLabel, int[] arrWeightLabelUser) {

        if(arrWeightLabel.length!= arrWeightLabelUser.length){
            throw new RuntimeException("Weight labels differ in size!");
        }
        int[] arrWeightLabelMerged = new int[arrWeightLabel.length];
        for (int i = 0; i < arrWeightLabel.length; i++) {
            int label = arrWeightLabel[i];

            if(arrWeightLabelUser[i]== LABEL_WEIGHT_LOW){
                label = arrWeightLabelUser[i];
            } else if(arrWeightLabel[i] == LABEL_WEIGHT_MANDATORY && arrWeightLabelUser[i] == LABEL_WEIGHT_MANDATORY){
                label = LABEL_WEIGHT_HIGH_USER;
            } else if(arrWeightLabelUser[i] == LABEL_WEIGHT_MANDATORY){
                label = LABEL_WEIGHT_HIGH_USER;
            }
            arrWeightLabelMerged[i]=label;
        }

        return arrWeightLabelMerged;
    }
    public static String toStringWeightLabels(int [] weightLabels){
        StringBuilder weightBuilder = new StringBuilder();
        for (int weightLabel : weightLabels) {
            weightBuilder.append((char)('0'+weightLabel));
        }
        return weightBuilder.toString();
    }

    /**
     *
     * @param s string with single digits '123456789'.
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
