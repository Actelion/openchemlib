package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.descriptor.DescriptorWeightsHelper;
import com.actelion.research.chem.descriptor.flexophore.generator.MultCoordFragIndex;
import com.actelion.research.util.ArrayUtils;

import java.util.ArrayList;
import java.util.List;

public class MolDistHistVizHelper {

    public static void centerNodes(MolDistHistViz mdhv){
        Coordinates barycenter = getBarycenter(mdhv);
        for (int i = 0; i < mdhv.getNumPPNodes(); i++) {
            PPNodeViz node = mdhv.getNode(i);
            node.getCoordinates().x -= barycenter.x;
            node.getCoordinates().y -= barycenter.y;
            node.getCoordinates().z -= barycenter.z;
        }

        Molecule3D m = mdhv.getMolecule();
        Coordinates [] arrCoordinates = m.getCoordinates();
        for (int i = 0; i < arrCoordinates.length; i++) {
            Coordinates c = arrCoordinates[i];
            c.x -= barycenter.x;
            c.y -= barycenter.y;
            c.z -= barycenter.z;
            m.setCoordinates(i, c);
        }
    }

    public static Coordinates getBarycenter(MolDistHistViz mdhv){
        Coordinates [] arrCoordinates = new Coordinates[mdhv.getNumPPNodes()];
        for (int i = 0; i < mdhv.getNumPPNodes(); i++) {
            PPNodeViz node = mdhv.getNode(i);
            arrCoordinates[i] = node.getCoordinates();
        }
        return Coordinates.createBarycenter(arrCoordinates);
    }

    public static MolDistHistViz create(List<PPNodeVizMultCoord> liPPNodeVizMultCoord){

        MolDistHistViz molDistHistViz = new MolDistHistViz(liPPNodeVizMultCoord.size(), null);
        List<PPNodeViz> liPPNodeViz = new ArrayList<>();

        for (PPNodeVizMultCoord ppNodeViz : liPPNodeVizMultCoord) {
            liPPNodeViz.add(ppNodeViz);
        }
        molDistHistViz.set(liPPNodeViz);
        for (int i = 0; i < liPPNodeVizMultCoord.size(); i++) {
            for (int j = i+1; j < liPPNodeVizMultCoord.size(); j++) {
                byte [] arrDistHist =
                        MultCoordFragIndex.getDistHist(
                                liPPNodeVizMultCoord.get(i).multCoordFragIndex,
                                liPPNodeVizMultCoord.get(j).multCoordFragIndex);
                // System.out.println(StringFunctions.toString(arrDistHist));
                molDistHistViz.setDistHist(i,j,arrDistHist);
            }
        }
        molDistHistViz.realize();
        return molDistHistViz;
    }

    /**
     * Sets the weights in MolDistHistViz. Rule based unification of weight labels for a pharmacophore node.
     * @param mdhv
     * @param arrWeightLabel Array dimension must equal number of atoms in the molecule in mdhv.
     */
    public static void setWeights(MolDistHistViz mdhv, int [] arrWeightLabel){

        // The molecule in the descriptor contains the pharmacophore points as additional single atoms.
        Molecule3D m3D = new Molecule3D(mdhv.getMolecule());
        m3D.ensureHelperArrays(Molecule.cHelperRings);
        m3D.stripSmallFragments();

        if(m3D.getAtoms()!=arrWeightLabel.length){
            throw new RuntimeException("Weight vector differs in dimension to number of atoms!");
        }

        //
        // Rule based unification of weight labels for a pharmacophore node
        //
        for (PPNodeViz ppNodeViz : mdhv.getNodes()) {
            int [] a = ppNodeViz.getArrayIndexOriginalAtoms();
            int [] w = new int[a.length];
            for (int i = 0; i < a.length; i++) {
                w[i]=arrWeightLabel[a[i]];
            }

            int maxWeightLabel = ArrayUtils.max(w);

            // If label is not high, low is king. Means the standard wight label is overruled.
            if(maxWeightLabel< DescriptorWeightsHelper.LABEL_WEIGHT_MANDATORY){
                maxWeightLabel = ArrayUtils.min(w);
            }

            if(DescriptorWeightsHelper.LABEL_WEIGHT_HIGH_USER == maxWeightLabel){
                mdhv.addMandatoryPharmacophorePoint(ppNodeViz.getIndex());
                mdhv.setNodeWeight(ppNodeViz.getIndex(), ConstantsFlexophore.VAL_WEIGHT_HIGH_USER);
            } else if(DescriptorWeightsHelper.LABEL_WEIGHT_MANDATORY ==maxWeightLabel){
                mdhv.addMandatoryPharmacophorePoint(ppNodeViz.getIndex());
                mdhv.setNodeWeight(ppNodeViz.getIndex(), ConstantsFlexophore.VAL_WEIGHT_HIGH);
            } else if(DescriptorWeightsHelper.LABEL_WEIGHT_LOW ==maxWeightLabel){
                mdhv.setNodeWeight(ppNodeViz.getIndex(), ConstantsFlexophore.VAL_WEIGHT_LOW);
            }
        }
    }
}
