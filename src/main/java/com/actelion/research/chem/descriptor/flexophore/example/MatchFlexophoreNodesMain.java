package com.actelion.research.chem.descriptor.flexophore.example;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerFlexophore;
import com.actelion.research.chem.descriptor.flexophore.ModelSolutionSimilarity;
import com.actelion.research.chem.descriptor.flexophore.MolDistHist;
import com.actelion.research.chem.descriptor.flexophore.MolDistHistViz;
import com.actelion.research.chem.descriptor.flexophore.PPNodeViz;
import com.actelion.research.chem.descriptor.flexophore.completegraphmatcher.ObjectiveBlurFlexophoreHardMatchUncovered;
import com.actelion.research.gui.table.ChemistryRenderPanel;
import com.actelion.research.util.ArrayUtils;
import com.actelion.research.util.Formatter;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;


public class MatchFlexophoreNodesMain {

    // DUD_E hs90a crystal
    public static final String IDCODE_100 = "ffchb@DXXJWb`dLbbRbrrTTRRMNRUGL@QSPPQQ@@@";
    // DUD_E hs90a seed CHEMBL192145
    public static final String IDCODE_101 = "fi{iA@IJcdn`Xa@cHhhdmDeDdmDjfBM[MXHHZj`bJp@@@";


    public static void main(String[] args) {

        String idcode = IDCODE_100;
        String idcodeQuery = IDCODE_101;

        pairMatch(idcode, idcodeQuery);
    }

    public static void pairMatch(String idcodeBase, String idcodeQuery) {
        DescriptorHandlerFlexophore dhFlexophore = new DescriptorHandlerFlexophore();

        ObjectiveBlurFlexophoreHardMatchUncovered objectiveBlurFlexophoreHardMatchUncovered = dhFlexophore.getObjectiveCompleteGraph();

        IDCodeParser parser = new IDCodeParser();

        MolDistHistViz mdhvBase = create(parser, dhFlexophore, idcodeBase);

        // Fetch original tom indices
        System.out.println("Num pp nodes base " + mdhvBase.getNumPPNodes());
        List<int[]> liBaseArrayIndexAtom = new ArrayList<>();
        for (int i = 0; i < mdhvBase.getNumPPNodes(); i++) {
            int [] arrIndexAt = mdhvBase.getNode(i).getArrayIndexOriginalAtoms();
            liBaseArrayIndexAtom.add(arrIndexAt);
        }

        // Fetch original tom indices
        MolDistHistViz mdhvQuery = create(parser, dhFlexophore, idcodeQuery);
        System.out.println("Num pp nodes query " + mdhvQuery.getNumPPNodes());
        List<int[]> liQueryArrayIndexAtom = new ArrayList<>();
        for (int i = 0; i < mdhvQuery.getNumPPNodes(); i++) {
            int [] arrIndexAt = mdhvQuery.getNode(i).getArrayIndexOriginalAtoms();
            liQueryArrayIndexAtom.add(arrIndexAt);
        }

        System.out.println(mdhvBase.toString());
        System.out.println(mdhvQuery.toString());

        //
        // Just to show that it works with MolDistHist as with MolDistHistViz
        //
        MolDistHist mdhBase = mdhvBase.getMolDistHist();
        MolDistHist mdhQuery = mdhvQuery.getMolDistHist();

        ModelSolutionSimilarity modelSolutionSimilarity = dhFlexophore.getBestMatch(mdhBase, mdhQuery);
        int heap = modelSolutionSimilarity.getSizeHeap();

        // System.out.println(Formatter.format3(sim) + "\t" + Formatter.format3(simDH));
        // System.out.println(Formatter.format3(sim) +"\t"+ Formatter.format3(simNormalized));

        objectiveBlurFlexophoreHardMatchUncovered.setBase(mdhBase);
        objectiveBlurFlexophoreHardMatchUncovered.setQuery(mdhQuery);

        for (int i = 0; i <heap; i++) {

            int indexQuery = modelSolutionSimilarity.getIndexQueryFromHeap(i);
            int indexBase = modelSolutionSimilarity.getIndexBaseFromHeap(i);

            PPNodeViz ppvBase = mdhvBase.getNode(indexBase);

            int [] arrAtomIndexBase = liBaseArrayIndexAtom.get(indexBase);

            PPNodeViz ppvQuery = mdhvQuery.getNode(indexQuery);

            int [] arrAtomIndexQuery = liQueryArrayIndexAtom.get(indexQuery);

            System.out.println(ppvBase.toString());
            System.out.println(ArrayUtils.toString(arrAtomIndexBase));

            System.out.println(ppvQuery.toString());
            System.out.println(ArrayUtils.toString(arrAtomIndexQuery));
            // System.out.println(Formatter.format3((double)ppv.getSimilarityMappingNodes()) + "\t" + Formatter.format3((double)ppvQuery.getSimilarityMappingNodes()));

            float simHistogram = objectiveBlurFlexophoreHardMatchUncovered.getSimilarityHistogramsForNode(modelSolutionSimilarity, i);

            System.out.println("Node similarity " + Formatter.format3((double)modelSolutionSimilarity.getSimilarityNode(indexQuery)) + ", histogram similarity " + Formatter.format3((double)simHistogram) + ".");

            System.out.println();
        }

        System.out.println("Un-normalized similarity " + Formatter.format3(modelSolutionSimilarity.getSimilarity()));

    }


    public static MolDistHistViz create(IDCodeParser parser, DescriptorHandlerFlexophore dh, String idcode){
        StereoMolecule mol = parser.getCompactMolecule(idcode);
        mol.ensureHelperArrays(Molecule.cHelperRings);
        return dh.createVisualDescriptor(mol);
    }

}
