package com.actelion.research.chem.descriptor.flexophore.example;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerFlexophore;
import com.actelion.research.chem.descriptor.flexophore.ModelSolutionSimilarity;
import com.actelion.research.chem.descriptor.flexophore.MolDistHistViz;
import com.actelion.research.chem.descriptor.flexophore.PPNodeViz;
import com.actelion.research.util.Formatter;


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

    public static void pairMatch(String idcode, String idcodeQuery) {
        DescriptorHandlerFlexophore dhFlexophore = new DescriptorHandlerFlexophore();
        IDCodeParser parser = new IDCodeParser();

        MolDistHistViz mdh = create(parser, dhFlexophore, idcode);

        System.out.println("Num pp nodes base " + mdh.getNumPPNodes());

        MolDistHistViz mdhQuery = create(parser, dhFlexophore, idcodeQuery);
        System.out.println("Num pp nodes query " + mdhQuery.getNumPPNodes());

        System.out.println(mdh.toString());
        System.out.println(mdhQuery.toString());

        ModelSolutionSimilarity modelSolutionSimilarity = dhFlexophore.getBestMatch(mdh, mdhQuery);
        int heap = modelSolutionSimilarity.getSizeHeap();

        // System.out.println(Formatter.format3(sim) + "\t" + Formatter.format3(simDH));
        // System.out.println(Formatter.format3(sim) +"\t"+ Formatter.format3(simNormalized));

        for (int i = 0; i <heap; i++) {

            int indexQuery = modelSolutionSimilarity.getIndexQueryFromHeap(i);
            int indexBase = modelSolutionSimilarity.getIndexBaseFromHeap(i);

            PPNodeViz ppv = mdh.getNode(indexBase);
            PPNodeViz ppvQuery = mdhQuery.getNode(indexQuery);
            System.out.println(ppv.toString());
            System.out.println(ppvQuery.toString());
            // System.out.println(Formatter.format3((double)ppv.getSimilarityMappingNodes()) + "\t" + Formatter.format3((double)ppvQuery.getSimilarityMappingNodes()));
            System.out.println(Formatter.format3((double)modelSolutionSimilarity.getSimilarityNode(indexQuery)));

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
