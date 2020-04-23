package com.actelion.research.chem.descriptor.flexophore.generator;

import com.actelion.research.calc.filter.SlidingWindow;
import com.actelion.research.chem.descriptor.flexophore.MolDistHist;
import com.actelion.research.chem.descriptor.flexophore.PPNode;
import com.actelion.research.chem.descriptor.sphere.ScaleClasses;
import com.actelion.research.util.BurtleHasher;
import com.actelion.research.util.ScaleLabel;
import com.actelion.research.util.datamodel.ByteVec;
import com.actelion.research.util.datamodel.IntegerDouble;

import java.util.List;

public class VectophoreGenerator {


    public static final double THRESH_SIMILARITY = 0.75;

    private static final double [] FILTER = {0.25,0.5,0.25};

    private static final int HASH_BITS = 10;
    private static final int HASH_INIT = 13;
    public static final int DESCRIPTOR_SIZE = (1 << HASH_BITS);


    private InteractionTypeMap interactionTypeMap;

    private SlidingWindow slidingWindow;

    public VectophoreGenerator() {

        interactionTypeMap = new InteractionTypeMap(THRESH_SIMILARITY);

        slidingWindow = new SlidingWindow(FILTER);
    }

    public float [] create(MolDistHist mdh) {

        float [] arr = new float[DESCRIPTOR_SIZE];

        int nodes = mdh.getNumPPNodes();

        for (int i = 0; i < nodes; i++) {
            PPNode node1 = mdh.getNode(i);
            int ats1 = node1.getInteractionTypeCount();

            for (int j = i+1; j < nodes; j++) {

                PPNode node2 = mdh.getNode(j);
                int ats2 = node2.getInteractionTypeCount();

                byte [] arrDistHist = mdh.getDistHist(i,j);

                byte [] arrDistHistFilt = slidingWindow.filter(arrDistHist);

                for (int k = 0; k < ats1; k++) {

                    int atomType1 = node1.getInteractionType(k);

                    List<IntegerDouble> liAtomTypeSimilarity1 = interactionTypeMap.getSimilars(atomType1);

                    for (int l = 0; l < ats2; l++) {

                        int atomType2 = node2.getInteractionType(l);
                        List<IntegerDouble> liAtomTypeSimilarity2 = interactionTypeMap.getSimilars(atomType2);

                        for (IntegerDouble atomTypeSim1 : liAtomTypeSimilarity1) {

                            for (IntegerDouble atomTypeSim2 : liAtomTypeSimilarity2) {

                                int hash = 0;
                                byte [] a, b = null;
                                if(atomTypeSim1.getInt() > atomTypeSim2.getInt()){
                                    a = ByteVec.toByteArray(atomTypeSim1.getInt());
                                    b = ByteVec.toByteArray(atomTypeSim2.getInt());
                                } else {
                                    b = ByteVec.toByteArray(atomTypeSim1.getInt());
                                    a = ByteVec.toByteArray(atomTypeSim2.getInt());
                                }

                                byte [] c = new byte[6];
                                System.arraycopy(a, 1, c, 0, 3);
                                System.arraycopy(b, 1, c, 3, 3);

                                double simProduct = atomTypeSim1.getDouble() * atomTypeSim2.getDouble();

                                int startIndexInVector = BurtleHasher.hashlittle(c, HASH_INIT);
                                startIndexInVector = (startIndexInVector & BurtleHasher.hashmask(HASH_BITS));

                                for (int m = 0; m < arrDistHistFilt.length; m++) {

                                    double v = arrDistHistFilt[m]*simProduct;

                                    int index = startIndexInVector + m;

                                    if(index>=arr.length){
                                        index = index - arr.length;
                                    }
                                    arr[index]+=v;
                                }
                            }
                        }
                    }
                }
            }
        }

//        for (int i = 0; i < arr.length; i++) {
//            arr[i] /= 10000;
//        }

        return arr;
    }
}
