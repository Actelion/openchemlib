package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.descriptor.flexophore.generator.ConstantsFlexophoreGenerator;
import com.actelion.research.chem.descriptor.flexophore.generator.SubFlexophoreGenerator;
import com.actelion.research.util.ArrayUtils;
import com.actelion.research.util.ByteArray;

import java.util.*;

public class MolDistHistHelper {


    public static int getNumNodesNotEmpty(MolDistHist mdh){

        if(mdh==null)
            return 0;

        if(MolDistHistHelper.isEmptyMolDistHist(mdh))
            return 0;

        return mdh.getNumPPNodes();
    }

    /**
     * Creates a single {@link MolDistHist} object from several MolDistHist. Missing descriptor histograms are added
     * and filled with ones.
     * @param arr
     * @return
     */
    public static MolDistHist assemble (MolDistHist ... arr){

        int nNodesSum = 0;
        int maxNumNodes = 0;
        for (MolDistHist mdhFrag : arr) {
            int nNodes = mdhFrag.getNumPPNodes();
            nNodesSum += nNodes;
            if(nNodes>maxNumNodes)
                maxNumNodes=nNodes;
        }
        MolDistHist mdh = new MolDistHist(nNodesSum);
        int [] [] arrMapIndexNew = new int[arr.length][maxNumNodes];
        int indexNew = 0;
        for (int i = 0; i < arr.length; i++) {
            int n = arr[i].getNumPPNodes();
            for (int j = 0; j < n; j++) {
                arrMapIndexNew[i][j]=indexNew;
                indexNew++;
                mdh.addNode(arr[i].getNode(j));
            }
        }

        for (int i = 0; i < arr.length; i++) {
            int n = arr[i].getNumPPNodes();
            MolDistHist mdhFrag = arr[i];

            for (int j = 0; j < n; j++) {
                for (int k = j+1; k < n; k++) {
                    byte [] arrDistHist = mdhFrag.getDistHist(j,k);
                    int index1 = arrMapIndexNew[i][j];
                    int index2 = arrMapIndexNew[i][k];
                    mdh.setDistHist(index1, index2, arrDistHist);
                }
            }
        }

        // Fill the missing distance histograms
        for (int i = 0; i < nNodesSum; i++) {
            for (int j = i+1; j < nNodesSum; j++) {
                byte [] arrDistHist = mdh.getDistHist(i,j);
                if(isZero(arrDistHist)){
                    Arrays.fill(arrDistHist, (byte)1);
                    mdh.setDistHist(i, j, arrDistHist);
                }
                // System.out.println(ArrayUtils.toString(arrDistHist));
            }
        }

        return mdh;
    }

    /**
     *
     * @param arr
     * @return canonical representation of nodes without distance histogram
     */
    public static MolDistHist assembleNoDistHist (MolDistHist ... arr){

        List<PPNode> liPPNode = new ArrayList<>();
        for (MolDistHist mdhFrag : arr) {
            int nNodes = mdhFrag.getNumPPNodes();
            for (int i = 0; i < nNodes; i++) {
                PPNode node = mdhFrag.getNode(i);
                if(node.getInteractionTypeCount()>0){
                    liPPNode.add(node);
                }
            }
        }

        if(liPPNode.size()==0){
            return null;
        }

        Collections.sort(liPPNode);
        MolDistHist mdh = new MolDistHist(liPPNode.size());
        for (PPNode ppNode : liPPNode) {
            mdh.addNode(ppNode);
        }
        mdh.realize();
        return mdh;
    }

    public static boolean isZero(byte [] b){
        boolean z=true;
        for (byte v : b) {
            if(v!=0){
                z=false;
                break;
            }
        }
        return z;
    }

    public static String toStringDistHist(MolDistHist mdh){
        StringBuilder sb = new StringBuilder();
        int nNodes = mdh.getNumPPNodes();
        for (int i = 0; i < nNodes; i++) {
            for (int j = i+1; j < nNodes; j++) {
                byte [] arrDistHist = mdh.getDistHist(i,j);

                sb.append(ArrayUtils.toString(arrDistHist));
                sb.append("\n");
            }
        }
        return sb.toString();
    }

    public static void setDistHistToOne(MolDistHist mdh){
        int nNodes = mdh.getNumPPNodes();
        for (int i = 0; i < nNodes; i++) {
            for (int j = i+1; j < nNodes; j++) {
                byte [] arrDistHist = mdh.getDistHist(i,j);
                Arrays.fill(arrDistHist, (byte) 1);
                mdh.setDistHist(i,j, arrDistHist);
            }
        }
    }

    public static MolDistHist getEmptyMolDistHist(){
        PPNode ppNode0 = new PPNode();
        ppNode0.realize();

        MolDistHist mdhEmpty = new MolDistHist(1);
        mdhEmpty.addNode(ppNode0);
        mdhEmpty.realize();

        return mdhEmpty;
    }
    public static boolean isEmptyMolDistHist(MolDistHist mdh){

        boolean empty = true;
        if(mdh.getNumPPNodes()>1){
            empty = false;
        } else if(mdh.getNumPPNodes()==1){
            if(mdh.getNode(0).getInteractionTypeCount()>0){
                empty = false;
            }
        }

        return empty;
    }

    public static MolDistHist getMostDistantPairOfNodes (MolDistHist mdh){

        int n = mdh.getNumPPNodes();

        int maxMedBinIndex = 0;
        int indexNode1=-1;
        int indexNode2=-1;

        for (int i = 0; i < n; i++) {
            for (int j = i+1; j < n; j++) {
                byte [] arr = mdh.getDistHist(i, j);
                int medBinInd = DistHistHelper.getMedianBin(arr);
                if(medBinInd>maxMedBinIndex){
                    maxMedBinIndex=medBinInd;
                    indexNode1 = i;
                    indexNode2 = j;
                }
            }
        }
        int [] arrIndexNodes = new int[2];
        arrIndexNodes[0]=indexNode1;
        arrIndexNodes[1]=indexNode2;

        MolDistHist mdhSub = SubFlexophoreGenerator.getSubFragment(mdh, arrIndexNodes);

        return mdhSub;
    }
    public static MolDistHist getMostDistantPairOfNodesOneHeteroAtom (MolDistHist mdh){

        int n = mdh.getNumPPNodes();

        int maxMedBinIndex = 0;
        int indexNode1=-1;
        int indexNode2=-1;

        boolean [] arrHetero = new boolean[n];

        for (int i = 0; i < n; i++) {
            arrHetero[i]=mdh.getNode(i).containsHetero();
        }

        for (int i = 0; i < n; i++) {
            for (int j = i+1; j < n; j++) {
                if (arrHetero[i] || arrHetero[j]) {
                    byte[] arr = mdh.getDistHist(i, j);
                    int medBinInd = DistHistHelper.getMedianBin(arr);
                    if (medBinInd > maxMedBinIndex) {
                        maxMedBinIndex = medBinInd;
                        indexNode1 = i;
                        indexNode2 = j;
                    }
                }
            }
        }
        int [] arrIndexNodes = new int[2];
        arrIndexNodes[0]=indexNode1;
        arrIndexNodes[1]=indexNode2;

        MolDistHist mdhSub = SubFlexophoreGenerator.getSubFragment(mdh, arrIndexNodes);

        return mdhSub;
    }

    public static boolean areNodesEqual(MolDistHist mdh1, MolDistHist mdh2){
        boolean eq = true;
        if(mdh1.getNumPPNodes() != mdh2.getNumPPNodes()){
            eq = false;
        } else {
            for (int i = 0; i < mdh1.getNumPPNodes(); i++) {
                if(!mdh1.getNode(i).equals(mdh2.getNode(i))){
                    eq = false;
                    break;
                }
            }
        }
        return eq;
    }


    /**
     *
     * @param col
     * @return
     */
    public static MolDistHist createFromNodes (Collection<PPNode> col){
        MolDistHistViz molDistHistViz = new MolDistHistViz(col.size());
        for (PPNode ppNode : col) {
            PPNodeViz ppNodeViz = new PPNodeViz(ppNode);
            molDistHistViz.addNode(ppNodeViz);
        }
        molDistHistViz.realize();
        return molDistHistViz.getMolDistHist();
    }

    public static MolDistHist [] toArray(List<MolDistHist> li){
        MolDistHist [] a = new MolDistHist[li.size()];
        for (int i = 0; i < li.size(); i++) {
            a[i]=li.get(i);
        }
        return a;
    }

    public static void reNormalizeDistHist(MolDistHist mdh, int margin){
        int n = mdh.getNumPPNodes();
        for (int i = 0; i < n; i++) {
            for (int j = i+1; j < n; j++) {
                byte [] b = mdh.getDistHist(i,j);
                int c = DistHistHelper.count(b);
                if(Math.abs(ConstantsFlexophoreGenerator.SUM_VAL_HIST-c) > margin){
                    byte [] distHistNew = DistHistHelper.normalize(b);
                    mdh.setDistHist(i,j,distHistNew);
//                    int c2 = DistHistHelper.count(distHistNew);
//                    System.out.println(c2);
                }
            }
        }
    }

}
