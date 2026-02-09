package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.descriptor.flexophore.generator.ConstantsFlexophoreGenerator;
import com.actelion.research.chem.descriptor.flexophore.generator.SubFlexophoreGenerator;
import com.actelion.research.util.ArrayUtils;
import com.actelion.research.util.ByteArray;

import java.util.*;
/**
 * Copyright (c) 2026
 * Alipheron AG
 * Hochbergerstrasse 60C
 * CH-4057 Basel
 * Switzerland
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
 * 3. Neither the name of the copyright holder nor the
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
 * @author Modest v. Korff
 */
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

    public static Comparator<MolDistHist> getComparator() {

        return new Comparator<MolDistHist>() {
            @Override
            public int compare(MolDistHist o1, MolDistHist o2) {
                int cmp = 0;

                if(o1.getNumPPNodes()>o2.getNumPPNodes()){
                    cmp=1;
                } else if(o1.getNumPPNodes()<o2.getNumPPNodes()){
                    cmp=-1;
                } else {
                    for (int i = 0; i < o1.getNumPPNodes(); i++) {
                        PPNode n1 = o1.getNode(i);
                        PPNode n2 = o2.getNode(i);
                        cmp = n1.compareTo(n2);
                        if(cmp!=0){
                            break;
                        }
                    }
                }
                return cmp;
            }
        };
    }

}
