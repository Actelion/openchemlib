package com.actelion.research.chem.descriptor.flexophore.generator;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.histogram.MatrixBasedHistogram;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.descriptor.flexophore.DistHistHelper;
import com.actelion.research.chem.descriptor.flexophore.redgraph.SubGraphIndices;
import com.actelion.research.util.Formatter;
import com.actelion.research.util.datamodel.DoubleArray;

import java.util.ArrayList;
import java.util.List;

public class MultCoordFragIndex {

    // Indices of the atoms in the molecule.
    private int [] arrIndexFrag;

    // List with center of gravity for multiple conformations.
    // One conformation describes one center of gravity for the fragments indicated by the index array.
    private List<Coordinates> liCoord;

    public MultCoordFragIndex(int[] arrIndexFrag) {
        this.arrIndexFrag = arrIndexFrag;
        liCoord = new ArrayList<>();
    }


    public MultCoordFragIndex getDeepClone(){
        int [] a = new int[arrIndexFrag.length];
        System.arraycopy(arrIndexFrag, 0, a, 0, a.length);
        MultCoordFragIndex m = new MultCoordFragIndex(a);
        for (Coordinates coordinates : liCoord) {
            m.liCoord.add(new Coordinates(coordinates));
        }
        return m;
    }

    public boolean isEqualIndices(SubGraphIndices sgi){
        return sgi.equalIndices(arrIndexFrag);
    }

    public void addCoord(Coordinates coordinates){
        liCoord.add(coordinates);
    }

    public List<Coordinates> getCoordinates() {
        return liCoord;
    }

    public int[] getArrIndexFrag() {
        return arrIndexFrag;
    }

    public static byte [] getDistHist(MultCoordFragIndex m1, MultCoordFragIndex m2) {

        final double minRangeDistanceHistogram =0;

        final double maxRangeDistanceHistogram = ConstantsFlexophoreGenerator.RANGE_HISTOGRAM;

        final int bins = ConstantsFlexophoreGenerator.BINS_HISTOGRAM;

        int n1 = m1.liCoord.size();
        int n2 = m2.liCoord.size();

        if(n1!=n2){
            throw new RuntimeException("Number of distances differ for two fragments sets!");
        }

        DoubleArray daDistances = new DoubleArray(n1);

        for (int i = 0; i < n1; i++) {
            double dist = m1.liCoord.get(i).distance(m2.liCoord.get(i));
            if(dist > ConstantsFlexophoreGenerator.RANGE_HISTOGRAM) {
                throw new RuntimeException("Distance " + Formatter.format1(dist) + " between two pharmacophore points exceeded maximum histogram range.");
            }
            daDistances.add(dist);
        }

        Matrix maBins = MatrixBasedHistogram.getHistogramBins(minRangeDistanceHistogram,maxRangeDistanceHistogram, bins);

        Matrix maHist =  MatrixBasedHistogram.getHistogram(daDistances.get(), maBins);

        double [] arrHistRaw = maHist.getRow(2);

//            System.out.println(StringFunctions.toString(daDistances.get()));
//            System.out.println(StringFunctions.toString(arrHist));

        return DistHistHelper.normalize(arrHistRaw);

    }

}
