package com.actelion.research.chem.descriptor.flexophore.generator;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.histogram.MatrixBasedHistogram;
import com.actelion.research.chem.*;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.descriptor.DescriptorHandlerFlexophore;
import com.actelion.research.chem.descriptor.flexophore.*;
import com.actelion.research.chem.descriptor.flexophore.Molecule3DFunctions;
import com.actelion.research.chem.descriptor.flexophoreV5.*;
import com.actelion.research.chem.descriptor.flexophore.redgraph.SubGraphExtractor;
import com.actelion.research.chem.descriptor.flexophore.redgraph.SubGraphIndices;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import com.actelion.research.util.Formatter;
import com.actelion.research.util.datamodel.DoubleArray;
import org.openmolecules.chem.conf.gen.ConformerGenerator;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;

/**
 * CreatorMolDistHistViz
 * <p>Copyright: Actelion Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * Created by korffmo1 on 19.02.16.
 */
public class CreatorMolDistHistViz {

    private static final boolean DEBUG = DescriptorHandlerFlexophore.DEBUG;

    private static final long SEED = 123456789;

    // Maximum number of tries to generate conformers with the torsion rule based conformer generator from Thomas Sander
    private static final int MAX_NUM_TRIES = 10000;

    private static final int CONF_GEN_TS = 0;

    public static final int CONF_GIVEN_SINGLE_CONFORMATION = 1;

    private static CreatorMolDistHistViz INSTANCE;


    private SubGraphExtractor subGraphExtractor;

    private ConformerGenerator conformerGenerator;

    private int conformationMode;

    private MoleculeStandardizer moleculeStandardizer;

    private long seed;

    // for debugging
    private boolean onlyOneConformer;


    public CreatorMolDistHistViz() {

        seed = SEED;

        subGraphExtractor = new SubGraphExtractor();

        conformerGenerator = new ConformerGenerator(seed, false);

        moleculeStandardizer = new MoleculeStandardizer();

        conformationMode = CONF_GEN_TS;

        // System.out.println("CreatorCompleteGraph conformationMode " + conformationMode);

    }

    public void setConformationMode(int conformationMode) {
        this.conformationMode = conformationMode;
    }

    public MolDistHistViz create(StereoMolecule molOrig) throws Exception {

        MolDistHistViz mdhv = null;

        switch (conformationMode) {
            case CONF_GEN_TS:
                mdhv = createMultipleConformations(molOrig);
                break;
            case CONF_GIVEN_SINGLE_CONFORMATION:
                mdhv = createFromGivenConformation(molOrig);
                break;
            default:
                throw new RuntimeException("Invalid conformation mode");

        }

        return mdhv;
    }


    /**
     * Conformation generator of Thomas Sander
     * @param molOrig
     * @return
     * @throws Exception
     */
    public MolDistHistViz createMultipleConformations(StereoMolecule molOrig) throws Exception {

        int nConformations = DescriptorHandlerFlexophore.NUM_CONFORMATIONS;

        StereoMolecule molStand = moleculeStandardizer.getStandardized(molOrig);

        molStand.ensureHelperArrays(Molecule.cHelperRings);

        Molecule3D molInPlace = new Molecule3D(molStand);

        molInPlace.ensureHelperArrays(Molecule.cHelperRings);

        InteractionAtomTypeCalculator.setInteractionTypes(molInPlace);

        List<SubGraphIndices> liFragment = subGraphExtractor.extract(molInPlace);

        if(DEBUG) {
            injectNewSeed();
        }

        conformerGenerator.initializeConformers(molInPlace, ConformerGenerator.STRATEGY_LIKELY_RANDOM, MAX_NUM_TRIES, false);

        int nAtoms = molInPlace.getAtoms();

        List<MultCoordFragIndex> liMultCoordFragIndex = new ArrayList<>();
        for (SubGraphIndices subGraphIndices : liFragment) {
            liMultCoordFragIndex.add(new MultCoordFragIndex(subGraphIndices.getAtomIndices()));
        }

        int ccConformationsGenerated = 0;

        Molecule3D molViz = null;
        for (int i = 0; i < nConformations; i++) {

            boolean conformerGenerated = generateConformerAndSetCoordinates(conformerGenerator, nAtoms, molInPlace);

            if(!conformerGenerated){
                break;
            }

            ccConformationsGenerated++;

            calcFragmentCenter(molInPlace, liMultCoordFragIndex);

            if(i==0){
                molViz = createMoleculeForVisualization(molInPlace, liMultCoordFragIndex);
            }
        }

        if(ccConformationsGenerated==0){
            throw new ExceptionConformationGenerationFailed("Impossible to generate one conformer!");
        }

        int nPotentialConformers = conformerGenerator.getPotentialConformerCount();

        onlyOneConformer = false;

        if((nPotentialConformers > 1) && (ccConformationsGenerated==1)){

            if(DEBUG) {
                System.out.println("CreatorCompleteGraph: only one conformer generated.");

                System.out.println("Seed " + seed);

                System.out.println("Potential conformer count " + nPotentialConformers);

                Canonizer can = new Canonizer(molInPlace);

                System.out.println(can.getIDCode());
            }

            onlyOneConformer = true;

        }

        MolDistHistViz mdhv = create(liMultCoordFragIndex, molViz);

        return mdhv;
    }

    private static MolDistHistViz create(List<MultCoordFragIndex> liMultCoordFragIndex, Molecule3D molecule3D){

        MolDistHistViz molDistHistViz = new MolDistHistViz(liMultCoordFragIndex.size(), molecule3D);

        List<PPNodeViz> liPPNodeViz = new ArrayList<>();
        for (int i = 0; i < liMultCoordFragIndex.size(); i++) {

            int [] arrIndexAtomFrag = liMultCoordFragIndex.get(i).getArrIndexFrag();

            PPNodeViz ppNodeViz = new PPNodeViz();
            ppNodeViz.setIndex(i);
            ppNodeViz.setCoordinates(liMultCoordFragIndex.get(i).liCoord.get(0));

            for (int index : arrIndexAtomFrag) {
                int interactionType = molecule3D.getInteractionAtomType(index);
                ppNodeViz.add(interactionType);
                ppNodeViz.addIndexOriginalAtom(index);
            }
            liPPNodeViz.add(ppNodeViz);
        }

        molDistHistViz.set(liPPNodeViz);

        for (int i = 0; i < liMultCoordFragIndex.size(); i++) {
            for (int j = i+1; j < liMultCoordFragIndex.size(); j++) {
                byte [] arrDistHist = MultCoordFragIndex.getDistHist(liMultCoordFragIndex.get(i), liMultCoordFragIndex.get(j));
                // System.out.println(StringFunctions.toString(arrDistHist));
                molDistHistViz.setDistHist(i,j,arrDistHist);
            }
        }

        molDistHistViz.realize();

        return molDistHistViz;

    }

    private Molecule3D createMoleculeForVisualization(Molecule3D molecule3D, List<MultCoordFragIndex> liMultCoordFragIndex) {

        Molecule3D molCenter = new Molecule3D(molecule3D);

        molCenter.ensureHelperArrays(Molecule.cHelperRings);

        for (MultCoordFragIndex multCoordFragIndex : liMultCoordFragIndex) {

            int [] arrAtomIndexList = multCoordFragIndex.getArrIndexFrag();

            // Calculate center coordinates.
            Coordinates coordCenter = Molecule3DFunctions.getCenterGravity(molCenter, arrAtomIndexList);

            for (int at = 0; at < arrAtomIndexList.length; at++) {

                // Atom type.
                int interactionType = molCenter.getInteractionAtomType(arrAtomIndexList[at]);

                if(interactionType == -1){
                    continue;
                }

                int iAtomicNo = molCenter.getAtomicNo(arrAtomIndexList[at]);

                molCenter.setAtomFlag(arrAtomIndexList[at], Molecule3D.FLAG1, true);

                int indexOriginalAtom = arrAtomIndexList[at];

                int indexAtm = molCenter.addAtom(iAtomicNo);

                molCenter.setInteractionAtomType(indexAtm, interactionType);

                String sOrigIndex = Integer.toString(indexOriginalAtom);

                molCenter.setAtomChainId(indexAtm, sOrigIndex);

                // Set the center coordinates
                molCenter.setCoordinates(indexAtm, coordCenter);

                molCenter.setAtomFlag(indexAtm, Molecule3DFunctions.FLAG_CENTER_ATOM, true);

                molCenter.setPPP(indexAtm, arrAtomIndexList);

            }
        }

        molCenter.ensureHelperArrays(Molecule.cHelperRings);

        return molCenter;
    }

    public void injectNewSeed(){

        seed = new Date().getTime();

        conformerGenerator = new ConformerGenerator(seed, false);
    }


    public boolean isOnlyOneConformer() {
        return onlyOneConformer;
    }


    public MolDistHistViz createFromGivenConformation(StereoMolecule molOrig) {

        // Just to make sure we will not change something in the original molecule.
        Molecule3D molStart = new Molecule3D(molOrig);

        molStart.ensureHelperArrays(Molecule.cHelperRings);

        InteractionAtomTypeCalculator.setInteractionTypes(molStart);

        StereoMolecule molInPlace = new Molecule3D(molStart);

        molInPlace.ensureHelperArrays(Molecule.cHelperRings);

        List<SubGraphIndices> liFragment = subGraphExtractor.extract(molInPlace);

        List<MultCoordFragIndex> liMultCoordFragIndex = new ArrayList<>();
        for (SubGraphIndices subGraphIndices : liFragment) {
            liMultCoordFragIndex.add(new MultCoordFragIndex(subGraphIndices.getAtomIndices()));
        }

        calcFragmentCenter(molStart, liMultCoordFragIndex);

        MolDistHistViz mdhv = create(liMultCoordFragIndex, molStart);

        return mdhv;
    }

    /**
     * 08.03.2017 Method set to public for debugging purposes.
     * @param conformerGenerator
     * @param nAtoms
     * @param molInPlace
     * @return
     */
    public static boolean generateConformerAndSetCoordinates(ConformerGenerator conformerGenerator, int nAtoms, Molecule3D molInPlace){

        boolean nextConformerAvailable = false;

        Conformer conformer = conformerGenerator.getNextConformer();

        if(conformer != null){

            // System.out.println(ConformerUtils.toStringDistances(ConformerUtils.getDistanceMatrix(conformer)));

            for (int i = 0; i < nAtoms; i++) {
                double x = conformer.getX(i);
                double y = conformer.getY(i);
                double z = conformer.getZ(i);

                molInPlace.setAtomX(i,x);
                molInPlace.setAtomY(i,y);
                molInPlace.setAtomZ(i,z);
            }
            nextConformerAvailable = true;
        }

        return nextConformerAvailable;
    }

    /**
     * Calculates the center of the fragments and stores the coordinates in MultCoordFragIndex.
     * @param molecule3D
     * @param liMultCoordFragIndex
     */
    private static void calcFragmentCenter(Molecule3D molecule3D, List<MultCoordFragIndex> liMultCoordFragIndex) {

        for (MultCoordFragIndex multCoordFragIndex : liMultCoordFragIndex) {

            int [] arrAtomIndexList = multCoordFragIndex.getArrIndexFrag();

            // Calculate center coordinates.
            Coordinates coordCenter = Molecule3DFunctions.getCenterGravity(molecule3D, arrAtomIndexList);

            multCoordFragIndex.addCoord(coordCenter);

        }
    }

    public static int getConformations2Generate(StereoMolecule mol, int maxNumConf) {

        int conf = 0;

        int rot = mol.getRotatableBondCount();

        double dConf = Math.pow(3,rot);

        if(dConf > maxNumConf){

            conf = maxNumConf;

        } else {

            conf = (int)dConf;

        }

        return conf;
    }

    private static class InteractionType {

        int type;

        String description;

        String comment;

        public InteractionType(int type, String description, String comment) {
            this.type = type;
            this.description = description;
            this.comment = comment;
        }
    }

    public static CreatorMolDistHistViz getInstance(){

        if(INSTANCE==null){

            INSTANCE = new CreatorMolDistHistViz();

        }

        return INSTANCE;
    }

    private static class MultCoordFragIndex {

        private int [] arrIndexFrag;

        private List<Coordinates> liCoord;

        public MultCoordFragIndex(int[] arrIndexFrag) {
            this.arrIndexFrag = arrIndexFrag;
            liCoord = new ArrayList<>();
        }

        public void addCoord(Coordinates coordinates){
            liCoord.add(coordinates);
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

            double [] arrHist = maHist.getRow(2);

//            System.out.println(StringFunctions.toString(daDistances.get()));
//            System.out.println(StringFunctions.toString(arrHist));

            int countValuesInHistogram = 0;

            for (int i = 0; i < arrHist.length; i++) {
                countValuesInHistogram += arrHist[i];
            }

            // Here, the percentage values for the histograms are calculated.
            byte [] arrHistPercent = new byte [maHist.getColDim()];

            for (int i = 0; i < arrHist.length; i++) {
                arrHistPercent[i]= (byte)  (((arrHist[i] / countValuesInHistogram) * 100.0) + 0.5);
            }

//            System.out.println(StringFunctions.toString(arrHistPercent));

            return arrHistPercent;

        }

    }

}
