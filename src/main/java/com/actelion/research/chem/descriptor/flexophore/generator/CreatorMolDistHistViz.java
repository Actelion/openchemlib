/*
 * Copyright (c) 1997 - 2016
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
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
 * 3. Neither the name of the the copyright holder nor the
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

package com.actelion.research.chem.descriptor.flexophore.generator;

import com.actelion.research.calc.ArrayUtilsCalc;
import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.*;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.descriptor.DescriptorHandlerFlexophore;
import com.actelion.research.chem.descriptor.flexophore.*;
import com.actelion.research.chem.descriptor.flexophore.redgraph.SubGraphExtractor;
import com.actelion.research.chem.descriptor.flexophore.redgraph.SubGraphIndices;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import com.actelion.research.util.StringFunctions;

import java.util.*;

/**
 * CreatorMolDistHistViz
 * Created by korffmo1 on 19.02.16.
 */
public class CreatorMolDistHistViz {

    private static final boolean DEBUG = DescriptorHandlerFlexophore.DEBUG;


    /**
     * Aromatic imide structure with an exocyclic N. Two N in the aromatic ring. Results in an extreme electron poor
     * exocyclic N. which is not making any interactions. Consequently, it is removed from the subgraph lists.
     * The electron poor N is the non-aromatic N in the fragment definitions!
     */
    public static final String IDCODE_EXO_N_AROM_IMIDE = "eMPARVCjK|X`";

    // Imide structure separated by one bond in the aromatic ring
    public static final String IDCODE_EXO_N_AROM_IMIDE_ALPHA = "gO|@AfeJih@PA@";

    public static String [] ARR_EXO_N_AROM_IMIDE = {IDCODE_EXO_N_AROM_IMIDE, IDCODE_EXO_N_AROM_IMIDE_ALPHA};


    /**
     * Similarity 0.977, similarity for identical molecule and a timeout of 5 min. So timeout of 6 min should be fine.
     *
     CreatorMolDistHistViz: ExceptionTimeOutConformerGeneration for idcode enY\JH@@amaNe`ZPICHhdhThdleEEEDhYThddZFKGLRX@J`@jjiijjZjAPbbT@@, hence generated 8 conformers.
     5 Minutes 52 Seconds 79 Millisec
     CreatorMolDistHistViz: ExceptionTimeOutConformerGeneration for idcode enY\JH@@amaNe`ZPICHhdhThdleEEEDhYThddZFKGLRX@J`@jjiijjZjAPbRR@@, hence generated 25 conformers.
     * 5 Minutes 52 Seconds 79 Millisec
     * [(390*2) (391) (392) (4358*6) (9,4358*6) (4358*4,4488) (590088,598407) (590088,598407)]
     */

    // public static final long TIMEOUT_CONFORMER_CALCULATION_MS = TimeDelta.MS_SECOND * 30;

    // Maximum number of tries to generate conformers with the torsion rule based conformer generator from Thomas Sander



    private static final int MAX_NUM_ATOMS = 1000;

    private static final int CONF_GEN_TS = 0;
    public static final int CONF_GIVEN_SINGLE_CONFORMATION = 1;
    public static final int SINGLE_CONFORMATION = 2;

    private static CreatorMolDistHistViz INSTANCE;


    private SubGraphExtractor subGraphExtractor;

    private int conformationMode;

    private long seed;


    // for debugging
    private boolean onlyOneConformer;

    private Exception recentException = null;

    // Calling SSSearcher frequently generates errors.
    // private SSSearcher ssSearcher;

    private StereoMolecule [] arrElectronPoorN;
    private ConformerGeneratorStageTries conformerGeneratorStageTries;

    public CreatorMolDistHistViz() {

        subGraphExtractor = new SubGraphExtractor();

        conformationMode = CONF_GEN_TS;

        conformerGeneratorStageTries = new ConformerGeneratorStageTries();

        IDCodeParser parser = new IDCodeParser();

        arrElectronPoorN = new StereoMolecule[ARR_EXO_N_AROM_IMIDE.length];

        for (int i = 0; i < ARR_EXO_N_AROM_IMIDE.length; i++) {
            StereoMolecule frag = parser.getCompactMolecule(ARR_EXO_N_AROM_IMIDE[i]);
            frag.ensureHelperArrays(Molecule.cHelperRings);
            arrElectronPoorN[i]=frag;
        }

    }

    public void setThreadMaster(ThreadMaster threadMaster) {
        conformerGeneratorStageTries.setThreadMaster(threadMaster);
    }

    public MolDistHistViz create(StereoMolecule molOrig) throws Exception {
        MolDistHistViz mdhv = null;
        switch (conformationMode) {
            case CONF_GEN_TS:
                mdhv = createMultipleConformations(molOrig, DescriptorHandlerFlexophore.NUM_CONFORMATIONS);
                break;
            case CONF_GIVEN_SINGLE_CONFORMATION:
                mdhv = createFromGivenConformation(molOrig);
                break;
            case SINGLE_CONFORMATION:
                mdhv = createMultipleConformations(molOrig, 1);
                break;
            default:
                throw new RuntimeException("Invalid conformation mode");
        }
        return mdhv;
    }

    public Exception getRecentException() {
        return recentException;
    }

    /**
     * If initializing with new molecule call resetInitializationStage() before!
     * @param molInPlace
     * @return
     */


    /**
     * Conformation generator of Thomas Sander
     * The molecule is standardized first.
     * @param molOrig
     * @return
     * @throws Exception
     */
    public MolDistHistViz createMultipleConformations(StereoMolecule molOrig, int nConformations) throws Exception {

        StereoMolecule molStand = molOrig.getCompactCopy();
        MoleculeStandardizer.standardize(molStand, MoleculeStandardizer.MODE_GET_PARENT);
        molStand.ensureHelperArrays(Molecule.cHelperRings);
        Molecule3D molInPlace = new Molecule3D(molStand);
        molInPlace.ensureHelperArrays(Molecule.cHelperRings);

        InteractionAtomTypeCalculator.setInteractionTypes(molInPlace);
        boolean successfulInitialization = conformerGeneratorStageTries.setMolecule(molInPlace);

        if(!successfulInitialization){
            return null;
        }

        List<SubGraphIndices> liSubGraphIndices = getSubGraphIndices(molInPlace);

        List<MultCoordFragIndex> liMultCoordFragIndex = new ArrayList<>();
        for (SubGraphIndices subGraphIndices : liSubGraphIndices) {
            liMultCoordFragIndex.add(new MultCoordFragIndex(subGraphIndices.getAtomIndices()));
        }

        Molecule3D molViz = createConformations(molInPlace, liMultCoordFragIndex, nConformations);
        int nPotentialConformers = conformerGeneratorStageTries.getPotentialConformerCount();
        onlyOneConformer = false;
        if((nPotentialConformers > 1) && (liMultCoordFragIndex.get(0).getCoordinates().size()==1)){

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
    public MolDistHistViz createFromConformerSet(ConformerSet conformerSet)  {

        Molecule3D mol = new Molecule3D(conformerSet.iterator().next().toMolecule());
        mol.ensureHelperArrays(Molecule.cHelperRings);
        InteractionAtomTypeCalculator.setInteractionTypes(mol);

        List<SubGraphIndices> liSGI = getSubGraphIndices(mol);
        List<MultCoordFragIndex> liMultCoordFragIndex = new ArrayList<>();
        for (SubGraphIndices subGraphIndices : liSGI) {
            liMultCoordFragIndex.add(new MultCoordFragIndex(subGraphIndices.getAtomIndices()));
        }

        Iterator<Conformer> it = conformerSet.iterator();
        while (it.hasNext()){
            Conformer c = it.next();
            Molecule3D molConf = new Molecule3D(c.toMolecule());
            molConf.ensureHelperArrays(Molecule.cHelperRings);
            calcFragmentCenter(molConf, liMultCoordFragIndex);
        }

        Molecule3D molViz = createPharmacophorePoints(mol, liMultCoordFragIndex);
        return create(liMultCoordFragIndex, molViz);
    }

    public List<SubGraphIndices> getSubGraphIndices(Molecule3D molInPlace){
        //
        // Handle carbon atoms connected to hetero atoms
        //
        List<SubGraphIndices> liSubGraphIndices = subGraphExtractor.extract(molInPlace);
        liSubGraphIndices = handleCarbonConnected2Hetero(liSubGraphIndices, molInPlace);
        liSubGraphIndices = removeExoCyclicElectronPoorN(liSubGraphIndices, molInPlace);

        return liSubGraphIndices;
    }

    public int getPotentialConformerCount(){
	    return conformerGeneratorStageTries.getPotentialConformerCount();
    }


    /**
     * This method must be called before:
     *  initializeConformers(molInPlace);
     *
     * Time in nanoseconds for a small molecule with idcode fegPb@JByH@QdbbbarTTbb^bRIRNQsjVZjjjh@J@@@
     * 75491600 first conformation
     * 16700 sec conformation
     * 37200 ...
     * 40100 ...
     *
     * @param molInPlace
     * @param liMultCoordFragIndex
     * @param nConformations
     * @return
     */
    public Molecule3D createConformations(Molecule3D molInPlace, List<MultCoordFragIndex> liMultCoordFragIndex, int nConformations) throws ExceptionConformationGenerationFailed  {
        int nAtoms = molInPlace.getAtoms();
        int ccConformationsGenerated = 0;
        Molecule3D molViz = null;

        for (int i = 0; i < nConformations; i++) {
            boolean conformerGenerated = false;
            try {
                conformerGenerated = conformerGeneratorStageTries.generateConformerAndSetCoordinates(nAtoms, molInPlace);
            } catch (ExceptionTimeOutConformerGeneration e) {
                System.err.println(
                    "CreatorMolDistHistViz: ExceptionTimeOutConformerGeneration for idcode " + molInPlace.getIDCode( )+ ", hence generated " + ccConformationsGenerated + " conformers.");
                // e.printStackTrace();
                break;
            }

            if(!conformerGenerated){
                break;
            }
            ccConformationsGenerated++;
            calcFragmentCenter(molInPlace, liMultCoordFragIndex);
            if(i==0){
                molViz = createPharmacophorePoints(molInPlace, liMultCoordFragIndex);
            }
        }

        if(ccConformationsGenerated==0){
            throw new ExceptionConformationGenerationFailed("Impossible to generate one conformer!");
        }

        return molViz;
    }

    /**
     *
     * @param liSubGraphIndices
     * @param molInPlace
     */
    public  static List<SubGraphIndices> handleCarbonConnected2Hetero(List<SubGraphIndices> liSubGraphIndices, StereoMolecule molInPlace){
        List<SubGraphIndices> liSubGraphIndicesProcessed = new ArrayList<>();
        for (SubGraphIndices sgi : liSubGraphIndices) {
            int [] arrIndexAtomFragment = sgi.getAtomIndices();
            HashSet<Integer> hsIndexAtom2Remove = new HashSet<>();
            for (int indexAtFrag : arrIndexAtomFragment) {
                if (ExtendedMoleculeFunctions.isCarbonConnected2Hetero(molInPlace, indexAtFrag)) {
                    // Is isolated carbon?
                    if(ExtendedMoleculeFunctions.isIsolatedCarbon(molInPlace, indexAtFrag, arrIndexAtomFragment)){
                        hsIndexAtom2Remove.add(indexAtFrag);
                    }
                }
            }
            SubGraphIndices sgiProcessed = new SubGraphIndices();
            if(hsIndexAtom2Remove.size()>0) {
                for (int indexAtFrag : arrIndexAtomFragment) {
                    if(!hsIndexAtom2Remove.contains(indexAtFrag)){
                        sgiProcessed.addIndex(indexAtFrag);
                    }
                }
            } else {
                sgiProcessed.addIndex(arrIndexAtomFragment);
            }
            if(sgiProcessed.getNumIndices()>0)
                liSubGraphIndicesProcessed.add(sgiProcessed);
        }
        return liSubGraphIndicesProcessed;
    }

    public List<SubGraphIndices> removeExoCyclicElectronPoorN(List<SubGraphIndices> liSubGraphIndices, StereoMolecule molInPlace){

        List<Integer> liElectronPoorN = getElectronPoorN(molInPlace);

        boolean [] arrMatchAtom = new boolean[molInPlace.getAtoms()];

        for (int indexAt : liElectronPoorN) {
            arrMatchAtom[indexAt]=true;
        }

        List<SubGraphIndices> liSubGraphIndicesProcessed = new ArrayList<>();
        for (SubGraphIndices sgi : liSubGraphIndices) {
            int[] arrIndexAtomFragment = sgi.getAtomIndices();
            HashSet<Integer> hsIndexAtom2Remove = new HashSet<>();
            for (int indexAtFrag : arrIndexAtomFragment) {
                if (arrMatchAtom[indexAtFrag]) {
                    hsIndexAtom2Remove.add(indexAtFrag);
                }
            }

            SubGraphIndices sgiProcessed = new SubGraphIndices();
            if (hsIndexAtom2Remove.size() > 0) {
                for (int indexAtFrag : arrIndexAtomFragment) {
                    if (!hsIndexAtom2Remove.contains(indexAtFrag)) {
                        sgiProcessed.addIndex(indexAtFrag);
                    }
                }
            } else {
                sgiProcessed.addIndex(arrIndexAtomFragment);
            }
            if (sgiProcessed.getNumIndices() > 0)
                liSubGraphIndicesProcessed.add(sgiProcessed);
        }
        return liSubGraphIndicesProcessed;
    }

    /**
     * The electron poor N is the non-aromatic N!
     * @param molInPlace
     * @return
     */
    private List<Integer> getElectronPoorN(StereoMolecule molInPlace){
        List<Integer> liElectronPoorN = new ArrayList<>();

        // Calling SSSearcher frequently generates errors.
        SSSearcher ssSearcher = new SSSearcher();
        ssSearcher.setMolecule(molInPlace);
        for (int i = 0; i < arrElectronPoorN.length; i++) {

            ssSearcher.setFragment(arrElectronPoorN[i]);
            int numFrags = ssSearcher.findFragmentInMolecule(
                    SSSearcher.cCountModeOverlapping,
                    SSSearcher.cMatchDBondToDelocalized | SSSearcher.cMatchAromDBondToDelocalized );
            if(numFrags>0) {
                // System.out.println("Found!");
                List<int[]> li = ssSearcher.getMatchList();
                for (int[] arrIndex : li) {
                    for (int j = 0; j < arrIndex.length; j++) {
                        int indAt = arrIndex[j];
                        int atNo = molInPlace.getAtomicNo(indAt);
                        if(atNo==7){
                            if(!molInPlace.isAromaticAtom(indAt)){
                                liElectronPoorN.add(indAt);
                            }
                        }
                    }
                }
            }
        }
        return liElectronPoorN;
    }


    /**
     * Creates the descriptor from the coordinates.
     * @param liMultCoordFragIndex contains the coordinates and the related atom indices of the molecule
     * @param molecule3D, must contain the interaction types.
     * @return
     */
    public static MolDistHistViz create(List<MultCoordFragIndex> liMultCoordFragIndex, Molecule3D molecule3D){

        MolDistHistViz molDistHistViz = new MolDistHistViz(liMultCoordFragIndex.size(), molecule3D);

        List<PPNodeViz> liPPNodeViz = new ArrayList<>();
        for (int i = 0; i < liMultCoordFragIndex.size(); i++) {
            MultCoordFragIndex mcfi = liMultCoordFragIndex.get(i);

            int [] arrIndexAtomFrag = mcfi.getArrIndexFrag();

            PPNodeViz ppNodeViz = new PPNodeViz();
            ppNodeViz.setIndex(i);
            ppNodeViz.setCoordinates(liMultCoordFragIndex.get(i).getCoordinates().get(0));

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

    public void addFlexophoreDistanceHistograms(MolDistHistViz mdhv){
        Molecule3D molecule = mdhv.getMolecule();

        List<MultCoordFragIndex> liMultCoordFragIndex = new ArrayList<>();
        for (PPNodeViz node : mdhv.getNodes()) {
            liMultCoordFragIndex.add(new MultCoordFragIndex(node.getArrayIndexOriginalAtoms()));
        }
        conformerGeneratorStageTries.setMolecule(molecule);
        int ccConformationsGenerated=0;
        for (int i = 0; i < DescriptorHandlerFlexophore.NUM_CONFORMATIONS; i++) {
            boolean conformerGenerated = false;
            try {
                conformerGenerated = conformerGeneratorStageTries.generateConformerAndSetCoordinates(molecule.getAtoms(), molecule);
            } catch (ExceptionTimeOutConformerGeneration e) {
                System.err.println(
                        "CreatorMolDistHistViz: ExceptionTimeOutConformerGeneration for idcode " + molecule.getIDCode( )+ ", hence generated " + ccConformationsGenerated + " conformers.");
                // e.printStackTrace();
                break;
            }

            if(!conformerGenerated){
                break;
            }
            ccConformationsGenerated++;
            calcFragmentCenter(molecule, liMultCoordFragIndex);

        }

        for (int i = 0; i < liMultCoordFragIndex.size(); i++) {
            for (int j = i+1; j < liMultCoordFragIndex.size(); j++) {
                byte [] arrDistHist = MultCoordFragIndex.getDistHist(liMultCoordFragIndex.get(i), liMultCoordFragIndex.get(j));
                // System.out.println(StringFunctions.toString(arrDistHist));
                mdhv.setDistHist(i,j,arrDistHist);
            }
        }

    }
    public MolDistHistViz createWithoutCoordinates(Molecule3D molecule3D){
        InteractionAtomTypeCalculator.setInteractionTypes(molecule3D);
        List<SubGraphIndices> sgis = getSubGraphIndices(molecule3D);
        return  createWithoutCoordinates(sgis, molecule3D);
    }
    public MolDistHistViz createWithoutCoordinates(Molecule3D molecule3D, int [] indices2consider){
        InteractionAtomTypeCalculator.setInteractionTypes(molecule3D);
        List<SubGraphIndices> sgis = getSubGraphIndices(molecule3D);
        List<SubGraphIndices> sgis2consider = new ArrayList<>();
        for (SubGraphIndices sgi : sgis) {
            if(ArrayUtilsCalc.containsAll(indices2consider, sgi.getAtomIndices()))    {
                sgis2consider.add(sgi);
            }
        }
        return  createWithoutCoordinates(sgis2consider, molecule3D);
    }


    public static MolDistHistViz createWithoutCoordinates(List<SubGraphIndices> liMultCoordFragIndex, Molecule3D molecule3D){

        MolDistHistViz molDistHistViz = new MolDistHistViz(liMultCoordFragIndex.size(), molecule3D);
        List<PPNodeViz> liPPNodeViz = new ArrayList<>();
        for (int i = 0; i < liMultCoordFragIndex.size(); i++) {
            SubGraphIndices sgi = liMultCoordFragIndex.get(i);

            int [] arrIndexAtomFrag = sgi.getAtomIndices();

            PPNodeViz ppNodeViz = new PPNodeViz();
            ppNodeViz.setIndex(i);

            for (int index : arrIndexAtomFrag) {
                int interactionType = molecule3D.getInteractionAtomType(index);
                ppNodeViz.add(interactionType);
                ppNodeViz.addIndexOriginalAtom(index);
            }
            liPPNodeViz.add(ppNodeViz);
        }

        molDistHistViz.set(liPPNodeViz);

        byte [] arrHistPercent = new byte [ConstantsFlexophoreGenerator.BINS_HISTOGRAM];
        Arrays.fill(arrHistPercent, (byte)1);

        for (int i = 0; i < liMultCoordFragIndex.size(); i++) {
            for (int j = i+1; j < liMultCoordFragIndex.size(); j++) {
                molDistHistViz.setDistHist(i,j,arrHistPercent);
            }
        }
        molDistHistViz.realize();
        return molDistHistViz;
    }

    private static Molecule3D createPharmacophorePoints(Molecule3D molecule3D, List<MultCoordFragIndex> liMultCoordFragIndex) {

        Molecule3D molCenter = new Molecule3D(molecule3D);
        molCenter.ensureHelperArrays(Molecule.cHelperRings);
        for (MultCoordFragIndex multCoordFragIndex : liMultCoordFragIndex) {
            int [] arrAtomIndexList = multCoordFragIndex.getArrIndexFrag();
            // Calculate center coordinates.
            Coordinates coordCenter = molCenter.getCenterOfGravity(arrAtomIndexList);
            for (int at = 0; at < arrAtomIndexList.length; at++) {
                // Atom type.
                int interactionType = molCenter.getInteractionAtomType(arrAtomIndexList[at]);
                if(interactionType == ConstantsFlexophoreGenerator.INTERACTION_TYPE_NONE){
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
                molCenter.setAtomFlag(indexAtm, ConstantsFlexophore.FLAG_CENTER_ATOM, true);
                molCenter.setPPP(indexAtm, arrAtomIndexList);
            }
        }
        molCenter.ensureHelperArrays(Molecule.cHelperRings);
        return molCenter;
    }




    /**
     * Use valid only with createMultipleConformations(StereoMolecule molOrig, int nConformations)
     * @return
     */
    public boolean isOnlyOneConformer() {
        return onlyOneConformer;
    }


    public MolDistHistViz createFromGivenConformation(StereoMolecule molOrig) {

        // Just to make sure we will not change something in the original molecule.
        Molecule3D molStart = new Molecule3D(molOrig);
        molStart.ensureHelperArrays(Molecule.cHelperRings);
        InteractionAtomTypeCalculator.setInteractionTypes(molStart);
        Molecule3D molInPlace = new Molecule3D(molStart);
        molInPlace.ensureHelperArrays(Molecule.cHelperRings);
        List<SubGraphIndices> liSubGraphIndices = getSubGraphIndices(molInPlace);
        List<MultCoordFragIndex> liMultCoordFragIndex = new ArrayList<>();
        for (SubGraphIndices subGraphIndices : liSubGraphIndices) {
            liMultCoordFragIndex.add(new MultCoordFragIndex(subGraphIndices.getAtomIndices()));
        }
        calcFragmentCenter(molStart, liMultCoordFragIndex);
        MolDistHistViz mdhv = create(liMultCoordFragIndex, molStart);
        return mdhv;
    }

    /**
     * Calculates the center of the fragments and stores the coordinates in MultCoordFragIndex.
     * @param molecule3D
     * @param liMultCoordFragIndex, the coordinates are added to this list.
     */
    public static void calcFragmentCenter(Molecule3D molecule3D, List<MultCoordFragIndex> liMultCoordFragIndex) {

        for (MultCoordFragIndex multCoordFragIndex : liMultCoordFragIndex) {
            int [] arrAtomIndexList = multCoordFragIndex.getArrIndexFrag();
            // Calculate center coordinates.
            Coordinates coordCenter = molecule3D.getCenterOfGravity(arrAtomIndexList);
            multCoordFragIndex.addCoord(coordCenter);
        }
    }

    public static CreatorMolDistHistViz getInstance(){
        if(INSTANCE==null){
            INSTANCE = new CreatorMolDistHistViz();
        }
        return INSTANCE;
    }
}
