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

import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.*;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.descriptor.DescriptorHandlerFlexophore;
import com.actelion.research.chem.descriptor.flexophore.*;
import com.actelion.research.chem.descriptor.flexophore.redgraph.SubGraphExtractor;
import com.actelion.research.chem.descriptor.flexophore.redgraph.SubGraphIndices;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import org.openmolecules.chem.conf.gen.ConformerGenerator;
import org.openmolecules.chem.conf.gen.RigidFragmentCache;

import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;

/**
 * CreatorMolDistHistViz
 * Created by korffmo1 on 19.02.16.
 */
public class CreatorMolDistHistViz {

    private static final boolean DEBUG = DescriptorHandlerFlexophore.DEBUG;

    public static final long SEED = 123456789;

    // Maximum number of tries to generate conformers with the torsion rule based conformer generator from Thomas Sander
    private static final int MAX_NUM_TRIES = 10000;

    private static final int MAX_NUM_ATOMS = 1000;

    private static final int CONF_GEN_TS = 0;

    public static final int CONF_GIVEN_SINGLE_CONFORMATION = 1;
    public static final int SINGLE_CONFORMATION = 2;

    private static CreatorMolDistHistViz INSTANCE;


    private SubGraphExtractor subGraphExtractor;

    private ConformerGenerator conformerGenerator;

    private int conformationMode;

    private long seed;

    // for debugging
    private boolean onlyOneConformer;

    private int [] arrIndexAtomNewTmp;

    public CreatorMolDistHistViz() {

        seed = SEED;

        subGraphExtractor = new SubGraphExtractor();

        conformerGenerator = new ConformerGenerator(seed, false);
        RigidFragmentCache.getDefaultInstance().loadDefaultCache();

        conformationMode = CONF_GEN_TS;

        arrIndexAtomNewTmp = new int[MAX_NUM_ATOMS];

        // System.out.println("CreatorCompleteGraph conformationMode " + conformationMode);

    }

    public void setThreadMaster(ThreadMaster threadMaster) {
        conformerGenerator.setThreadMaster(threadMaster);
    }

    public void setConformationMode(int conformationMode) {
        this.conformationMode = conformationMode;
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


    /**
     * Conformation generator of Thomas Sander
     * The molecule is standardized first.
     * @param molOrig
     * @return
     * @throws Exception
     */
    public MolDistHistViz createMultipleConformations(StereoMolecule molOrig, int nConformations) throws Exception {

        // int nConformations = DescriptorHandlerFlexophore.NUM_CONFORMATIONS;

        StereoMolecule molStand = molOrig.getCompactCopy();

        MoleculeStandardizer.standardize(molStand, MoleculeStandardizer.MODE_GET_PARENT);

        molStand.ensureHelperArrays(Molecule.cHelperRings);

        Molecule3D molInPlace = new Molecule3D(molStand);

        molInPlace.ensureHelperArrays(Molecule.cHelperRings);

        conformerGenerator.initializeConformers(molInPlace, ConformerGenerator.STRATEGY_LIKELY_RANDOM, MAX_NUM_TRIES, false);
        
        InteractionAtomTypeCalculator.setInteractionTypes(molInPlace);

        //
        // Handle carbon atoms connected to hetero atoms
        //
        List<SubGraphIndices> liSubGraphIndices = subGraphExtractor.extract(molInPlace);
        liSubGraphIndices = handleCarbonConnected2Hetero(liSubGraphIndices, molInPlace);

        if(DEBUG) {
            injectNewSeed();
        }

        int nAtoms = molInPlace.getAtoms();

        List<MultCoordFragIndex> liMultCoordFragIndex = new ArrayList<>();
        for (SubGraphIndices subGraphIndices : liSubGraphIndices) {
            liMultCoordFragIndex.add(new MultCoordFragIndex(subGraphIndices.getAtomIndices()));
        }

        Molecule3D molViz = createConformations(molInPlace, liMultCoordFragIndex, nConformations, conformerGenerator);

        int nPotentialConformers = conformerGenerator.getPotentialConformerCount();

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


    /**
     * This method must be called before:
     *  conformerGenerator.initializeConformers(molInPlace, ConformerGenerator.STRATEGY_LIKELY_RANDOM, MAX_NUM_TRIES, false);
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
     * @param conformerGenerator
     * @return
     */
    public static Molecule3D createConformations(Molecule3D molInPlace, List<MultCoordFragIndex> liMultCoordFragIndex, int nConformations, ConformerGenerator conformerGenerator){
        int nAtoms = molInPlace.getAtoms();
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

    /**
     * Creates the descriptor from the coordinates.
     * @param liMultCoordFragIndex contains the ccordinates and the related atom indices of the molecule
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

    private static Molecule3D createPharmacophorePoints(Molecule3D molecule3D, List<MultCoordFragIndex> liMultCoordFragIndex) {

        Molecule3D molCenter = new Molecule3D(molecule3D);

        molCenter.ensureHelperArrays(Molecule.cHelperRings);

        for (MultCoordFragIndex multCoordFragIndex : liMultCoordFragIndex) {

            int [] arrAtomIndexList = multCoordFragIndex.getArrIndexFrag();

            // Calculate center coordinates.
            Coordinates coordCenter = ExtendedMolecule.getCenterGravity(molCenter, arrAtomIndexList);

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

        List<SubGraphIndices> liSubGraphIndices = subGraphExtractor.extract(molInPlace);
        liSubGraphIndices = handleCarbonConnected2Hetero(liSubGraphIndices, molInPlace);

        List<MultCoordFragIndex> liMultCoordFragIndex = new ArrayList<>();
        for (SubGraphIndices subGraphIndices : liSubGraphIndices) {
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
    public static void calcFragmentCenter(Molecule3D molecule3D, List<MultCoordFragIndex> liMultCoordFragIndex) {

        for (MultCoordFragIndex multCoordFragIndex : liMultCoordFragIndex) {

            int [] arrAtomIndexList = multCoordFragIndex.getArrIndexFrag();

            // Calculate center coordinates.
            Coordinates coordCenter = ExtendedMolecule.getCenterGravity(molecule3D, arrAtomIndexList);
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
