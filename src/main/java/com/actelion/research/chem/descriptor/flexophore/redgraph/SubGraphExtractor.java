package com.actelion.research.chem.descriptor.flexophore.redgraph;

import com.actelion.research.chem.*;
import com.actelion.research.util.hash.HashSetInt;

import java.util.*;


/**
 * SubGraphExtractor
 *
 * Creates substructures by fragmentation and summarizes them.
 *
 * <p>Copyright: Actelion Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * Created by korffmo1 on 17.02.16.
 * 10.03.2017 Bug fix in isEnclosingRing() for multi-bridged ring systems.
 */
public class SubGraphExtractor {

    private static final int MAX_NUM_ATOMS = 10000;

    private static final int MAX_RING_SIZE_TO_SUMMARIZE_HETERO_RINGS = 6;

    private static final int MAX_DEPTH_END_STANDING_ALIPHATIC_GROUP = 1;

    private static final int MIN_SIZE_ALIPHATIC_CHAIN_IN_RING = 3;

    private static final int DISTANCE_THRESH_SPLIT_GROUP = 6;

    private static final int GROUP_SIZE_SPLITTED = 3;


    private boolean [] arrAtomInSmallRing;

    private boolean [] arrAtomInLargeRing;

    private boolean [] arrRingAtom;

    private HashSet<Integer> hsAtomicNumberExcludeAliphatic;


    public SubGraphExtractor() {
        arrRingAtom = new boolean[MAX_NUM_ATOMS];

        hsAtomicNumberExcludeAliphatic = new HashSet<>();
        hsAtomicNumberExcludeAliphatic.add(7);
        hsAtomicNumberExcludeAliphatic.add(8);
    }

    /**
     * Delivers a list of indices that describes sub structures of the given molecule.
     * @param molOrig
     * @return
     */
    public List<SubGraphIndices> extract(StereoMolecule molOrig) {

        int nAtoms = molOrig.getAtoms();

        for (int i = 0; i < nAtoms; i++) {
            arrRingAtom[i] = false;
        }

        StereoMolecule mol = new StereoMolecule(molOrig);

        mol.ensureHelperArrays(Molecule.cHelperRings);

        createRingAtomIndexMap(mol);

        //
        // Detect end standing atoms
        //
        List<Integer> liEndStandingAtoms = new ArrayList<>();
        for (int i = 0; i < nAtoms; i++) {

            int nConnAtms = mol.getAllConnAtoms(i);
            int ccConnNonHydrogen = 0;
            for (int j = 0; j < nConnAtms; j++) {
                int indexConnAtm = mol.getConnAtom(i,j);
                int atNo = mol.getAtomicNo(indexConnAtm);
                if(atNo>1){
                    ccConnNonHydrogen++;
                }
            }

            if(ccConnNonHydrogen==1){
                liEndStandingAtoms.add(i);
            }
        }

        List<SubGraphIndices> liFragment = new ArrayList<>();

        //
        // End standing hetero atom groups
        //
        List<SubGraphIndices> liFragmentHetero = getEndStandingHeteroGroups(mol, liEndStandingAtoms);

        liFragment.addAll(liFragmentHetero);

        HashSetInt hsAtomIndicesInFragment = new HashSetInt();

        SubGraphIndices.addAtomIndices(hsAtomIndicesInFragment, liFragmentHetero);

        //
        // Small rings
        //
        List<SubGraphIndices> liFragmentRings = getSmallRings(mol);

        liFragment.addAll(liFragmentRings);

        HashSetInt hsAtomIndicesInSmallRings = new HashSetInt();

        SubGraphIndices.addAtomIndices(hsAtomIndicesInSmallRings, liFragmentRings);

        hsAtomIndicesInFragment.add(hsAtomIndicesInSmallRings.getValues());

        //
        // Remaining hetero atoms
        // Hetero atoms that were not covered by the end standing groups and that are not enclosed in the small rings.
        //
        List<SubGraphIndices> liFragmentRemainingHetero = getRemainingHeteroGroups(mol, hsAtomIndicesInFragment);

        liFragment.addAll(liFragmentRemainingHetero);

        SubGraphIndices.addAtomIndices(hsAtomIndicesInFragment, liFragmentRemainingHetero);

        //
        // End standing aliphatic group
        //
        List<SubGraphIndices> liFragmentEndStandingAliphaticGroup = getEndStandingAliphaticGroups(mol, liEndStandingAtoms, hsAtomIndicesInFragment);

        liFragment.addAll(liFragmentEndStandingAliphaticGroup);

        HashSetInt hsAtomIndicesInFragmentEndStandingAliphaticGroup = new HashSetInt();

        SubGraphIndices.addAtomIndices(hsAtomIndicesInFragmentEndStandingAliphaticGroup, liFragmentEndStandingAliphaticGroup);

        //
        // Aliphatic groups in large rings.
        //
        List<SubGraphIndices> liFragmentLargeRings = getAliphaticGroupsInLargeRings(mol, hsAtomIndicesInFragmentEndStandingAliphaticGroup);

        liFragment.addAll(liFragmentLargeRings);

        //
        // Merge bridged aliphatic rings
        // Only in place for heavily bridged rings.
        mergeBridgedAliphaticRings(mol, liFragment);

        //
        // Split long hetero groups into subgroups
        //
        splitLongHeteroGroups(mol, liFragment);

        return  liFragment;
    }

    public List<SubGraphIndices> extractAliphaticRingsAndEndStandingAliphaticGroups(StereoMolecule molOrig) {

        StereoMolecule mol = new StereoMolecule(molOrig);

        mol.ensureHelperArrays(Molecule.cHelperRings);

        int nAtoms = mol.getAtoms();

        createRingAtomIndexMap(mol);

        //
        // Detect end standing atoms
        //
        List<Integer> liEndStandingAtoms = new ArrayList<>();
        for (int i = 0; i < nAtoms; i++) {

            int nConnAtms = mol.getAllConnAtoms(i);
            int ccConnNonHydrogen = 0;
            for (int j = 0; j < nConnAtms; j++) {
                int indexConnAtm = mol.getConnAtom(i,j);
                int atNo = mol.getAtomicNo(indexConnAtm);
                if(atNo>1){
                    ccConnNonHydrogen++;
                }
            }

            if(ccConnNonHydrogen==1){
                liEndStandingAtoms.add(i);
            }
        }

        HashSetInt hsAtomNonAliphatic = new HashSetInt();
        for (int indexEndStandingAtom : liEndStandingAtoms) {
            if(ExtendedMoleculeFunctions.isHetero(mol, indexEndStandingAtom)){
                hsAtomNonAliphatic.add(indexEndStandingAtom);
            }
        }

        List<SubGraphIndices> liFragment = new ArrayList<>();

        //
        // Small aliphatic rings
        //
        List<SubGraphIndices> liSGIRings = getSmallRings(mol);
        for (SubGraphIndices sgiRing : liSGIRings) {

            int [] arrIndAt =  sgiRing.getAtomIndices();

            int ccAromatic = 0;
            for (int indAt : arrIndAt) {
                if(mol.isAromaticAtom(indAt)){
                    ccAromatic++;
                }
            }

            if(ccAromatic/arrIndAt.length<0.5){
                if(containsOnlyCarbon(mol, sgiRing)){
                    liFragment.add(sgiRing);
                }
            }
        }

        //
        // End standing aliphatic group
        //
        List<SubGraphIndices> liFragmentEndStandingAliphaticGroup = getEndStandingAliphaticGroups(mol, liEndStandingAtoms, hsAtomNonAliphatic);

        liFragment.addAll(liFragmentEndStandingAliphaticGroup);

        return  liFragment;
    }

    /**
     * Splits long hetero atom chains into smaller sub groups.
     * @param mol
     * @param liFragment
     */
    private void splitLongHeteroGroups(StereoMolecule mol, List<SubGraphIndices> liFragment){

        final int nAtomsSubGroup = GROUP_SIZE_SPLITTED;

        int [][] arrTopologicalDistanceMatrix = ExtendedMoleculeFunctions.getTopologicalDistanceMatrix(mol);

        for (int i = liFragment.size()-1; i >= 0; i--) {

            SubGraphIndices sgi = liFragment.get(i);

            int [] arrIndAtm = sgi.getAtomIndices();

            MaximumTopologicalDist maxTopoDist = getMaximumDistance(arrTopologicalDistanceMatrix, arrIndAtm);

            if(maxTopoDist.dist >= DISTANCE_THRESH_SPLIT_GROUP){

                liFragment.remove(i);

                int [] arrIndAtmPathMaxDist = new int[maxTopoDist.dist+1];

                mol.getPath(arrIndAtmPathMaxDist, maxTopoDist.at1, maxTopoDist.at2, maxTopoDist.dist+1, null);

                boolean ringAtm = false;

                for (int indAtmPath : arrIndAtmPathMaxDist) {

                    if(mol.isRingAtom(indAtmPath)){

                        ringAtm = true;

                        break;

                    }

                }

                if(ringAtm){

                    continue;

                }

                HashSetInt hsIndAtmSubPath = new HashSetInt(arrIndAtmPathMaxDist);

                int nAtomsInPath = arrIndAtmPathMaxDist.length;

                for (int j = 0; j < nAtomsInPath; j+= nAtomsSubGroup) {

                    int size = Math.min(nAtomsSubGroup, nAtomsInPath-j);

                    int [] arrAtmSubPath = new int[size];

                    for (int k = 0; k < size; k++) {

                        int index = j+k;

                        arrAtmSubPath[k] = arrIndAtmPathMaxDist[index];
                    }


                    SubGraphIndices sgiSub = new SubGraphIndices();

                    sgiSub.addIndex(arrAtmSubPath);

                    addConnAtoms(mol, sgiSub, sgi, hsIndAtmSubPath);

                    liFragment.add(sgiSub);

                }
            }
        }
    }


    private static void addConnAtoms(StereoMolecule mol, SubGraphIndices sgiSub, SubGraphIndices sgi, HashSetInt hsIndAtmSubPath){

        int [] arrAtmSubPath = sgiSub.getAtomIndices();


        LinkedList<Integer> queue = new LinkedList<>();

        for (int indAtm : arrAtmSubPath) {

            queue.add(indAtm);

        }

        while (!queue.isEmpty()){

            int indAtm = queue.poll();

            int nConn = mol.getConnAtoms(indAtm);

            for (int i = 0; i < nConn; i++) {

                int indAtmConn = mol.getConnAtom(indAtm, i);

                if(hsIndAtmSubPath.contains(indAtmConn)){

                    continue;

                }

                if(sgiSub.contains(indAtmConn)){

                    continue;

                }

                if(!sgi.contains(indAtmConn)) {

                    continue;

                }

                sgiSub.addIndex(indAtmConn);

                queue.add(indAtmConn);

            }

        }

    }

    /**
     * Maximum topological distance in bonds.
     * @param arrTopologicalDistanceMatrix
     * @param arrIndAtm
     * @return
     */

    private static MaximumTopologicalDist getMaximumDistance(int [][] arrTopologicalDistanceMatrix, int [] arrIndAtm) {

        MaximumTopologicalDist maxTopoDist = new MaximumTopologicalDist();

        for (int i = 0; i < arrIndAtm.length; i++) {

            for (int j = i+1; j < arrIndAtm.length; j++) {

                int d = arrTopologicalDistanceMatrix[arrIndAtm[i]][arrIndAtm[j]];

                if(d > maxTopoDist.dist){

                    maxTopoDist.dist = d;

                    maxTopoDist.at1 = arrIndAtm[i];

                    maxTopoDist.at2 = arrIndAtm[j];
                }
            }
        }

        return maxTopoDist;
    }


    private static void mergeBridgedAliphaticRings(StereoMolecule mol, List<SubGraphIndices> liFragment){

        // If set to two also condensed rings will be merged.
        final int minNumOverLappingIndices = 3;

        boolean merged = true;

        RingCollection ringCollection = mol.getRingSet();

        while (merged){

            merged = false;

            mergedB:
            for (int i = 0; i < liFragment.size(); i++) {

                SubGraphIndices sgi0 = liFragment.get(i);

                if(!containsOnlyCarbon(mol, sgi0)){

                    continue;

                }

                for (int j = i+1; j < liFragment.size(); j++) {

                    SubGraphIndices sgi1 = liFragment.get(j);

                    if(!containsOnlyCarbon(mol, sgi1)){

                        continue;

                    }

                    int nOverLapRings = getNumOverlappingRingIndices(mol, ringCollection, sgi0, sgi1);

                    if(nOverLapRings >= minNumOverLappingIndices){

                        sgi0.merge(sgi1);

                        liFragment.remove(j);

                        merged = true;

                        break mergedB;
                    }
                }
            }
        }
    }

    private static int getNumOverlappingRingIndices(StereoMolecule mol, RingCollection ringCollection, SubGraphIndices sgi1, SubGraphIndices sgi2){

        int maxOverlap = 0;

        int rings = ringCollection.getSize();

        for (int i = 0; i < rings; i++) {

            int [] arrIndAtRing1 = ringCollection.getRingAtoms(i);

            HashSetInt hs1 = new HashSetInt();

            for (int indAtRing1 : arrIndAtRing1) {

                if(sgi1.contains(indAtRing1)){

                    hs1.add(indAtRing1);

                }
            }

            int [] arrIndAtHs1 = hs1.getValues();

            for (int j = 0; j < rings; j++) {

                int[] arrIndAtRing2 = ringCollection.getRingAtoms(j);

                HashSetInt hs2 = new HashSetInt();

                for (int indAtRing2 : arrIndAtRing2) {

                    if (sgi2.contains(indAtRing2)) {

                        hs2.add(indAtRing2);

                    }
                }

                int nOverlap = 0;

                for (int indAtHs1 : arrIndAtHs1) {

                    if(hs2.contains(indAtHs1)){

                        nOverlap++;
                    }

                }

                if(nOverlap > maxOverlap){

                    maxOverlap = nOverlap;
                }
            }
        }

        return maxOverlap;

    }


    private static boolean containsOnlyCarbon(StereoMolecule mol, SubGraphIndices sgi) {

        boolean carbon = true;

        int [] arrAtInd = sgi.getAtomIndices();

        for (int indAt : arrAtInd) {

            if(mol.getAtomicNo(indAt) != 6) {

                carbon = false;

                break;

            }
        }

        return carbon;
    }




    /**
     * Gets the remaining hetero groups that are in a chain without exo-atoms. like -S- , -O-, -N=-.
     * @param mol
     * @param hsAtomIndicesUsed
     * @return
     */
    private List<SubGraphIndices> getRemainingHeteroGroups(StereoMolecule mol, HashSetInt hsAtomIndicesUsed){

        List<SubGraphIndices> liFragmentRemainingHetero = new ArrayList<>();

        int atoms = mol.getAtoms();

        boolean [] arrAtomIndicesUsedMap = new boolean[atoms];

        int [] arrAtomIndicesUsed = hsAtomIndicesUsed.getValues();

        for (int indexAtmUsed : arrAtomIndicesUsed) {
            arrAtomIndicesUsedMap[indexAtmUsed] = true;
        }

        //
        // First layer hetero atoms
        //
        for (int indAtm = 0; indAtm < atoms; indAtm++) {

            if(arrAtomIndicesUsedMap[indAtm]){
                continue;
            }

            if(!ExtendedMoleculeFunctions.isHetero(mol, indAtm)){
                continue;
            }

            if(ExtendedMoleculeFunctions.isEtherOxygenAtAromatic(mol, indAtm)){
                continue;
            }

            SubGraphIndices fragment = new SubGraphIndices();

            fragment.addIndex(indAtm);

            liFragmentRemainingHetero.add(fragment);

        }

        //
        // Second layer
        //
        for (SubGraphIndices fragment : liFragmentRemainingHetero) {

            int indexAtmHetero1 = fragment.getAtomIndices()[0];

            int nAtmsConn = mol.getConnAtoms(indexAtmHetero1);

            for (int i = 0; i < nAtmsConn; i++) {

                int indexAtmConn = mol.getConnAtom(indexAtmHetero1, i);

                if(arrAtomIndicesUsedMap[indexAtmConn]){

                    continue;

                }

                if(ExtendedMoleculeFunctions.isHetero(mol, indexAtmConn)) {

                    fragment.addIndex(indexAtmConn);

                } else if(mol.getConnAtoms(indexAtmConn)==1) { // Add end standing carbon atoms.

                    fragment.addIndex(indexAtmConn);

                }
            }
        }

        SubGraphIndices.merge(liFragmentRemainingHetero);

        return liFragmentRemainingHetero;

    }

    /**
     * Also takes carbon atoms from large rings.
     * @param mol
     * @param liEndStandingAtoms
     * @param hsAtomIndicesUsed
     * @return
     */
    private List<SubGraphIndices> getEndStandingAliphaticGroups(StereoMolecule mol, List<Integer> liEndStandingAtoms, HashSetInt hsAtomIndicesUsed){

        List<SubGraphIndices> liFragmentEndStandingAliphaticGroup = new ArrayList<>();

        int atoms = mol.getAtoms();

        boolean [] arrAtomIndicesUsedMap = new boolean[atoms];

        int [] arrAtomIndicesUsed = hsAtomIndicesUsed.getValues();

        for (int indexAtmUsed : arrAtomIndicesUsed) {
            arrAtomIndicesUsedMap[indexAtmUsed] = true;
        }

        //
        // First layer end standing carbon atoms
        //

        for (int indexEndStandingAtom : liEndStandingAtoms) {
            if(arrAtomIndicesUsedMap[indexEndStandingAtom]){
                continue;
            }

            if(ExtendedMoleculeFunctions.isHetero(mol, indexEndStandingAtom)){ // All end standing hetero atoms should already be used here.
                throw new RuntimeException("This should not happen.");
            }

            if(areAtomicNoConnectedInList(mol, indexEndStandingAtom, hsAtomicNumberExcludeAliphatic)){
                continue;
            }

            int nConnected = mol.getConnAtoms(indexEndStandingAtom);

            boolean addEndStanding=true;

            for (int i = 0; i < nConnected; i++) {
                int indAtConn = mol.getConnAtom(indexEndStandingAtom, i);

                if(mol.isRingAtom(indAtConn)){
                    addEndStanding=false;
                    break;
                }

                if(areAtomicNoConnectedInList(mol, indAtConn, hsAtomicNumberExcludeAliphatic)){
                    addEndStanding=false;
                    break;
                }
            }


            if(addEndStanding) {
                SubGraphIndices fragment = new SubGraphIndices();
                fragment.addIndex(indexEndStandingAtom);
                liFragmentEndStandingAliphaticGroup.add(fragment);
            }



        }

        //
        // Next layers
        //
        for (SubGraphIndices fragment : liFragmentEndStandingAliphaticGroup) {
            broadFirstForNonRingCarbon(mol, fragment, hsAtomIndicesUsed, MAX_DEPTH_END_STANDING_ALIPHATIC_GROUP);
        }

        SubGraphIndices.merge(liFragmentEndStandingAliphaticGroup);

        return liFragmentEndStandingAliphaticGroup;

    }

    private static boolean areAtomicNoConnectedInList(StereoMolecule mol, int indexAtCenter, HashSet<Integer> hsAtomicNumber) {

        int nConnected = mol.getConnAtoms(indexAtCenter);

        boolean contains=false;

        for (int i = 0; i < nConnected; i++) {
            int indAtConn = mol.getConnAtom(indexAtCenter, i);

            int atomicNoConn = mol.getAtomicNo(indAtConn);
            if(hsAtomicNumber.contains(atomicNoConn)){
                contains=true;
                break;
            }
        }

        return contains;
    }

    /**
     * Broad first search to add carbon atoms to the end standing carbon given with the fragment.
     * @param mol
     * @param fragment
     * @param hsAtomIndicesUsed
     * @param maxDepth
     */
    private void broadFirstForNonRingCarbon(StereoMolecule mol, SubGraphIndices fragment, HashSetInt hsAtomIndicesUsed, int maxDepth){

        int atoms = mol.getAtoms();

        boolean [] arrAtomIndicesUsedMap = new boolean[atoms];

        int [] arrAtomIndicesUsed = hsAtomIndicesUsed.getValues();

        for (int indexAtmUsed : arrAtomIndicesUsed) {

            arrAtomIndicesUsedMap[indexAtmUsed] = true;

        }

        LinkedList<Integer> liIndAtm = new LinkedList<>();

        int [] arrIndAtmLayer = new int[atoms];

        // Contains only one index.
        int [] arrIndAtm = fragment.getAtomIndices();

        for (int indAtm : arrIndAtm) {
            liIndAtm.add(indAtm);
            arrIndAtmLayer[indAtm] = 0;
        }

        while (!liIndAtm.isEmpty()){

            int indAtm = liIndAtm.poll();

            if(arrIndAtmLayer[indAtm]>maxDepth){

                continue;

            }

            fragment.addIndex(indAtm);

            int nConnAtms = mol.getConnAtoms(indAtm);

            for (int i = 0; i < nConnAtms; i++) {

                int indAtmConn = mol.getConnAtom(indAtm, i);

                if(arrAtomIndicesUsedMap[indAtmConn]){

                    continue;

                }

                if(arrAtomInSmallRing[indAtmConn]) {

                    continue;

                }

                if(mol.getAtomicNo(indAtmConn)==6) {

                    liIndAtm.add(indAtmConn);

                    arrIndAtmLayer[indAtmConn] = arrIndAtmLayer[indAtm] + 1;

                }
            }
        }
    }

    /**
     * Checks for electronegative neighbors. Carbon with N, O or F as neighbors do not count for aliphatic.
     * @param mol
     * @param arrIndAtm
     * @return null if the indices refer only to non carbon. Also null if the minimum
     */
    private SubGraphIndices largestAliphaticFragment(StereoMolecule mol, int [] arrIndAtm, int minSizeFragment){

        List<SubGraphIndices> liFragment = new ArrayList<>();

        int atoms = mol.getAtoms();

        boolean [] arrIndAtmMap = new boolean[atoms];

        for (int indAtm : arrIndAtm) {

            arrIndAtmMap[indAtm] = true;

        }

        boolean [] arrIndAtmUsedMap = new boolean[atoms];

        LinkedList<Integer> liIndAtm = new LinkedList<>();

        boolean atomFound = true;

        while(atomFound) {

            atomFound = false;

            liIndAtm.clear();

            for (int indAtm : arrIndAtm) {

                if(!arrIndAtmUsedMap[indAtm]){

                    if(mol.getAtomicNo(indAtm)==6) {

                        liIndAtm.add(indAtm);

                        atomFound = true;

                        break;
                    }
                }
            }

            SubGraphIndices fragment = new SubGraphIndices();

            while (!liIndAtm.isEmpty()){

                int indAtm = liIndAtm.poll();

                fragment.addIndex(indAtm);

                arrIndAtmUsedMap[indAtm] = true;

                int nConnAtms = mol.getConnAtoms(indAtm);

                for (int i = 0; i < nConnAtms; i++) {

                    int indAtmConn = mol.getConnAtom(indAtm, i);

                    if(arrIndAtmUsedMap[indAtmConn]){

                        continue;

                    }

                    if(!arrIndAtmMap[indAtmConn]){

                        continue;

                    }

                    if(mol.getAtomicNo(indAtmConn)==6) {

                        liIndAtm.add(indAtmConn);

                    }
                }
            }

            liFragment.add(fragment);
        }

        //
        // Remove fragments too small or with to many electronegative neighbors.
        //
        for (int i = liFragment.size()-1; i >= 0; i--) {

            SubGraphIndices fragment = liFragment.get(i);

            if(fragment.getNumIndices()  < minSizeFragment) {

                liFragment.remove(i);

                continue;

            }

            int [] arrIndAtmFrag = fragment.getAtomIndices();

            int nCarbonWithElectronegativeNeighbors = 0;

            for (int indAtmFrag : arrIndAtmFrag) {

                int nAtmConn = mol.getConnAtoms(indAtmFrag);

                for (int j = 0; j < nAtmConn; j++) {

                    int indAtmConn = mol.getConnAtom(indAtmFrag, j);

                    int atomicNo = mol.getAtomicNo(indAtmConn);

                    if((atomicNo == PeriodicTable.Nitrogen) ||
                            (atomicNo == PeriodicTable.Oxygen)  ||
                            (atomicNo == PeriodicTable.Fluorine)) {

                        nCarbonWithElectronegativeNeighbors++;

                        break;

                    }
                }
            }

            if(arrIndAtmFrag.length - nCarbonWithElectronegativeNeighbors < minSizeFragment) {

                liFragment.remove(i);

            }
        }

        SubGraphIndices fragment = null;

        if(liFragment.size() > 0){

            Collections.sort(liFragment, SubGraphIndices.getComparatorNumIndices());

            fragment = liFragment.get(liFragment.size()-1);

        }

        return fragment;
    }


    private void createRingAtomIndexMap(StereoMolecule mol) {

        RingCollection ringCollection = mol.getRingSet();

        int rings = ringCollection.getSize();

        arrAtomInSmallRing = new boolean[mol.getAtoms()];

        arrAtomInLargeRing = new boolean[mol.getAtoms()];

        for (int ringNo = 0; ringNo < rings; ringNo++) {

            int ringSize = ringCollection.getRingSize(ringNo);

            int [] arrIndexRingAtoms = ringCollection.getRingAtoms(ringNo);

            if (ringSize <= MAX_RING_SIZE_TO_SUMMARIZE_HETERO_RINGS) {

                for (int indexRingAtom : arrIndexRingAtoms) {
                    arrAtomInSmallRing[indexRingAtom] = true;
                }

            } else {

                for (int indexRingAtom : arrIndexRingAtoms) {
                    arrAtomInLargeRing[indexRingAtom] = true;
                }
            }
        }
    }




    private boolean isSmallRingAtom(int indexAtom){

        return arrAtomInSmallRing[indexAtom];

    }

    /**
     * Takes the atoms in rings <= MAX_RING_SIZE_TO_SUMMARIZE_HETERO_RINGS
     * Also takes end standing atoms which are attached to the ring.
     *
     * @param mol
     * @return
     */
    private List<SubGraphIndices> getSmallRings(StereoMolecule mol){

        List<SubGraphIndices> liFragmentRings = new ArrayList<>();

        RingCollection ringCollection = mol.getRingSet();

        int rings = ringCollection.getSize();

        for (int ringNo = 0; ringNo < rings; ringNo++) {

            int ringSize = ringCollection.getRingSize(ringNo);

            // Exclude enclosing rings
            if(isEnclosingRing(ringCollection, ringNo)){
                continue;
            }

            if(ringSize <= MAX_RING_SIZE_TO_SUMMARIZE_HETERO_RINGS) {

                SubGraphIndices fragment = new SubGraphIndices();

                int [] arrIndexRingAtoms = ringCollection.getRingAtoms(ringNo);

                fragment.addIndex(arrIndexRingAtoms);

                liFragmentRings.add(fragment);
            }
        }

        // Merge bridged rings.
        List<SubGraphIndices> liFragmentRingsMerged = SubGraphIndices.mergeOverlapping(liFragmentRings, 3);

        // Add all end standing exo groups to the ring.
        for (SubGraphIndices subGraphIndices : liFragmentRingsMerged) {

            int [] arrIndexRingAtoms = subGraphIndices.getAtomIndices();

            for (int indexRingAtom : arrIndexRingAtoms) {

                int nConnAtms = mol.getConnAtoms(indexRingAtom);

                if(nConnAtms < 3) {

                    continue;

                }

                for (int i = 0; i < nConnAtms; i++) {

                    int indAtmConn = mol.getConnAtom(indexRingAtom, i);

                    int nConnAtomsChild = mol.getConnAtoms(indAtmConn);

                    if(nConnAtomsChild == 1){

                        subGraphIndices.addIndex(indAtmConn);

                    }
                }
            }
        }

        return liFragmentRingsMerged;
    }


    /**
     * Checks for enclosing ring.
     * An enclosing ring is a ring that contains all atom indices of two or more rings. And the enclosing ring contains
     * no other atom index.
     *
     * @param ringCollection
     * @param ringNo2Check
     * @return
     */
    private boolean isEnclosingRing(RingCollection ringCollection, int ringNo2Check){

        boolean enclosingRing = false;

        int rings = ringCollection.getSize();

        int [] arrIndexRingAtomsRing2Check = ringCollection.getRingAtoms(ringNo2Check);

        for (int indexRingAt : arrIndexRingAtomsRing2Check) {
            arrRingAtom[indexRingAt] = true;
        }

        int sizeRing2Check = arrIndexRingAtomsRing2Check.length;

        LinkedList<LinkedList> liliIndexRing = new LinkedList<>();

        /**
         * Detect all rings that have the same indices as the potentially enclosing ring.
         */
        for (int ringNo = 0; ringNo < rings; ringNo++) {

            if(ringNo==ringNo2Check){
                continue;
            }

            int [] arrIndexRingAtoms = ringCollection.getRingAtoms(ringNo2Check);

            if(arrIndexRingAtoms.length >= sizeRing2Check) {
                continue;
            }

            boolean allIndicesMatch = true;

            for (int indexRingAtom : arrIndexRingAtoms) {

                if(!arrRingAtom[indexRingAtom]){
                    allIndicesMatch = false;
                    break;
                }
            }

            if(allIndicesMatch) {

                LinkedList<Integer> li = new LinkedList<>();

                li.add(ringNo);

                liliIndexRing.add(li);

            }
        }

        while (!liliIndexRing.isEmpty()){

            LinkedList<Integer> liIndexRingParent = liliIndexRing.poll();

            HashSetInt hsIndexAtom = new HashSetInt();

            for (int ringNo : liIndexRingParent) {

                int [] arrIndexRingAtoms = ringCollection.getRingAtoms(ringNo);

                hsIndexAtom.add(arrIndexRingAtoms);

            }

            if(hsIndexAtom.size() == sizeRing2Check) {

                int [] arrIndexRingAtoms = hsIndexAtom.getValues();

                Arrays.sort(arrIndexRingAtoms);

                boolean match = true;

                for (int i = 0; i < arrIndexRingAtoms.length; i++) {

                    if(arrIndexRingAtomsRing2Check[i] != arrIndexRingAtoms[i]){

                        match = false;

                        break;

                    }

                }

                if(match) {

                    enclosingRing = true;

                    break;

                }

            } else if(hsIndexAtom.size() < sizeRing2Check) {

                for (int ringNo = 0; ringNo < rings; ringNo++) {

                    if(ringNo == ringNo2Check) {

                        continue;

                    }

                    if(!liIndexRingParent.contains(ringNo)){

                        int [] arrIndexRingAtoms = ringCollection.getRingAtoms(ringNo);

                        boolean allIndicesMatch = true;

                        for (int indexRingAtom : arrIndexRingAtoms) {
                            if(!arrRingAtom[indexRingAtom]){
                                allIndicesMatch = false;
                                break;
                            }
                        }

                        if(allIndicesMatch) {

                            LinkedList<Integer> liIndexRingChild = new LinkedList<>(liIndexRingParent);

                            liIndexRingChild.add(ringNo);

                            liliIndexRing.add(liIndexRingChild);
                        }
                    }
                }
            }
        }

        // Reset the array
        for (int indexRingAt : arrIndexRingAtomsRing2Check) {
            arrRingAtom[indexRingAt] = false;
        }

        return enclosingRing;
    }

//  10.03.2017 commented out. Not working for multi bridged boron containing rings.
//    private static boolean isEnclosingRing(RingCollection ringCollection, int ringNo){
//
//        boolean enclosingRing = false;
//
//        int rings = ringCollection.getSize();
//
//        int [] arrIndexRingAtoms0 = ringCollection.getRingAtoms(ringNo);
//
//        Arrays.sort(arrIndexRingAtoms0);
//
//        int size0 = arrIndexRingAtoms0.length;
//
//        LinkedList<LinkedList> liliIndexRing = new LinkedList<>();
//
//        for (int i = 0; i < rings; i++) {
//
//            if(i==ringNo){
//
//                continue;
//
//            }
//
//            if(ringCollection.getRingAtoms(i).length < size0) {
//
//                LinkedList<Integer> li = new LinkedList<>();
//
//                li.add(i);
//
//                liliIndexRing.add(li);
//
//            }
//
//        }
//
//        while (!liliIndexRing.isEmpty()){
//
//            LinkedList<Integer> liIndexRingParent = liliIndexRing.poll();
//
//            HashSetInt hsIndexAtom = new HashSetInt();
//
//            for (int indRing : liIndexRingParent) {
//
//                int [] arrIndexRingAtoms = ringCollection.getRingAtoms(indRing);
//
//                hsIndexAtom.add(arrIndexRingAtoms);
//
//            }
//
//            if(hsIndexAtom.size() == size0) {
//
//                int [] arrIndexRingAtoms = hsIndexAtom.getValues();
//
//                Arrays.sort(arrIndexRingAtoms);
//
//                boolean match = true;
//
//                for (int i = 0; i < arrIndexRingAtoms.length; i++) {
//
//                    if(arrIndexRingAtoms0[i] != arrIndexRingAtoms[i]){
//
//                        match = false;
//
//                        break;
//
//                    }
//
//                }
//
//                if(match) {
//
//                    enclosingRing = true;
//
//                    break;
//
//                }
//
//            } else if(hsIndexAtom.size() < size0) {
//
//                for (int i = 0; i < rings; i++) {
//
//                    if(i == ringNo) {
//
//                        continue;
//
//                    }
//
//                    if(!liIndexRingParent.contains(i)){
//
//                        LinkedList<Integer> liIndexRingChild = new LinkedList<>(liIndexRingParent);
//
//                        liIndexRingChild.add(i);
//
//                        liliIndexRing.add(liIndexRingChild);
//
//                    }
//                }
//            }
//        }
//
//        return enclosingRing;
//    }




    /**
     * The largest carbon chain in the ring is taken. The chain may have a minimum length.
     * Rings are skipped that contain atoms that were already used for and end standing aliphatic chain.
     * Rings are skipped that contain atoms already used in a small ring. This prevents the
     *
     * @param mol
     * @return
     */
    private List<SubGraphIndices> getAliphaticGroupsInLargeRings(StereoMolecule mol, HashSetInt hsAtomIndicesInFragmentEndStandingAliphaticGroup) {

        List<SubGraphIndices> liFragmentLargeRings = new ArrayList<>();

        RingCollection ringCollection = mol.getRingSet();

        int rings = ringCollection.getSize();

        boolean [] arrIndAtmUsedInEndStandingAliphaticGroupMap = getMapFromHashSetOfIndices(mol, hsAtomIndicesInFragmentEndStandingAliphaticGroup);

        for (int ringNo = 0; ringNo < rings; ringNo++) {

            int ringSize = ringCollection.getRingSize(ringNo);

            if(ringSize > MAX_RING_SIZE_TO_SUMMARIZE_HETERO_RINGS) {

                int [] arrIndexRingAtoms = ringCollection.getRingAtoms(ringNo);

                boolean atomsInRingAlreadyUsed = false;

                for (int indAtmRing : arrIndexRingAtoms) {

                    if(arrIndAtmUsedInEndStandingAliphaticGroupMap[indAtmRing]){

                        atomsInRingAlreadyUsed = true;

                        break;

                    }
                }

                if(atomsInRingAlreadyUsed) {

                    continue;
                }

                SubGraphIndices fragment = largestAliphaticFragment(mol, arrIndexRingAtoms, MIN_SIZE_ALIPHATIC_CHAIN_IN_RING);

                if(fragment != null) {

                    liFragmentLargeRings.add(fragment);

                }
            }
        }

        //
        // Remove fragments that overlap on more than two atoms with a small ring
        //
        for (int i = liFragmentLargeRings.size()-1; i >= 0; i--) {

            SubGraphIndices fragment = liFragmentLargeRings.get(i);

            for (int ringNo = 0; ringNo < rings; ringNo++) {

                int ringSize = ringCollection.getRingSize(ringNo);

                if (ringSize <= MAX_RING_SIZE_TO_SUMMARIZE_HETERO_RINGS) {

                    boolean [] arrIndAtmMap = getMapFromArrayOfIndices(mol, ringCollection.getRingAtoms(ringNo));

                    int [] arrIndAtmFrag = fragment.getAtomIndices();

                    int nOverlap = 0;

                    for (int indAtmFrag : arrIndAtmFrag) {

                        if(arrIndAtmMap[indAtmFrag]){

                            nOverlap++;
                        }
                    }

                    if(nOverlap > 2) {

                        liFragmentLargeRings.remove(i);

                        break;

                    }
                }
            }
        }

        return liFragmentLargeRings;
    }

    private static boolean [] getMapFromHashSetOfIndices(StereoMolecule mol, HashSetInt hsAtomIndices){

        boolean [] arrIndAtmMap = new boolean[mol.getAtoms()];

        int [] arrIndAtmUsed = hsAtomIndices.getValues();

        for (int indAtmUsed : arrIndAtmUsed) {

            arrIndAtmMap[indAtmUsed] = true;

        }

        return arrIndAtmMap;
    }

    private static boolean [] getMapFromArrayOfIndices(StereoMolecule mol, int [] arrAtomIndices){

        boolean [] arrIndAtmMap = new boolean[mol.getAtoms()];

        for (int indAtmUsed : arrAtomIndices) {

            arrIndAtmMap[indAtmUsed] = true;

        }

        return arrIndAtmMap;
    }


    /**
     * Detects end standing hetero groups. End standing is already a group with =O, -NH2 etc.
     * Also takes atoms into account which are in rings > MAX_RING_SIZE_TO_SUMMARIZE_HETERO_RINGS
     *
     * @param mol
     * @param liEndStandingAtoms
     * @return
     */
    private List<SubGraphIndices> getEndStandingHeteroGroups(StereoMolecule mol, List<Integer> liEndStandingAtoms){

        //
        // Get second layer end standing atoms
        // The parent- or the child atom may be of hetero type.
        //
        List<SubGraphIndices> liFragmentHetero = new ArrayList<>();

        for (int indexEndStandingAtom : liEndStandingAtoms) {

            int atomicNo = mol.getAtomicNo(indexEndStandingAtom);

            int indAtmConn = mol.getConnAtom(indexEndStandingAtom, 0);

            if(isSmallRingAtom(indAtmConn)){

                // If the end standing atom is of this type it will become a own pharmacophore point.
                // Even if it is directly connected to a ring.
                // Additionally it can be included to the atom types of the ring.
                if(atomicNo != PeriodicTable.Nitrogen &&
                        atomicNo != PeriodicTable.Oxygen &&
                        atomicNo != PeriodicTable.Sulfur &&
                        atomicNo != PeriodicTable.Chlorine &&
                        atomicNo != PeriodicTable.Bromine)

                continue;

            }

            if (ExtendedMoleculeFunctions.isHetero(mol, indexEndStandingAtom) || ExtendedMoleculeFunctions.isHetero(mol, indAtmConn)){

                SubGraphIndices fragment = new SubGraphIndices();

                fragment.addIndex(indexEndStandingAtom);
                fragment.addIndex(indAtmConn);

                liFragmentHetero.add(fragment);
            }
        }

        SubGraphIndices.merge(liFragmentHetero);

        //
        // Get third layer
        // Only if the new atom or the connecting atom is a hetero atom.
        //
        for (SubGraphIndices frag : liFragmentHetero) {

            int [] arrAtomIndices = frag.getAtomIndices();

            for (int atmIndex : arrAtomIndices) {

                int nConnAtms = mol.getConnAtoms(atmIndex);

                if(nConnAtms > 1) {

                    boolean heteroParent = ExtendedMoleculeFunctions.isHetero(mol, atmIndex);

                    for (int i = 0; i < nConnAtms; i++) {

                        int indexAtmConn = mol.getConnAtom(atmIndex, i);

                        if(frag.contains(indexAtmConn)) {

                            continue;

                        }

                        boolean heteroChild = ExtendedMoleculeFunctions.isHetero(mol, indexAtmConn);

                        // A hetero atom in a ring is allowed in this layer. This gets the N of an amide in the ring etc.
//                        if (!heteroChild && isSmallRingAtom(indexAtmConn)) {
//                            continue;
//                        }

                        // 07.05.2020 Changed, no hetero atom in ring allowed any more.
                        // Result was a coarse granulation. i.e. for (Ring)N-C-S
                        if (isSmallRingAtom(indexAtmConn)) {
                            continue;
                        }


                        if (heteroParent || heteroChild){
                            frag.addIndex(indexAtmConn);
                        }
                    }
                }
            }
        }

        SubGraphIndices.merge(liFragmentHetero);

        //
        // Get fourth layer
        // Only if the new atom is a hetero atom.
        //

        for (SubGraphIndices frag : liFragmentHetero) {

            int [] arrAtomIndices = frag.getAtomIndices();

            for (int atmIndex : arrAtomIndices) {

                int nConnAtms = mol.getConnAtoms(atmIndex);

                if(nConnAtms > 1) {

                    for (int i = 0; i < nConnAtms; i++) {

                        int indexAtmConn = mol.getConnAtom(atmIndex, i);

                        if (mol.isRingAtom(indexAtmConn)) {

                            continue;

                        } else if(frag.contains(indexAtmConn)) {

                            continue;

                        }

                        if (ExtendedMoleculeFunctions.isHetero(mol, indexAtmConn)){

                            frag.addIndex(indexAtmConn);

                        }
                    }
                }
            }
        }

        SubGraphIndices.merge(liFragmentHetero);

        //
        // Add end standing attached carbon atoms to the group.
        //
        for (SubGraphIndices frag : liFragmentHetero) {

            int [] arrAtomIndices = frag.getAtomIndices();

            for (int atmIndex : arrAtomIndices) {

                int nConnAtms = mol.getConnAtoms(atmIndex);

                if (nConnAtms == 1) {

                    continue;

                }

                for (int i = 0; i < nConnAtms; i++) {

                    int indexAtmConn = mol.getConnAtom(atmIndex, i);

                    if(frag.contains(indexAtmConn)) {

                        continue;

                    }

                    boolean endStandingChild = (mol.getConnAtoms(indexAtmConn) == 1) ? true : false;

                    if(!endStandingChild){

                        continue;

                    }

                    frag.addIndex(indexAtmConn);

                    if(mol.getAtomicNo(indexAtmConn) != 6) {

                        throw new RuntimeException("This should not happen!");

                    }
                }
            }
        }

        return  liFragmentHetero;
    }

    private static class MaximumTopologicalDist {

        int at1;

        int at2;

        int dist;


        @Override
        public String toString() {

            final StringBuilder sb = new StringBuilder();

            sb.append(at1);
            sb.append(" ");
            sb.append(at2);
            sb.append(" ");
            sb.append(dist);

            return sb.toString();
        }
    }

}
