package com.actelion.research.chem;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

public class Molecule3DFunctions {

    public static final Molecule3D removeAllAtomsWithoutNeighbours(Molecule3D ffMol) {
        Molecule3D molecule3D = new Molecule3D(ffMol);
        molecule3D.ensureHelperArrays(Molecule.cHelperRings);

        HashSet<Integer> hsAt2Del = new HashSet<Integer>();
        for (int i = 0; i < molecule3D.getAllAtoms(); i++) {
            if(molecule3D.getConnAtoms(i)==0)
                hsAt2Del.add(i);
        }

        List<Integer> liAt2Del = new ArrayList<Integer>(hsAt2Del);
        Collections.sort(liAt2Del);
        Collections.reverse(liAt2Del);

        for (Integer at : liAt2Del) {
            molecule3D.deleteAtom(at);
        }

        molecule3D.ensureHelperArrays(Molecule.cHelperRings);

        return molecule3D;
    }


}
