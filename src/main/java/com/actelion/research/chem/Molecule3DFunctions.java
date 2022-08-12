package com.actelion.research.chem;

import java.io.*;
import java.util.*;

public class Molecule3DFunctions {

    public static final Molecule3D removeAllAtomsWithoutNeighbours(Molecule3D mol) {
        Molecule3D molecule3D = new Molecule3D(mol);
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

    public static final String toStringSerialized(Molecule3D mol) throws IOException {

        ByteArrayOutputStream byteArrayOutputStream = new ByteArrayOutputStream();
        ObjectOutputStream oos = new ObjectOutputStream(byteArrayOutputStream);
        oos.writeObject(mol);
        oos.close();

        return Base64.getEncoder().encodeToString(byteArrayOutputStream.toByteArray());
    }
    public static final Molecule3D readSerialized(String serMol) throws IOException, ClassNotFoundException {
        byte [] arrByte = Base64.getDecoder().decode(serMol);
        ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(arrByte));
        Molecule3D molecule3D  = (Molecule3D)is.readObject();
        is.close();
        return molecule3D;
    }

}
