package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule3D;

public class MolDistHistVizHelper {

    public static void centerNodes(MolDistHistViz mdhv){
        Coordinates barycenter = getBarycenter(mdhv);
        for (int i = 0; i < mdhv.getNumPPNodes(); i++) {
            PPNodeViz node = mdhv.getNode(i);
            node.getCoordinates().x -= barycenter.x;
            node.getCoordinates().y -= barycenter.y;
            node.getCoordinates().z -= barycenter.z;
        }

        Molecule3D m = mdhv.getMolecule();
        Coordinates [] arrCoordinates = m.getCoordinates();
        for (int i = 0; i < arrCoordinates.length; i++) {
            Coordinates c = arrCoordinates[i];
            c.x -= barycenter.x;
            c.y -= barycenter.y;
            c.z -= barycenter.z;
            m.setCoordinates(i, c);
        }
    }

    public static Coordinates getBarycenter(MolDistHistViz mdhv){
        Coordinates [] arrCoordinates = new Coordinates[mdhv.getNumPPNodes()];
        for (int i = 0; i < mdhv.getNumPPNodes(); i++) {
            PPNodeViz node = mdhv.getNode(i);
            arrCoordinates[i] = node.getCoordinates();
        }
        return Coordinates.createBarycenter(arrCoordinates);
    }


}
