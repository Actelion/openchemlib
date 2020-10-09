package com.actelion.research.chem.descriptor.flexophore.generator;

import com.actelion.research.chem.interactionstatistics.InteractionDistanceStatistics;

import java.util.HashMap;
import java.util.List;

public class AtomTypeMap {


    private static AtomTypeMap INSTANCE;

    private HashMap<Integer, Integer> hmAtomType_Index;

    private int [] arrAtomType;

    private AtomTypeMap() {
        init();
    }

    private void init() {

        InteractionDistanceStatistics interactionDistanceStatistics = InteractionDistanceStatistics.getInstance();
        List<Integer> liAtomTypes = interactionDistanceStatistics.getAtomTypes();

        hmAtomType_Index = new HashMap<>(liAtomTypes.size());

        arrAtomType = new int[liAtomTypes.size()];

        for (int i = 0; i < liAtomTypes.size(); i++) {

            int atType = liAtomTypes.get(i);

            hmAtomType_Index.put(atType, i);

            arrAtomType[i]=atType;

        }
    }

    public int getIndex(int atomType){
        return hmAtomType_Index.get(atomType);
    }

    public int getAtomType(int index){
        return arrAtomType[index];
    }


    public static AtomTypeMap getInstance() {

        if(INSTANCE==null){
            synchronized(AtomTypeMap.class) {
                INSTANCE = new AtomTypeMap();
            }
        }

        return INSTANCE;
    }

    public static void main(String[] args) {

        InteractionDistanceStatistics interactionDistanceStatistics = InteractionDistanceStatistics.getInstance();
        List<Integer> liAtomTypes = interactionDistanceStatistics.getAtomTypes();

        AtomTypeMap atomTypeMap = AtomTypeMap.getInstance();

        for (int i = 0; i < liAtomTypes.size(); i++) {

            int atType = liAtomTypes.get(i);

            int index = atomTypeMap.getIndex(atType);
            int atType2 = atomTypeMap.getAtomType(index);

            if(atType!=atType2){
                throw new RuntimeException("Error in algorithm!");
            }
        }

        System.out.println("Done for " + liAtomTypes.size() + " atom types");



    }

}
