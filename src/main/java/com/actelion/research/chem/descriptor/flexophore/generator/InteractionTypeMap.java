package com.actelion.research.chem.descriptor.flexophore.generator;

import com.actelion.research.chem.interactionstatistics.InteractionDistanceStatistics;
import com.actelion.research.chem.interactionstatistics.InteractionSimilarityTable;
import com.actelion.research.util.datamodel.IntegerDouble;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class InteractionTypeMap {

    private double threshSimilarity;

    private HashMap<Integer, List<IntegerDouble>> hmAtomType_ListSimilars;

    private InteractionSimilarityTable interactionSimilarityTable;

    private List<Integer> liAtomTypes;

    public InteractionTypeMap(double threshSimilarity) {

        this.threshSimilarity = threshSimilarity;

        interactionSimilarityTable = InteractionSimilarityTable.getInstance();

        InteractionDistanceStatistics interactionDistanceStatistics = InteractionDistanceStatistics.getInstance();
        liAtomTypes = interactionDistanceStatistics.getAtomTypes();

        hmAtomType_ListSimilars = new HashMap<>();
    }

    /**
     *
     * @param interactionTypeQuery
     * @return descending sorted list by similarity.
     */
    public List<IntegerDouble> getSimilars(int interactionTypeQuery){

        List<IntegerDouble> li = hmAtomType_ListSimilars.get(interactionTypeQuery);

        if(li==null){

            synchronized (InteractionTypeMap.class){

                List<IntegerDouble> liInteractionType_Similarity = new ArrayList<>(liAtomTypes.size());

                for (int interactionTypeBase : liAtomTypes) {
                    double similarity = 1.0 - interactionSimilarityTable.getDistance(interactionTypeQuery, interactionTypeBase);
                    if(similarity >= threshSimilarity) {
                        liInteractionType_Similarity.add(new IntegerDouble(interactionTypeBase, similarity));
                    }
                }

                Collections.sort(liInteractionType_Similarity, IntegerDouble.getComparatorDouble());
                Collections.reverse(liInteractionType_Similarity);

                hmAtomType_ListSimilars.put(interactionTypeQuery, liInteractionType_Similarity);

                li = liInteractionType_Similarity;
            }
        }

        return li;
    }

}
