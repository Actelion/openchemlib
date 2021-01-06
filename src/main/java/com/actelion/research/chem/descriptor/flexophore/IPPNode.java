package com.actelion.research.chem.descriptor.flexophore;

public interface IPPNode {
    int getInteractionTypeCount();
    int getInteractionType(int i);
    boolean hasHeteroAtom();
}
