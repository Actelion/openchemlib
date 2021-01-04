package com.actelion.research.gui.dock;

import java.awt.Component;

public abstract class TreeContainer extends TreeElement {
    public abstract void removeWithLeaf(TreeLeaf leaf);
    public abstract void replaceChildElement(TreeElement oldElement, TreeElement newElement);
    public abstract void updateChildElement(Component oldContent, TreeElement childElement);
    }

