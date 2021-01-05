package com.actelion.research.gui.dock;

import javax.swing.*;

public abstract class TreeElement {
    protected TreeContainer mParent;
    protected JComponent mComponent;

    public TreeContainer getParent() {
        return mParent;
        }

    public void setParent(TreeContainer parent) {
        mParent = parent;
        }

    public JComponent getComponent() {
        return mComponent;
        }

    public boolean isSelected() {
        return false;
        }

    abstract protected void clearStateInfo();
    abstract protected void printStatus();
    }
