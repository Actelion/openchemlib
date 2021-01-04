package com.actelion.research.gui.dock;

import java.awt.BorderLayout;
import java.awt.Component;
import java.util.ArrayList;

import javax.swing.*;

public class TreeRoot extends TreeContainer {
    private TreeElement mChildElement;

    /**
     * Constructor to create a root element to which the first leaf should be connected.
     * @param parent
     */
    public TreeRoot(JComponent rootComponent, TreeElement child) {
        mComponent = rootComponent;
        mComponent.add(child.getComponent(), BorderLayout.CENTER);
        mChildElement = child;
        mChildElement.setParent(this);
        }

    public TreeElement getChild() {
        return mChildElement;
        }

    public void setParent(TreeContainer parent) {
        throw new IllegalArgumentException("no parent to root");
        }

    public void removeWithLeaf(TreeLeaf leaf) {
        leaf.setParent(null);
        }

    public void replaceChildElement(TreeElement oldElement, TreeElement newElement) {
        mComponent.remove(oldElement.getComponent());
        mComponent.add(newElement.getComponent(), BorderLayout.CENTER);
        newElement.setParent(this);
        mChildElement = newElement;
        }

    public void updateChildElement(Component oldContent, TreeElement childElement) {
        mComponent.remove(oldContent);
        mComponent.add(childElement.getComponent());
        }

    protected void clearStateInfo() {
        mChildElement.clearStateInfo();
        }

    public ArrayList<String> createStateInfo() {
        clearStateInfo();
        ArrayList<String> stateInfo = new ArrayList<String>();
        if (mChildElement instanceof TreeLeaf)
            ((TreeLeaf)mChildElement).addStateInfo(stateInfo, "root");
        if (mChildElement instanceof TreeFork)
            ((TreeFork)mChildElement).addStateInfo(stateInfo, "root");

        return stateInfo;
        }

    public void printStatus() {
        if (mComponent == null)
            System.out.println("Root childStatus: none");
        else {
            String status = (mComponent.getComponent(0) == mChildElement.getComponent()) ? "OK" : "failed";
            System.out.println("Root childStatus:"+status);
            mChildElement.printStatus();
            }
        }
    }

