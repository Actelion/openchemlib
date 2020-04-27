package com.actelion.research.util.graph.complete;

/**
 * 
 * 
 * AMemorizedObject
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Sep 28, 2012 MvK: Start implementation
 */
public abstract class AMemorizedObject {
	
	private int positionInContainer;
	
	/**
	 * Has to be a deep copy
	 * @param m
	 */
	abstract public void copyIntoThis(AMemorizedObject m);
	
	abstract public void reset();
	
	protected int getPositionInContainer() {
		return positionInContainer;
	}

	protected void setPositionInContainer(int positionInContainer) {
		this.positionInContainer = positionInContainer; 
	}
}
