package com.actelion.research.util.graph.complete;

/**
 * 
 * 
 * IFactory
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Sep 27, 2012 MvK: Start implementation
 */
public interface IFactory<S extends AMemorizedObject> {

	public S createObject();
	
}
