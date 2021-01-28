package com.actelion.research.util.datamodel;

import java.util.ArrayList;
import java.util.List;

/**
 * ModelXYColTags
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Aug 13, 2015 MvK Start implementation
 */
public class ModelXYColTags extends ModelXYClassLabel {

	
	public List<String> liTagX;
	
	public List<String> liTagY;


	public ModelXYColTags(ModelXYColTags modelXYColTags) {
		super(modelXYColTags);
		this.liTagX = new ArrayList<>(modelXYColTags.liTagX);
		this.liTagY = new ArrayList<>(modelXYColTags.liTagY);
	}

	/**
	 * 
	 */
	public ModelXYColTags() {
	}


	@Override
	public ModelXYColTags getDeepClone() {
		return new ModelXYColTags(this);
	}


}
