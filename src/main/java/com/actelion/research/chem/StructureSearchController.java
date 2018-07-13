package com.actelion.research.chem;

public interface StructureSearchController {
	/**
	 * If the structure search specification is part of a broader criteria list,
	 * then this method determines, whether all other, i.e. non-structure search
	 * related, criteria match. If other non-matching criteria exist, then this
	 * method returns false. If no other search criteria exist, then true is returned.
	 * This method's implementation must be thread safe.
	 * @param row
	 * @return
	 */
	public boolean rowQualifies(int row);
	}
