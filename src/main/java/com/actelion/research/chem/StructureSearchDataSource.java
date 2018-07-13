package com.actelion.research.chem;

public interface StructureSearchDataSource {
	public static final long SEARCH_TYPE_NOT_SUPPORTED = -1L;

	/**
	 * Returns true if idcodes and/or hash codes for this particular
	 * search specification are directly available by the source, i.e. without
	 * calculating idcodes, descriptors or hash codes on-the-fly.
	 * In case of search types other than substructure search the implementation needs
	 * to check specification.isLargestFragmentOnly(). In case of similarity
	 * searches the implementation needs to consider specification.getDescriptorShortName().
	 * @param specification
	 * @return
	 */
	public boolean isSupportedSearchType(StructureSearchSpecification specification);

	/**
	 * This is the total number of records, whether they qualify for a structure search or not.
	 * @return
	 */
	public int getRowCount();

	/**
	 * Returns a number to be used in getDescriptor() to address the
	 * the descriptor with the given short name.
	 * If a descriptor of the desired kind is not available, then return -1
	 * as indication that it needs to be generated on-the-fly.
	 * @param descriptorShortName
	 * @return descriptor column or -1, if the descriptor is not available
	 */
	public int getDescriptorColumn(String descriptorShortName);

	/**
	 * @param column
	 * @param row
	 * @param largestFragmentOnly
	 * @return null if this row has no descriptor
	 */
	public Object getDescriptor(int column, int row, boolean largestFragmentOnly);

	/**
	 * Returns the idcode of the structure or largest fragment.
	 * If the code is not available, null is returned.
	 * This method is not supposed to calculate idcodes on-the-fly.
	 * @param row
	 * @param largestFragmentOnly
	 * @return idcode of this row or null, if no structure is associated to this row
	 */
	public byte[] getIDCode(int row, boolean largestFragmentOnly);

	/**
	 * Returns a hash code representing the structure or its largest fragment
	 * without any stereo information. If this code is not available, SEARCH_TYPE_NOT_SUPPORTED
	 * is returned. This method is not supposed to calculate the code on-the-fly.
	 * @param row
	 * @param largestFragmentOnly
	 * @return
	 */
	public long getNoStereoCode(int row, boolean largestFragmentOnly);

	/**
	 * Returns a hash code representing the generic tautomer of the structure or its
	 * largest fragment. 
	 * If this code is not available, SEARCH_TYPE_NOT_SUPPORTED is returned.
	 * This method is not supposed to calculate the code on-the-fly.
	 * @param row
	 * @param largestFragmentOnly
	 * @return
	 */
	public long getTautomerCode(int row, boolean largestFragmentOnly);

	/**
	 * Returns a hash code representing the generic tautomer of the structure or its
	 * largest fragment with the stereo information removed before creating the generic tautomer.
	 * If this code is not available, SEARCH_TYPE_NOT_SUPPORTED is returned.
	 * This method is not supposed to calculate the code on-the-fly.
	 * @param row
	 * @param largestFragmentOnly
	 * @return
	 */
	public long getNoStereoTautomerCode(int row, boolean largestFragmentOnly);

	/**
	 * Returns a hash code representing the structure or its largest fragment
	 * without stereo information and with all unsaturated bonds converted to single bonds.
	 * SEARCH_TYPE_NOT_SUPPORTED is returned.
	 * This method is not supposed to calculate the code on-the-fly.
	 * @param row
	 * @param largestFragmentOnly
	 * @return
	 */
	public long getBackboneCode(int row, boolean largestFragmentOnly);
	}
