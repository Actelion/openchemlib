package com.actelion.research.chem.reaction;

public interface ReactionSearchDataSource {
	public static final long SEARCH_TYPE_NOT_SUPPORTED = -1L;

	/**
	 * Returns true if idcodes and/or hash codes for this particular
	 * search specification are directly available by the source, i.e. without
	 * calculating idcodes, descriptors or hash codes on-the-fly.
	 * In case of search types other than substructure search the implementation needs
	 * to check specification.isLargestFragmentOnly(). In case of similarity
	 * searches the implementation needs to consider specification.getDescriptorShortName().
	 *
	 * @param specification
	 * @return
	 */
	public boolean isSupportedSearchType(ReactionSearchSpecification specification);

	/**
	 * This is the total number of records, whether they qualify for a structure search or not.
	 *
	 * @return
	 */
	public int getRowCount();

	/**
	 * @param row
	 * @return null if this row has no reaction descriptor
	 */
	public long[] getReactionDescriptor(int row);

	/**
	 * @param row
	 * @return null if this row has no reactant descriptor
	 */
	public long[] getReactantDescriptor(int row);

	/**
	 * @param row
	 * @return null if this row has no product descriptor
	 */
	public long[] getProductDescriptor(int row);

	/**
	 * Returns the reaction idcode.
	 * If the code is not available, null is returned.
	 * This method is not supposed to calculate idcodes on-the-fly.
	 *
	 * @param row
	 * @return idcode of this row or null, if no reaction is associated to this row
	 */
	public byte[] getReactionCode(int row);

	/**
	 * Returns the encoded reaction coordinates.
	 *
	 * @param row
	 * @return coords of this row or null, if no reaction is associated to this row
	 */
	public byte[] getCoordinates(int row);

	/**
	 * Returns the encoded reaction mapping.
	 *
	 * @param row
	 * @return mapping of this row or null, if no reaction is associated to this row
	 */
	public byte[] getMapping(int row);

	/**
	 * Returns a hash code representing the exact reaction.
	 * If this code is not available, SEARCH_TYPE_NOT_SUPPORTED is returned.
	 * This method is not supposed to calculate the code on-the-fly.
	 *
	 * @param row
	 * @return
	 */
	public long getExactHash(int row);

	/**
	 * Returns a hash code representing the reaction without any stereo information.
	 * If this code is not available, SEARCH_TYPE_NOT_SUPPORTED is returned.
	 * This method is not supposed to calculate the code on-the-fly.
	 *
	 * @param row
	 * @return
	 */
	public long getNoStereoHash(int row);
}