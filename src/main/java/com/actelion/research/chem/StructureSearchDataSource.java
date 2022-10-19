/*
 * Copyright (c) 1997 - 2016
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the the copyright holder nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Thomas Sander
 */

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
	 * This is the number of individual structures within that row that a query will be checked against.
	 * Typically, this will be one. However, if a database contains empty structures, mixtures,
	 * or reactants/products of a reaction, then the row may contain multiple structures to be searched.
	 * @param row
	 * @return typically 1, but may be 0...n
	 */
	public int getStructureCount(int row);

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
	 * @param i structure index for that particular row
	 * @param largestFragmentOnly
	 * @return null if this row has no descriptor
	 */
	public Object getDescriptor(int column, int row, int i, boolean largestFragmentOnly);

	/**
	 * Returns the idcode of the structure or largest fragment.
	 * If the code is not available, null is returned.
	 * This method is not supposed to calculate idcodes on-the-fly.
	 * @param row
	 * @param i structure index for that particular row
	 * @param largestFragmentOnly
	 * @return idcode of this row or null, if no structure is associated to this row
	 */
	public byte[] getIDCode(int row, int i, boolean largestFragmentOnly);

	/**
	 * Returns a hash code representing the structure or its largest fragment
	 * without any stereo information. If this code is not available, SEARCH_TYPE_NOT_SUPPORTED
	 * is returned. This method is not supposed to calculate the code on-the-fly.
	 * @param row
	 * @param i structure index for that particular row
	 * @param largestFragmentOnly
	 * @return
	 */
	public long getNoStereoCode(int row, int i, boolean largestFragmentOnly);

	/**
	 * Returns a hash code representing the generic tautomer of the structure or its
	 * largest fragment. 
	 * If this code is not available, SEARCH_TYPE_NOT_SUPPORTED is returned.
	 * This method is not supposed to calculate the code on-the-fly.
	 * @param row
	 * @param i structure index for that particular row
	 * @param largestFragmentOnly
	 * @return
	 */
	public long getTautomerCode(int row, int i, boolean largestFragmentOnly);

	/**
	 * Returns a hash code representing the generic tautomer of the structure or its
	 * largest fragment with the stereo information removed before creating the generic tautomer.
	 * If this code is not available, SEARCH_TYPE_NOT_SUPPORTED is returned.
	 * This method is not supposed to calculate the code on-the-fly.
	 * @param row
	 * @param i structure index for that particular row
	 * @param largestFragmentOnly
	 * @return
	 */
	public long getNoStereoTautomerCode(int row, int i, boolean largestFragmentOnly);

	/**
	 * Returns a hash code representing the structure or its largest fragment
	 * without stereo information and with all unsaturated bonds converted to single bonds.
	 * SEARCH_TYPE_NOT_SUPPORTED is returned.
	 * This method is not supposed to calculate the code on-the-fly.
	 * @param row
	 * @param i structure index for that particular row
	 * @param largestFragmentOnly
	 * @return
	 */
	public long getBackboneCode(int row, int i, boolean largestFragmentOnly);
	}
