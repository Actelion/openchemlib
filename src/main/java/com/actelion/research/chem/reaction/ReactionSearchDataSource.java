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
	 * FFP of all reactants merged into one molecule
	 * @param row
	 * @return null if this row has no reactant descriptor
	 */
	public long[] getMergedReactantDescriptor(int row);

	/**
	 * FFP of all products merged into one molecule
	 * @param row
	 * @return null if this row has no product descriptor
	 */
	public long[] getMergedProductDescriptor(int row);

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