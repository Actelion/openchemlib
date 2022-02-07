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

import java.io.Serializable;

public class ReactionSearchSpecification implements Serializable {
    static final long serialVersionUID = 0x20200823;

	private static final int TYPE_MASK				= 0x0000FF;
	public static final int TYPE_NO_REACTION 		= 0x000000; // no structure search (other criteria only)
	public static final int TYPE_SUBREACTION        = 0x000001;
	public static final int TYPE_SIMILARITY			= 0x000002;
	public static final int TYPE_RETRON 			= 0x000003;
	public static final int TYPE_EXACT_STRICT		= 0x000004;
	public static final int TYPE_EXACT_NO_STEREO	= 0x000005;

	private int mSearchType;
	private String[] mQuery;
	private long[][] mReactionDescriptor,mReactantDescriptor,mProductDescriptor,mRetronDescriptor;
	private float mReactionCenterSimilarity,mPeripherySimilarity;

	/**
	 * Creates a complete specification for a reaction similarity search
	 * with one or more query reactions.
	 * @param searchType one of TYPE_...
	 * @param query list of encoded query reactions
	 * @param reactionDescriptor null or query reaction descriptors
	 * @param reactantDescriptor null or reactant FFP512 descriptors
	 * @param productDescriptor null or product FFP512 descriptors
	 * @param reactionCenterSimilarity in case of TYPE_SIMILARITY
	 * @param peripherySimilarity in case of TYPE_SIMILARITY
	 */
	public ReactionSearchSpecification(int searchType, String[] query,
	                                   long[][] reactionDescriptor, long[][] reactantDescriptor, long[][] productDescriptor,
	                                   float reactionCenterSimilarity, float peripherySimilarity) {
		mSearchType = searchType;
		mQuery = query;
		mReactionDescriptor = reactionDescriptor;
		mReactantDescriptor = reactantDescriptor;
		mProductDescriptor = productDescriptor;
		mReactionCenterSimilarity = reactionCenterSimilarity;
		mPeripherySimilarity = peripherySimilarity;
		}

	/**
	 * Creates a complete specification for a reaction substructure search
	 * with one or more query generic reactions.
	 * @param searchType one of TYPE_...
	 * @param query list of encoded query reactions
	 * @param reactionDescriptor null or query reaction descriptors
	 * @param reactantDescriptor null or reactant FFP512 descriptors
	 * @param productDescriptor null or product FFP512 descriptors
	 */
	public ReactionSearchSpecification(int searchType, String[] query,
	                                   long[][] reactionDescriptor, long[][] reactantDescriptor, long[][] productDescriptor) {
		mSearchType = searchType;
		mQuery = query;
		mReactionDescriptor = reactionDescriptor;
		mReactantDescriptor = reactantDescriptor;
		mProductDescriptor = productDescriptor;
	}

	/**
	 * Creates a complete specification for a retron search
	 * with one or more query retron substructures.
	 * @param query list of encoded retron structures
	 * @param retronDescriptor null or product FFP512 descriptors
	 */
	public ReactionSearchSpecification(String[] query, long[][] retronDescriptor) {
		mSearchType = TYPE_RETRON;
		mQuery = query;
		mRetronDescriptor = retronDescriptor;
		}

	/**
	 * Returns the search type as integer including mode flags.
	 * In case of TYPE_SIMILARITY use getDescriptorShortName() and
	 * getSimilarityThreshold() for a full search specification.
	 * @return one of TYPE_... and possibly MODE_LARGEST_FRAGMENT_ONLY
	 *
	public int getSearchType() {
		return mSearchType;
		}	*/

	public int getReactionCount() {
		return (mQuery == null) ? 0 : mQuery.length;
		}

	/**
	 * Returns the (or one of the) query structures encodes as idcode.
	 * @param index
	 * @return
	 */
	public String getEncodedQuery(int index) {
		return mQuery[index];
		}

	/**
	 * @param index
	 * @return
	 */
	public long[] getReactionDescriptor(int index) {
		return (mReactionDescriptor == null) ? null : mReactionDescriptor[index];
		}

	/**
	 * @param index
	 * @return
	 */
	public long[] getReactantDescriptor(int index) {
		return (mReactantDescriptor == null) ? null : mReactantDescriptor[index];
		}

	/**
	 * @param index
	 * @return
	 */
	public long[] getProductDescriptor(int index) {
		return (mProductDescriptor == null) ? null : mProductDescriptor[index];
		}

	/**
	 * @param index
	 * @return
	 */
	public long[] getRetronDescriptor(int index) {
		return (mReactionDescriptor == null) ? null : mReactionDescriptor[index];
		}

	public boolean isSimilaritySearch() {
		return (mSearchType & TYPE_MASK) == TYPE_SIMILARITY;
		}

	/**
	 * @return whether this search does not include a structure search component is uses exclusively external criteria
	 */
	public boolean isNoReactionSearch() {
		return (mSearchType & TYPE_MASK) == TYPE_NO_REACTION;
		}

	public boolean isSubreactionSearch() {
		return (mSearchType & TYPE_MASK) == TYPE_SUBREACTION;
		}

	public boolean isRetronSearch() {
		return (mSearchType & TYPE_MASK) == TYPE_RETRON;
	}

	/**
	 * An exact search is a comparison of idcodes of
	 * standardized molecules with full stereo features.
	 * @return
	 */
	public boolean isExactSearch() {
		return (mSearchType & TYPE_MASK) == TYPE_EXACT_STRICT;
		}

	/**
	 * A noStereo search is a hash code comparison from
	 * encoding stereo depleted structures.
	 * @return
	 */
	public boolean isNoStereoSearch() {
		return (mSearchType & TYPE_MASK) == TYPE_EXACT_NO_STEREO;
		}

	public void removeDescriptors() {
		mReactionDescriptor = null;
		mReactantDescriptor = null;
		mProductDescriptor = null;
		mRetronDescriptor = null;
		}

	public float getReactionCenterSimilarity() {
		return mReactionCenterSimilarity;
		}

	public float getPeripherySimilarity() {
		return mPeripherySimilarity;
		}

	/**
	 * Checks, whether this specification is correctly defining a search.
	 * If something is missing or inconsistent, an error message
	 * describing the problem is returned.
	 * @return null or error message
	 */
	public String validate() {
		int count = getReactionCount();
		if (count == 0)
			return "No query reactions defined.";

		for (int i=0; i<count; i++) {
			if (getEncodedQuery(i) == null || getEncodedQuery(i).length() == 0)   // TODO better validity checking
				return "Empty reaction among query reactions.";
			}

		if (isSimilaritySearch()) {
			if (getReactionCenterSimilarity() < 0.5 || getReactionCenterSimilarity() > 1.0)
				return "Reaction center similarity threshold out of allowed range.";
			if (getPeripherySimilarity() < 0.0 || getPeripherySimilarity() > 1.0)
				return "Periphery similarity threshold out of allowed range.";
			}

		return null;
		}

	@Override
	public String toString() {
		int type = mSearchType & TYPE_MASK;
		String typeString = ((type == TYPE_SUBREACTION) 	? "subreaction"
						   : (type == TYPE_SIMILARITY)		? "similarity("+mReactionCenterSimilarity+"/"+mPeripherySimilarity+")"
						   : (type == TYPE_RETRON)      	? "retron"
						   : (type == TYPE_EXACT_STRICT)	? "exact"
						   : (type == TYPE_EXACT_NO_STEREO)	? "noStereo" : "undefined");

		return "type:"+typeString
			 + (mQuery==null?" reaction:null": mQuery.length==1?" reaction:"+(mQuery[0]==null?"null": mQuery[0]):" reactionCount:"+ mQuery.length)
			 + (mReactionDescriptor==null?"":" reactionDescriptorCount:"+mReactionDescriptor.length)
			 + (mReactantDescriptor==null?"":" reactantDescriptorCount:"+mReactantDescriptor.length)
			 + (mProductDescriptor==null?"":" productDescriptorCount:"+mProductDescriptor.length)
			 + (mRetronDescriptor==null?"":" retronDescriptorCount:"+mRetronDescriptor.length);
		}
	}
