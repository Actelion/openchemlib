package com.actelion.research.chem;

import com.actelion.research.chem.descriptor.DescriptorHelper;

import java.io.Serializable;

public class StructureSearchSpecification implements Serializable {
    static final long serialVersionUID = 0x20120402;

	private static final int TYPE_MASK				= 0x0000FF;
	public static final int TYPE_NO_STRUCTURE		= 0x000000; // no structure search (other criteria only)
	public static final int TYPE_SUBSTRUCTURE		= 0x000001;
	public static final int TYPE_SIMILARITY			= 0x000002;
	public static final int TYPE_EXACT_STRICT		= 0x000003;
	public static final int TYPE_EXACT_NO_STEREO	= 0x000004;
	public static final int TYPE_TAUTOMER			= 0x000005;
	public static final int TYPE_TAUTOMER_NO_STEREO	= 0x000006;
	public static final int TYPE_BACKBONE_NO_STEREO	= 0x000007;

	public static final int MODE_LARGEST_FRAGMENT_ONLY		= 0x000100;

	private int mSearchType;
	private byte[][] mIDCode;
	private Object[] mDescriptor;
	private String mDescriptorShortName;
	private float mSimilarityThreshold;

	/**
	 * Creates a complete specification for a substructure or similarity search
	 * with one or more query fragments or molecules.
	 * @param searchType one of TYPE_... + optionally MODE_...
	 * @param idcode list of query fragments/molecules
	 * @param descriptor null or list descriptors (in case of substructure search this should be the long FFP512)
	 * @param descriptorShortName in case of TYPE_SIMILARITY
	 * @param similarityThreshold in case of TYPE_SIMILARITY
	 */
	public StructureSearchSpecification(int searchType, byte[][] idcode, Object[] descriptor, String descriptorShortName, float similarityThreshold) {
		mSearchType = searchType;
		mIDCode = idcode;
		mDescriptor = descriptor;
		mDescriptorShortName = descriptorShortName;
		mSimilarityThreshold = similarityThreshold;
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

	public int getStructureCount() {
		return (mIDCode == null) ? 0 : mIDCode.length;
		}

	/**
	 * Returns the (or one of the) query structures encodes as idcode.
	 * @param index
	 * @return
	 */
	public byte[] getIDCode(int index) {
		return mIDCode[index];
		}

	/**
	 * Depending on the mode isLargestFragmentOnly() the descriptor
	 * is expected to be represent the whole structure or the largest fragment.
	 * @param index
	 * @return
	 */
	public Object getDescriptor(int index) {
		return (mDescriptor == null) ? null : mDescriptor[index];
		}

	public boolean isSimilaritySearch() {
		return (mSearchType & TYPE_MASK) == TYPE_SIMILARITY;
		}

	/**
	 * @return whether this search does not include a structure search component is uses exclusively external criteria
	 */
	public boolean isNoStructureSearch() {
		return (mSearchType & TYPE_MASK) == TYPE_NO_STRUCTURE;
		}

	public boolean isSubstructureSearch() {
		return (mSearchType & TYPE_MASK) == TYPE_SUBSTRUCTURE;
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

	/**
	 * A tautomer search is a hash code comparison from
	 * encoding generic tautomer structures.
	 * @return
	 */
	public boolean isTautomerSearch() {
		return (mSearchType & TYPE_MASK) == TYPE_TAUTOMER;
		}

	/**
	 * A no-stereo-tautomer search is a hash code comparison from
	 * encoding stereo depleted generic tautomer structures.
	 * @return
	 */
	public boolean isNoStereoTautomerSearch() {
		return (mSearchType & TYPE_MASK) == TYPE_TAUTOMER_NO_STEREO;
		}

	/**
	 * A backbone search is a hash code comparison from
	 * encoding stereo and unsaturation depleted structures.
	 * @return
	 */
	public boolean isBackboneSearch() {
		return (mSearchType & TYPE_MASK) == TYPE_BACKBONE_NO_STEREO;
		}

	/**
	 * This setting is relevant for all search types except the substructure search.
	 * @return
	 */
	public boolean isLargestFragmentOnly() {
		return (mSearchType & MODE_LARGEST_FRAGMENT_ONLY) != 0;
		}

	public void removeDescriptors() {
		mDescriptor = null;
		}

	public void setLargestFragmentOnly(boolean b) {
		mSearchType &= ~MODE_LARGEST_FRAGMENT_ONLY;
		if (b)
			mSearchType |= MODE_LARGEST_FRAGMENT_ONLY;
		}

	public String getDescriptorShortName() {
		return mDescriptorShortName;
		}

	public float getSimilarityThreshold() {
		return mSimilarityThreshold;
		}

	/**
	 * Checks, whether this specification is correctly defining a search.
	 * If something is missing or inconsistent, an error message
	 * describing the problem is returned.
	 * @return null or error message
	 */
	public String validate() {
		int count = getStructureCount();
		if (count == 0)
			return "No query structures defined.";

		for (int i=0; i<count; i++)
			if (new IDCodeParser(false).getAtomCount(getIDCode(i), 0) == 0)
				return "Empty structure among query structures.";

		if (isSimilaritySearch()) {
			if (!DescriptorHelper.isDescriptorShortName(getDescriptorShortName()))
				return "No descriptor type defined.";
			if (getSimilarityThreshold() < 0.5 || getSimilarityThreshold() > 1.0)
				return "Similarity threshold out of allowed range.";
			}

		return null;
		}

	@Override
	public String toString() {
		int type = mSearchType & TYPE_MASK;
		String typeString = ((type == TYPE_SUBSTRUCTURE)	? "substructure"
						   : (type == TYPE_SIMILARITY)		? "similarity/"+mDescriptorShortName+"("+mSimilarityThreshold+")"
						   : (type == TYPE_EXACT_STRICT)	? "exact"
						   : (type == TYPE_EXACT_NO_STEREO)	? "noStereo"
						   : (type == TYPE_TAUTOMER)		? "tautomer"
						   : (type == TYPE_TAUTOMER_NO_STEREO)	? "no-stereo tautomer"
						   : (type == TYPE_BACKBONE_NO_STEREO)	? "backbone" : "undefined")
						  + (((mSearchType & MODE_LARGEST_FRAGMENT_ONLY) != 0) ? "/largestFragmentOnly":"");

		return "type:"+typeString
			 + (mIDCode==null?" idcodes:null":mIDCode.length==1?" idcode:"+(mIDCode[0]==null?"null":new String(mIDCode[0])):" idcodeCount:"+mIDCode.length)
			 + (mDescriptor==null?" descriptors:null":" descriptorCount:"+mDescriptor.length);
		}
	}
