package com.actelion.research.chem.descriptor;

import com.actelion.research.chem.SSSearcherWithIndex;
import com.actelion.research.chem.StereoMolecule;

public class DescriptorHandlerLongFFP512 extends AbstractDescriptorHandlerLongFP<StereoMolecule> {
	public static final String VERSION = SSSearcherWithIndex.cIndexVersion;
	private static DescriptorHandlerLongFFP512 sDefaultInstance;
	private static final int sLongCount = (SSSearcherWithIndex.getNoOfKeys() + 63) / 64;

	public static DescriptorHandlerLongFFP512 getDefaultInstance() {
		synchronized(DescriptorHandlerLongFFP512.class) {
			if (sDefaultInstance == null)
				sDefaultInstance = new DescriptorHandlerLongFFP512();
		}
		return sDefaultInstance;
	}

	public DescriptorInfo getInfo() {
		return DescriptorConstants.DESCRIPTOR_FFP512;
	}

	public String getVersion() {
		return VERSION;
	}

	public long[] decode(String s) {
		long[] descriptor = (s != null && s.length() == 128) ?
				SSSearcherWithIndex.getLongIndexFromHexString(s) : super.decode(s);
		return (descriptor != null && descriptor.length == sLongCount) ? descriptor : null;
	}

	public long[] decode(byte[] bytes) {
		long[] descriptor = (bytes != null && bytes.length == 128) ?
				SSSearcherWithIndex.getLongIndexFromHexString(bytes) : super.decode(bytes);
		return (descriptor != null && descriptor.length == sLongCount) ? descriptor : null;
	}

	public long[] createDescriptor(StereoMolecule mol) {
		long[] descriptor = new SSSearcherWithIndex().createLongIndex(mol);
		return (descriptor == null) ? FAILED_OBJECT : descriptor;
	}

	public DescriptorHandler<long[], StereoMolecule> getThreadSafeCopy() {
		return this;
	}
}
