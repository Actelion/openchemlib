package org.openmolecules.chem.conf.gen;

import java.util.concurrent.ConcurrentHashMap;

public class RigidFragmentCache {
	private int mMaxCacheSize;
	private volatile ConcurrentHashMap mConformerMap;

	public RigidFragmentCache() {
		this(10000);
	}

	public RigidFragmentCache(int maxCacheSize) {

	}
}
