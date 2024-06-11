/*
 * Copyright 2013-2020 Thomas Sander, openmolecules.org
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.gen;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.IDCodeParserWithoutCoordinateInvention;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.io.CompoundFileParser;
import com.actelion.research.gui.FileHelper;
import com.actelion.research.util.DoubleFormat;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.TreeSet;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.zip.ZipInputStream;

/**
 * This class implements a thread-save, concurrent cache of rigid fragments' 3D-atom-coordinates.
 * It is accessed by the RigidFragmentProvider instances, which serve RigidFragments to
 * ConformerGenerator instances when constructing 3D-coordinates for molecules by assembling
 * them from 3-dimensional rigid fragments and torsion tables.
 * Typically, ConformerGenerators start with an empty cache that fills over time or with
 * a default cache, which is prefilled with many common fragments from organic and medicinal
 * chemistry as well as with common building block fragments.<br>
 * The default cache is balanced in memory footprint and number of fragments it contains.
 * For special purposes you may consider creating an own custom cache file using the createCacheFiles() method.
 **/
public class RigidFragmentCache extends ConcurrentHashMap<String, RigidFragmentCache.CacheEntry> implements Serializable {
	private static final int DEFAULT_MAX_ENTRY_COUNT = 500000;
	private static final String DEFAULT_CACHE_FILE = "/resources/defaultRigidFragments.zip";
	private static RigidFragmentCache sInstance;
	private int mHitCount,mGetCount,mNonCachableCount,mMaxEntryCount;
	private boolean mDefaultCacheLoaded;
	private TreeSet<String> mSetOfLoadedCacheFiles;

	public static RigidFragmentCache getDefaultInstance() {
		if (sInstance != null)
			return sInstance;

		synchronized (RigidFragmentCache.class) {
			if (sInstance == null)
				sInstance = new RigidFragmentCache();
			return sInstance;
		}
	}

	public static RigidFragmentCache createInstance(String cacheFileName) {
		RigidFragmentCache cache = new RigidFragmentCache();
		if (cacheFileName != null)
			cache.loadCache(cacheFileName);
		return cache;
	}

	private RigidFragmentCache() {
		mMaxEntryCount = DEFAULT_MAX_ENTRY_COUNT;
	}

	@Override
	public void clear() {
		super.clear();
		mDefaultCacheLoaded = false;
	}

	@Override
	public RigidFragmentCache.CacheEntry get(Object key) {
		RigidFragmentCache.CacheEntry entry = super.get(key);
		mGetCount++;
		if (entry != null) {
			entry.incrementHitCount();
			mHitCount++;
			}
		return entry;
		}

	public double getHitQuote() {
		return (double)mHitCount/(double)mGetCount;
		}

	public int getHitCount() {
		return mHitCount;
		}

	public int getRequestCount() {
		return mGetCount;
		}

	public int getNonCachableCount() {
		return mNonCachableCount;
		}

	public void increaseNonCachableCount() {
		mNonCachableCount++;
		}

	public void resetAllCounters() {
		mNonCachableCount = 0;
		mHitCount = 0;
		mGetCount = 0;
		}

	public void setMaxEntryCount(int count) {
		mMaxEntryCount = count;
		}

	public boolean canAddEntry() {
		return size() < mMaxEntryCount;
	}

	@Override
	public RigidFragmentCache.CacheEntry put(String key, RigidFragmentCache.CacheEntry cacheEntry) {
		return (size() < mMaxEntryCount) ? super.put(key, cacheEntry) : null;
	}

	/**
	 * Writes for every distinct fragment: one idcode, multiple encoded coordinate sets, multiple conformer likelihoods
	 * @param cacheFileName
	 * @param minHits number of hits for a cache entry to be included in the cache file
	 */
	public boolean serializeCache(String cacheFileName, int minHits) {
		try{
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(cacheFileName),"UTF-8"));
			for (String key : keySet()) {
				RigidFragmentCache.CacheEntry cacheEntry = super.get(key);  // we need super to not increment hit counter

				if (cacheEntry.hitCount >= minHits) {
					bw.write(key);
					bw.newLine();

					bw.write(Integer.toString(cacheEntry.coordinates.length));
					bw.newLine();

					StereoMolecule mol = new IDCodeParserWithoutCoordinateInvention().getCompactMolecule(key);
					Canonizer canonizer = new Canonizer(mol, Canonizer.COORDS_ARE_3D);
					for (Coordinates[] coords:cacheEntry.coordinates) {
						bw.write(canonizer.getEncodedCoordinates(true, coords));
						bw.newLine();
						canonizer.invalidateCoordinates();
					}

					for (double likelihood:cacheEntry.likelihood) {
						bw.write(DoubleFormat.toString(likelihood));
						bw.newLine();
					}
				}
			}
			bw.close();
			return true;
		} catch (IOException ex) {
			ex.printStackTrace();
		}
		return false;
	}

	/**
	 * This loads the default cache file
	 */
	public synchronized void loadDefaultCache() {
		if (!mDefaultCacheLoaded) {
			try {
				InputStream is = RigidFragmentCache.class.getResourceAsStream(DEFAULT_CACHE_FILE);
				if (is != null) {
					ZipInputStream zipStream = new ZipInputStream(is);
					zipStream.getNextEntry();
					BufferedReader reader = new BufferedReader(new InputStreamReader(zipStream, StandardCharsets.UTF_8));
					loadCache(reader);
					reader.close();
					mDefaultCacheLoaded = true;
					}
				}
			catch (Exception e) {
				e.printStackTrace();
				}
			}
		}

	/**
	 * Loads pre-calculated rigid fragment coordinates from a cache file, which is either a text file
	 * created by the createCacheFiles() method, or a zip archive of the text file.
	 * This method can be called multiple times to add conformer data from multiple sources.
	 * If the method is called with a cacheFileNam, which was loaded before, then it is not loaded a second time.
	 * @param cacheFileName text file or zipped text file with extension .zip
	 */
	public void loadCache(String cacheFileName) {
		if (mSetOfLoadedCacheFiles == null)
			mSetOfLoadedCacheFiles = new TreeSet();

		if (!mSetOfLoadedCacheFiles.contains(cacheFileName)) {
			try {
				BufferedReader reader;
				if (cacheFileName.endsWith(".zip")) {
					ZipInputStream zipStream = new ZipInputStream(new FileInputStream(cacheFileName));
					zipStream.getNextEntry();
					reader = new BufferedReader(new InputStreamReader(zipStream, StandardCharsets.UTF_8));
				}
				else {
					reader = new BufferedReader(new FileReader(cacheFileName));
				}
				loadCache(reader);
				reader.close();
			}
			catch (Exception e) {
				e.printStackTrace();
			}
			mSetOfLoadedCacheFiles.add(cacheFileName);
		}
	}

	private void loadCache(BufferedReader br) throws Exception {
		String idcode;
		while ((idcode = br.readLine()) != null) {
			IDCodeParserWithoutCoordinateInvention parser = new IDCodeParserWithoutCoordinateInvention();
			StereoMolecule mol = parser.getCompactMolecule(idcode);

			int count = Integer.parseInt(br.readLine());

			Coordinates[][] coords = new Coordinates[count][mol.getAllAtoms()];
			for (int i=0; i<count; i++) {
				for (int j=0; j<coords[i].length; j++)
					coords[i][j] = new Coordinates();
				parser.parseCoordinates(br.readLine().getBytes(StandardCharsets.UTF_8), 0, mol, coords[i]);
			}

			double[] likelihood = new double[count];
			for (int i=0; i<count; i++)
				likelihood[i] = Double.parseDouble(br.readLine());

			put(idcode, new CacheEntry(coords, likelihood));
		}
	}

	/**
	 * This is a helper method to generate a custom cache and optionally a set of cache files from one or
	 * more compound files. You may use this function to create a custom fragment cache if the default cache
	 * file used by the ConformerGenerator in not adequate for your purpose. The default file covers many
	 * common fragments in organic and medicinal chemistry and common building block fragments.
	 * However, it is limited in size. You may consider using a custom cache file in these cases:<br>
	 * - To achieve a maximum of speed on the expense of memory, e.g. for a cloud based service that
	 * generates conformers on request.<br>
	 * - If you process molecules with limited diversity, e.g. combinatorial libraries as the Enamine REAL space.
	 * Then you may use a complete cache covering every existing fragment for maximum speed.<br>
	 * - If you store conformer sets as fragment references and torsion tables. Then your fragment cache
	 * needs a complete cache covering every existing fragment.<br>
	 * This method processes all input files, locates and all rigid fragments, produces one or more
	 * distinct conformers from the fragments and creates a new cache from them. Optionally, the
	 * fragment conformers can be energy minimized using the MMFF94s+ forcefield. Then multiple cache
	 * cache export files are written: with all cache entries, with entries used at least 2,3,5, and 10 times.
	 * The numbers 1,2,3,5,10 and .txt extention will be appended to the given cache file name.
	 * @param inputFileNames array of one or more input file paths (may be mixture of sdf and dwar)
	 * @param outputDirectory path to output directory ('_cache_n_.txt' will be added)
	 * @param threadCount if 1 then a single threaded approach is used; if 0 then all existing cores are used
	 * @param optimizeFragments whether to energy minimize fragments using MMFF94s+
	 * @param maxCompoundsPerFile if an input file contains more compounds than this, then the rest are skipped
	 * @param rfp null or custom RigidFragmentProvider if fragments shall be minimized with a different method
	 * @return created cache or null, if an input file could not be found
	 */
	public static RigidFragmentCache createCache(String[] inputFileNames, String outputDirectory, int threadCount,
                                boolean optimizeFragments, int maxCompoundsPerFile, RigidFragmentProvider rfp) {
		boolean notFound = false;
		for (String ifn:inputFileNames)
			if (!FileHelper.fileExists(new File(ifn), 1000)) {
				System.out.println("File not found: '"+ifn+"'");
				notFound = true;
			}
		if (notFound)
			return null;

		RigidFragmentCache cache = createInstance(null);
		if (rfp != null)
			rfp.setCache(cache);

		for (String ifn:inputFileNames) {
			long millis = (threadCount != 1) ?
					addFragmentsToCacheSMP(cache, optimizeFragments, ifn, maxCompoundsPerFile, rfp, threadCount)
				  : addFragmentsToCache(cache, optimizeFragments, ifn, maxCompoundsPerFile, rfp);

			System.out.println("File '"+ifn+"' processed in "+millis+" milliseconds.");
		}

		if (inputFileNames != null) {
			System.out.print("Writing cache files... ");
			String cacheFileName = outputDirectory.concat("/cache_");
			boolean success = cache.serializeCache(cacheFileName + "1.txt", 0)   // we have one hit less than usages
					&& cache.serializeCache(cacheFileName + "2.txt", 1)
					&& cache.serializeCache(cacheFileName + "3.txt", 2)
					&& cache.serializeCache(cacheFileName + "5.txt", 4)
					&& cache.serializeCache(cacheFileName + "10.txt", 9);
			System.out.println(success ? "done" : "failure !!!");
			}

		return cache;
	}

	private static long addFragmentsToCache(RigidFragmentCache cache, boolean optimizeFragments,
	                                        String inputFile, int maxCompounds, RigidFragmentProvider rfp) {
		long start_millis = System.currentTimeMillis();
		int compoundNo = 0;

		System.out.println("Processing '"+inputFile+"'... ('.' = 50 molecules)");

		CompoundFileParser parser = CompoundFileParser.createParser(inputFile);

		ConformerGenerator cg = (rfp == null) ?
				new ConformerGenerator(123L, cache, optimizeFragments)
				: new ConformerGenerator(123L, rfp);

		while (parser.next() && compoundNo < maxCompounds) {
			if (compoundNo % 50 == 49)
				System.out.print(".");
			if (compoundNo % 5000 == 4999) {
				System.out.println(" hit-rate:" + DoubleFormat.toString(cache.getHitQuote(), 5, false)
						+ " millis:" + (System.currentTimeMillis() - start_millis)
						+ " cacheSize:" + cache.size());
				cache.resetAllCounters();
			}

			cg.initialize(parser.getMolecule(), false);

			compoundNo++;
		}
		System.out.println();

		return System.currentTimeMillis() - start_millis;
	}

	private static long addFragmentsToCacheSMP(RigidFragmentCache cache, boolean optimizeFragments,
	                                        String inputFile, int maxCompounds, RigidFragmentProvider rfp, int threadCount) {
		long start_millis = System.currentTimeMillis();
		int compoundNo = 0;

		System.out.println("Processing '" + inputFile + "'... ('.' = 50 molecules)");

		CompoundFileParser parser = CompoundFileParser.createParser(inputFile);

		if (threadCount == 0)
			threadCount = Runtime.getRuntime().availableProcessors();

		ArrayBlockingQueue<StereoMolecule> queue = new ArrayBlockingQueue<>(2*threadCount);
		Thread[] t = new Thread[threadCount];
		for (int i = 0; i<threadCount; i++) {
			t[i] = new Thread(() -> consumeMoleculesToCacheFragments(queue, cache, optimizeFragments, rfp));
			t[i].setPriority(Thread.MIN_PRIORITY);
			t[i].start();
		}

		while (parser.next() && compoundNo<maxCompounds) {
			if (compoundNo % 50 == 49)
				System.out.print(".");
			if (compoundNo % 5000 == 4999) {
				System.out.println(" hit-rate:" + DoubleFormat.toString(cache.getHitQuote(), 5, false)
						+ " millis:" + (System.currentTimeMillis() - start_millis)
						+ " cacheSize:" + cache.size());
				cache.resetAllCounters();
			}

			try {
				StereoMolecule mol = parser.getMolecule();
				if (mol.getAllAtoms() != 0)
					queue.put(mol);
			}
			catch (InterruptedException ie) {}

			compoundNo++;
		}

		for (int i=0; i<threadCount; i++)
			t[i].interrupt();
		for (int i=0; i<threadCount; i++)
			try { t[i].join(); } catch (InterruptedException e) {}

		System.out.println();

		return System.currentTimeMillis() - start_millis;
	}

	private static void consumeMoleculesToCacheFragments(ArrayBlockingQueue<StereoMolecule> queue, RigidFragmentCache cache,
	                                                     boolean optimizeFragments, RigidFragmentProvider rfp) {
		ConformerGenerator cg = (rfp == null) ?
				new ConformerGenerator(123L, cache, optimizeFragments)
				: new ConformerGenerator(123L, rfp);

		try {
			while (true)
				cg.initialize(queue.take(), false);
		}
		catch (InterruptedException ie) {}  // spawning thread interrupts this after last molecule
	}

	/**
	 * Writes a TAB delimited text file that can be opened for debug or other purposes by DataWarrior containing
	 * idcode, idcoords,  multiple conformer likelihoods
	 * @param cacheFileName
	 */
	public boolean writeTabDelimitedTable(String cacheFileName) {
		try {
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(cacheFileName),"UTF-8"));
			bw.write("Fragment No\tConformer No\tConformer Count\tidcode\tidcoords\tLikelihood");
			bw.newLine();

			int fragmentCount = 0;
			for (String key : keySet()) {
				fragmentCount++;

				RigidFragmentCache.CacheEntry cacheEntry = super.get(key);  // we need super to not increment hit counter

				StereoMolecule mol = new IDCodeParserWithoutCoordinateInvention().getCompactMolecule(key);
				Canonizer canonizer = new Canonizer(mol, Canonizer.COORDS_ARE_3D);

				for (int i=0; i<cacheEntry.coordinates.length; i++) {
					bw.write(Integer.toString(fragmentCount));
					bw.write("\t");
					bw.write(Integer.toString(i+1));
					bw.write("\t");
					bw.write(Integer.toString(cacheEntry.coordinates.length));
					bw.write("\t");
					bw.write(key);
					bw.write("\t");
					bw.write(canonizer.getEncodedCoordinates(true, cacheEntry.coordinates[i]));
					bw.write("\t");
					bw.write(DoubleFormat.toString(cacheEntry.likelihood[i]));
					bw.newLine();

					canonizer.invalidateCoordinates();
				}
			}
			bw.close();
			return true;
		} catch (IOException ex) {
			ex.printStackTrace();
		}
		return false;
	}

	public static class CacheEntry implements Comparable<CacheEntry> {
		Coordinates[][] coordinates;
		double[] likelihood;
		int hitCount;

		public CacheEntry(Coordinates[][] coordinates, double[] likelihoods) {
			this.coordinates = coordinates;
			this.likelihood = likelihoods;
		}

		public void incrementHitCount() {
			hitCount++;
		}

		@Override
		public int compareTo(CacheEntry o) {
			if (hitCount != o.hitCount)
				return hitCount < o.hitCount ? -1 : 1;
			return 0;
		}
	}
}
