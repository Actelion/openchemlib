package org.openmolecules.chem.interaction.rf;

import com.actelion.research.chem.interactions.statistics.InteractionDistanceStatistics;

import java.io.*;
import java.net.URL;
import java.util.TreeMap;

public class RFKnowledgeBase extends TreeMap<Integer, RFKnowledgeBase.RFDetail> {
	private static final long serialVersionUID = 0x20260513;

	public static final String FILE_NAME = "rfInteractionDB.bin";

	private static RFKnowledgeBase sKnowledgeBase;

	public static RFKnowledgeBase getInstance() {
		ensureInitialization();
		return sKnowledgeBase;
	}

	public static int getLigandType(int combinedType) {
		return combinedType & 0xFFFF;
	}

	public static int getProteinType(int combinedType) {
		return combinedType >> 16;
	}

	public static RFKnowledgeBase load(File saved_file) {
		try (FileInputStream load = new FileInputStream(saved_file);
			 ObjectInputStream in = new ObjectInputStream(load)) {
			return (RFKnowledgeBase) in.readObject();
		} catch (Exception ioe) {
			return null;
		}
	}

	public static boolean save(String dirPath) {
		try {
			ObjectOutputStream os = new ObjectOutputStream(new FileOutputStream(dirPath + File.separator + FILE_NAME));
			os.writeObject(sKnowledgeBase);
			os.close();
			return true;
		}
		catch (Exception e) {
			e.printStackTrace();
			return false;
		}
	}

	public static void addRFValue(int type, RFDetail detail) {
		if (sKnowledgeBase == null)
			sKnowledgeBase = new RFKnowledgeBase();

		sKnowledgeBase.put(type, detail);
	}

	private static void ensureInitialization() {
		if (sKnowledgeBase == null) {
			synchronized (RFKnowledgeBase.class) {
				if (sKnowledgeBase == null) {
					try {
						URL url =  InteractionDistanceStatistics.class.getResource("/resources/"+FILE_NAME);
						if (url == null)
							throw new RuntimeException("Could not find file '"+FILE_NAME+"' in the classpath or resources.");
						ObjectInputStream is = new ObjectInputStream(url.openStream());
						sKnowledgeBase = (RFKnowledgeBase)is.readObject();
						is.close();
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
			}
		}
	}

	public static double getRFValue(int lType, int pType) {
		ensureInitialization();
		RFDetail rf = sKnowledgeBase.get((pType << 16) | lType);
		return rf == null ? Double.NaN : rf.getRF();
	}

	public static double getRFValue(int lType, int pType, double lAngle, double lDist, double pAngle, double pDist) {
		ensureInitialization();
		RFDetail rf = sKnowledgeBase.get((pType << 16) | lType);
		// TODO calculate from angle and dist parameters
		return rf == null ? Double.NaN : rf.getRF();
	}

	public static class RFDetail implements Serializable {
		private static final long serialVersionUID = 0x20260513;

		private final double rfL2P,rfP2L,uncertainty;

		public RFDetail(double rfL2P,double rfP2L, double uncertainty) {
			this.rfL2P = rfL2P;
			this.rfP2L = rfP2L;
			this.uncertainty = uncertainty;
		}

		public double getRFL2P() {
			return rfL2P;
		}

		public double getRFP2L() {
			return rfP2L;
		}

		public double getRF() {
			return Math.sqrt(rfL2P * rfP2L);
		}

		public double getUncertainty() {
			return uncertainty;
		}
	}
}
