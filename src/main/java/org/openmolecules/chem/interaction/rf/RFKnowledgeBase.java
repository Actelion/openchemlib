package org.openmolecules.chem.interaction.rf;

import com.actelion.research.util.DoubleFormat;

import java.io.*;
import java.net.URL;
import java.util.Set;
import java.util.TreeMap;

public class RFKnowledgeBase implements Serializable {
	private static final long serialVersionUID = 0x20260611;

	public static final String FILE_NAME = "rfInteractionDB.bin";
	private static final int DENSITY_BINS = 36;
	private static final double DENSITY_BIN_SIZE = Math.PI / DENSITY_BINS;
	private static final double GEOMETRY_INFLUENCE = 0.33;	// larger values increase RF reduction with bad geometries

	private static RFKnowledgeBase sKnowledgeBase;

	private TreeMap<Integer, RFKnowledgeBase.RFDetail> mRFDetailMap;
	private TreeMap<Integer, DensityMapsWithDistances> mLigandGeometryMap;
	private TreeMap<Integer, DensityMapsWithDistances> mProteinGeometryMap;

//	public static void createOldKnowledgeBaseFile() {
//		ensureInitialization();
//		RFKnowledgeBase kb = RFKnowledgeBase.createEmptyInstance();
//		for (int key : sKnowledgeBase.mRFDetailMap.keySet()) {
//			NewRFKnowledgeBase.RFDetail detail = sKnowledgeBase.mRFDetailMap.get(key);
//			kb.addRFValue(key, new RFKnowledgeBase.RFDetail(detail.rfL2P, detail.rfP2L, detail.uncertainty));
//		}
//		kb.save("/home/thomas/dev/local/xtal/pdbinteractions");
//	}

	private RFKnowledgeBase() {
		mRFDetailMap = new TreeMap<>();
		mLigandGeometryMap = new TreeMap<>();
		mProteinGeometryMap = new TreeMap<>();
		mLigandGeometryMap = new TreeMap<>();
		mProteinGeometryMap = new TreeMap<>();
	}

	public static RFKnowledgeBase getInstance() {
		ensureInitialization();
		return sKnowledgeBase;
	}

	public void setLigandGeometryDensities(int geometryType, byte[][][] densities, double[] distances, double[] meanDensities) {
		mLigandGeometryMap.put(geometryType, new DensityMapsWithDistances(densities, distances, meanDensities));
	}

	public void setProteinGeometryDensities(int geometryType, byte[][][] densities, double[] distances, double[] meanDensities) {
		mProteinGeometryMap.put(geometryType, new DensityMapsWithDistances(densities, distances, meanDensities));
	}

	/**
	 * Returns the raw RF-value as a geometric mean from ligand and protein perspective raw RF-values.
	 * These are for the given combination of ligand and protein atoms types the count of experimentally
	 * found close, line-of-sight interaction divided by the number of expected interactions close
	 * line-of-sight interaction considering<br>
	 * - exposed surface area fraction of atoms of the considered type (ligand or protein side)<br>
	 * - all binding sites that contains considered ligand and protein atom types at least once<br>
	 * This method does not consider interaction geometries on ligand not protein side (angles or torsions).
	 * @param lType
	 * @param pType
	 * @return
	 */
	public double getRawRFValue(int lType, int pType) {
		RFDetail rf = mRFDetailMap.get((pType << 16) | lType);
		return rf == null ? Double.NaN : rf.getRF();
	}

	public double getRFValue(RFInteraction ia) {
		return getRFValue(ia.getLType(), ia.getL2PGeometryType(), ia.getL2PAngle(), ia.getL2PTorsion(),
						  ia.getPType(), ia.getP2LGeometryType(), ia.getP2LAngle(), ia.getP2LTorsion(),
						  ia.getRelDistance());
	}

	private double getRFValue(int lType, int lGeometry, double lAngle, double lTorsion,
							 int pType, int pGeometry, double pAngle, double pTorsion, double relDistance) {
		double rawRFValue = getRawRFValue(lType, pType);
		if (Double.isNaN(rawRFValue))
			return Double.NaN;
		double lFactor = mLigandGeometryMap.get(lGeometry).getDensity(lAngle, lTorsion, relDistance);
		double pFactor = mProteinGeometryMap.get(pGeometry).getDensity(pAngle, pTorsion, relDistance);
		return rawRFValue * Math.pow(lFactor * pFactor, GEOMETRY_INFLUENCE);
	}

	public String getFullRFDetails(RFInteraction ia) {
		return "rf:"+DoubleFormat.toString(getRFValue(ia))
			+ " rawRF:"+DoubleFormat.toString(getRawRFValue(ia.getLType(), ia.getPType()))
			+" LIG(geo:"+ia.getL2PGeometryName()+" "+mLigandGeometryMap.get(ia.getL2PGeometryType()).getFullInteractionDetails(ia.getL2PAngle(), ia.getL2PTorsion(), ia.getRelDistance())+")"
			+" CAV(geo:"+ia.getP2LGeometryName()+" "+mProteinGeometryMap.get(ia.getP2LGeometryType()).getFullInteractionDetails(ia.getP2LAngle(), ia.getP2LTorsion(), ia.getRelDistance())+")";
	}

	public static RFKnowledgeBase createEmptyInstance() {
		return new RFKnowledgeBase();
	}

	public static int getLigandType(int combinedType) {
		return combinedType & 0xFFFF;
	}

	public static int getProteinType(int combinedType) {
		return combinedType >> 16;
	}

	public static RFKnowledgeBase load(String filePath) {
		RFKnowledgeBase knowledgeBase = new RFKnowledgeBase();
		try (FileInputStream fis = new FileInputStream(filePath);
			ObjectInputStream ois = new ObjectInputStream(fis)) {
			knowledgeBase.readObject(ois);
		} catch (Exception ioe) {}
		return knowledgeBase;
	}

	public boolean save(String filePath) {
		try {
			ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(filePath));
			writeObject(oos);
			oos.close();
			return true;
		}
		catch (Exception e) {
			e.printStackTrace();
			return false;
		}
	}

	public void addRFValue(int type, RFDetail detail) {
		mRFDetailMap.put(type, detail);
	}

	private static void ensureInitialization() {
		if (sKnowledgeBase == null) {
			synchronized (RFKnowledgeBase.class) {
				if (sKnowledgeBase == null) {
					try {
						URL url =  RFKnowledgeBase.class.getResource("/resources/"+FILE_NAME);
						if (url == null)
							throw new RuntimeException("Could not find file '"+FILE_NAME+"' in the classpath or resources.");
						ObjectInputStream ois = new ObjectInputStream(url.openStream());
						sKnowledgeBase = new RFKnowledgeBase();
						sKnowledgeBase.readObject(ois);
						ois.close();
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
			}
		}
	}

	public static class RFDetail implements Serializable {
		private static final long serialVersionUID = 0x20260513;

		private double rfL2P,rfP2L,uncertainty;

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

		private void writeObject(ObjectOutputStream stream) throws IOException {
			stream.writeDouble(rfL2P);
			stream.writeDouble(rfP2L);
			stream.writeDouble(uncertainty);
		}

		private void readObject(ObjectInputStream stream) throws ClassNotFoundException,IOException {
			rfL2P = stream.readDouble();
			rfP2L = stream.readDouble();
			uncertainty = stream.readDouble();
		}
	}

	public Set<Integer> getInteractionKeySet() {
		return mRFDetailMap.keySet();
	}

	public RFDetail getInteractionDetail(int key) {
		return mRFDetailMap.get(key);
	}

	private void writeObject(ObjectOutputStream stream) throws IOException {
		stream.writeObject(mRFDetailMap);
		stream.writeObject(mLigandGeometryMap);
		stream.writeObject(mProteinGeometryMap);
	}

	private void readObject(ObjectInputStream stream) throws ClassNotFoundException,IOException {
		mRFDetailMap = (TreeMap<Integer, RFKnowledgeBase.RFDetail>)stream.readObject();
		mLigandGeometryMap = (TreeMap<Integer, RFKnowledgeBase.DensityMapsWithDistances>)stream.readObject();
		mProteinGeometryMap = (TreeMap<Integer, RFKnowledgeBase.DensityMapsWithDistances>)stream.readObject();
	}

	/**
	 * Two or more distances with associates smooth density grid that describe
	 */
	public static class DensityMapsWithDistances implements Serializable {
		private static final long serialVersionUID = 0x20260611;

		byte[][][] densityGrids;
		double[] distances,meanDensities;

		DensityMapsWithDistances(byte[][][] densityGrids, double[] distances, double[] meanDensities) {
			this.densityGrids = densityGrids;
			this.distances = distances;
			this.meanDensities = meanDensities;
		}

		public double getDensity(double angle, double torsion, double distance) {
			int angleIndex = Math.min(DENSITY_BINS-1, (int)(angle/DENSITY_BIN_SIZE));
			int torsionIndex = Math.min(DENSITY_BINS-1, (int)(torsion/DENSITY_BIN_SIZE));
			double shortFactor = (double)(0xFF & densityGrids[0][angleIndex][torsionIndex]) / meanDensities[0];
			double longFactor = (double)(0xFF & densityGrids[1][angleIndex][torsionIndex]) / meanDensities[1];
			double distancePosition = (distance < distances[0]) ? 0.0 : (distance > distances[1]) ? 1.0
					: (distance - distances[0]) / (distances[1] - distances[0]);
			return (1.0 - distancePosition) * shortFactor + distancePosition * longFactor;
		}

		public String getFullInteractionDetails(double angle, double torsion, double distance) {
			int angleIndex = Math.min(DENSITY_BINS-1, (int)(angle/DENSITY_BIN_SIZE));
			int torsionIndex = Math.min(DENSITY_BINS-1, (int)(torsion/DENSITY_BIN_SIZE));
			double shortFactor = (0xFF & densityGrids[0][angleIndex][torsionIndex]) / meanDensities[0];
			double longFactor = (0xFF & densityGrids[1][angleIndex][torsionIndex]) / meanDensities[1];
			double distancePosition = (distance < distances[0]) ? 0.0 : (distance > distances[1]) ? 1.0
					: (distance - distances[0]) / (distances[1] - distances[0]);
			return "denF:"+DoubleFormat.toString((1.0 - distancePosition) * shortFactor + distancePosition * longFactor, 3)
				+" ang:"+Math.round(180/Math.PI*angle)
				+" tor:"+Math.round(180/Math.PI*torsion)
				+" dis:"+DoubleFormat.toString(distance,3)
				+" pos:"+DoubleFormat.toString(distancePosition,3)
				+" denS:"+(0xFF & densityGrids[0][angleIndex][torsionIndex])
				+" denL:"+(0xFF & densityGrids[1][angleIndex][torsionIndex]);
		}

		private void writeObject(ObjectOutputStream stream) throws IOException {
			stream.writeObject(densityGrids);
			stream.writeObject(distances);
			stream.writeObject(meanDensities);
		}

		private void readObject(ObjectInputStream stream) throws ClassNotFoundException,IOException {
			densityGrids = (byte[][][])stream.readObject();
			distances = (double[])stream.readObject();
			meanDensities = (double[])stream.readObject();
		}
	}
}
