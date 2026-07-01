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

	private static volatile RFKnowledgeBase sKnowledgeBase;

	private TreeMap<Integer, RFKnowledgeBase.RFDetail> mRFDetailMap;
	private TreeMap<Integer, DensityMapsWithDistances> mLigandGeometryMap;
	private TreeMap<Integer, DensityMapsWithDistances> mProteinGeometryMap;

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

	public double getRawUncertainty(int lType, int pType) {
		RFDetail rf = mRFDetailMap.get((pType << 16) | lType);
		return rf == null ? Double.NaN : rf.getUncertainty();
	}

	/**
	 * Calculates the geometry dependent RF-value for the given interaction.
	 * @param ia RF-interaction
	 * @param uncertaintyHolder null or double[1] to receive the uncertainty value
	 * @return RF-value or NaN in case of unknown atom types
	 */
	public double getRFValue(RFInteraction ia, double[] uncertaintyHolder) {
		return getRFValue(ia.getLType(), ia.getL2PGeometryType(), ia.getL2PAngle(), ia.getL2PTorsion(),
						  ia.getPType(), ia.getP2LGeometryType(), ia.getP2LAngle(), ia.getP2LTorsion(),
						  ia.getRelDistance(), uncertaintyHolder);
	}

	private double getRFValue(int lType, int lGeometry, double lAngle, double lTorsion,
							 int pType, int pGeometry, double pAngle, double pTorsion, double relDistance, double[] uncertaintyHolder) {
		double rawRFValue = getRawRFValue(lType, pType);
		if (Double.isNaN(rawRFValue))
			return Double.NaN;
		RFKnowledgeBase.DensityMapsWithDistances lMap = mLigandGeometryMap.get(lGeometry);
		RFKnowledgeBase.DensityMapsWithDistances pMap = mProteinGeometryMap.get(pGeometry);
		if (lMap == null || pMap == null) {
			if (uncertaintyHolder != null)
				uncertaintyHolder[0] = getRawUncertainty(lType, pType);
			return rawRFValue;
		}
		double lFactor = lMap.getDensity(lAngle, lTorsion, relDistance);
		double pFactor = pMap.getDensity(pAngle, pTorsion, relDistance);
		double geomFactor = Math.pow(lFactor * pFactor, GEOMETRY_INFLUENCE);
		if (uncertaintyHolder != null)
			uncertaintyHolder[0] = geomFactor * getRawUncertainty(lType, pType);
		return rawRFValue * geomFactor;
	}

	public String getFullRFDetails(RFInteraction ia) {
		return "rf:"+DoubleFormat.toString(getRFValue(ia, null), 3)
			+ " rawRF:"+DoubleFormat.toString(getRawRFValue(ia.getLType(), ia.getPType()), 3)
			+ "±"+DoubleFormat.toString(getRawUncertainty(ia.getLType(), ia.getPType()), 2)
			+ " relD:"+DoubleFormat.toString(ia.getRelDistance(), 3)
			+"\n     LIG(geo:"+ia.getL2PGeometryName()+" "+mLigandGeometryMap.get(ia.getL2PGeometryType()).getFullInteractionDetails(ia.getL2PAngle(), ia.getL2PTorsion(), ia.getRelDistance())+")"
			+"\n     CAV(geo:"+ia.getP2LGeometryName()+" "+mProteinGeometryMap.get(ia.getP2LGeometryType()).getFullInteractionDetails(ia.getP2LAngle(), ia.getP2LTorsion(), ia.getRelDistance())+")";
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
						RFKnowledgeBase kb = new RFKnowledgeBase();
						kb.readObject(ois);
						ois.close();
						sKnowledgeBase = kb;
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
			}
		}
	}

	public static class RFDetail implements Serializable {
		private static final long serialVersionUID = 0x20260513;

		private float rfL2P,rfP2L,uncertainty;

		public RFDetail(float rfL2P,float rfP2L, float uncertainty) {
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
			stream.writeFloat(rfL2P);
			stream.writeFloat(rfP2L);
			stream.writeFloat(uncertainty);
		}

		private void readObject(ObjectInputStream stream) throws ClassNotFoundException,IOException {
			rfL2P = stream.readFloat();
			rfP2L = stream.readFloat();
			uncertainty = stream.readFloat();
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
		private static final int SHORT = 0;
		private static final int LONG = 1;
		private static final double MIN_DENSITY = 0.01;	// we don't return lower density values to prevent catastrophic influence on scores

		byte[][][] densityGrids;
		double[] distances,meanDensities;

		DensityMapsWithDistances(byte[][][] densityGrids, double[] distances, double[] meanDensities) {
			this.densityGrids = densityGrids;
			this.distances = distances;
			this.meanDensities = meanDensities;
		}

// old lookup based version not adjusted by grid neighbours
//		public double getDensity(double angle, double torsion, double distance) {
//			int angleIndex = Math.min(DENSITY_BINS-1, (int)(angle/DENSITY_BIN_SIZE));
//			int torsionIndex = Math.min(DENSITY_BINS-1, (int)(torsion/DENSITY_BIN_SIZE));
//			double shortFactor = (double)(0xFF & densityGrids[SHORT][angleIndex][torsionIndex]) / meanDensities[SHORT];
//			double longFactor = (double)(0xFF & densityGrids[LONG][angleIndex][torsionIndex]) / meanDensities[LONG];
//			double distancePosition = (distance < distances[SHORT]) ? 0.0 : (distance > distances[LONG]) ? 1.0
//					: (distance - distances[SHORT]) / (distances[LONG] - distances[SHORT]);
//			return (1.0 - distancePosition) * shortFactor + distancePosition * longFactor;
//		}

		/**
		 * Determines an interaction frequency factor for the given angle/torsion combination considering also
		 * the relative interaction distance. The factor is determined by a lookup of the four closest angle/torsion
		 * pairs in the density map and then calculating a neighbor adjusted mean densities. This is done for both,
		 * the long and the short interaction density maps weighted accordingly.
		 * @param angle
		 * @param torsion
		 * @param distance
		 * @return neighbour adjusted density factor considering relative interaction distance
		 */
		public double getDensity(double angle, double torsion, double distance) {
			int angleI1 = (int)((angle-0.5*DENSITY_BIN_SIZE)/DENSITY_BIN_SIZE);
			int angleI2 = Math.min(DENSITY_BINS-1, angleI1+1);
			if (angleI1 == -1)
				angleI1 = 0;
			double angleF2 = ((angle / DENSITY_BIN_SIZE) + 0.5) % 1.0;
			double angleF1 = 1.0 - angleF2;

			int torsionI1 = (int)((torsion-0.5*DENSITY_BIN_SIZE)/DENSITY_BIN_SIZE);
			int torsionI2 = Math.min(DENSITY_BINS-1, torsionI1+1);
			if (torsionI1 == -1)
				torsionI1 = 0;
			double torsionF2 = ((torsion / DENSITY_BIN_SIZE) + 0.5) % 1.0;
			double torsionF1 = 1.0 - torsionF2;

			double shortDensity = angleF1 * torsionF1 * (double)(0xFF & densityGrids[SHORT][angleI1][torsionI1])
								+ angleF1 * torsionF2 * (double)(0xFF & densityGrids[SHORT][angleI1][torsionI2])
								+ angleF2 * torsionF1 * (double)(0xFF & densityGrids[SHORT][angleI2][torsionI1])
								+ angleF2 * torsionF2 * (double)(0xFF & densityGrids[SHORT][angleI2][torsionI2]);

			double longDensity =  angleF1 * torsionF1 * (double)(0xFF & densityGrids[LONG][angleI1][torsionI1])
								+ angleF1 * torsionF2 * (double)(0xFF & densityGrids[LONG][angleI1][torsionI2])
								+ angleF2 * torsionF1 * (double)(0xFF & densityGrids[LONG][angleI2][torsionI1])
								+ angleF2 * torsionF2 * (double)(0xFF & densityGrids[LONG][angleI2][torsionI2]);

			double distanceF2 = (distance < distances[SHORT]) ? 0.0 : (distance > distances[LONG]) ? 1.0
							  : (distance - distances[SHORT]) / (distances[LONG] - distances[SHORT]);
			double distanceF1 = 1.0 - distanceF2;

			return Math.max(MIN_DENSITY, distanceF1 * shortDensity / meanDensities[SHORT]
									   + distanceF2 * longDensity / meanDensities[LONG]);
		}

		public String getFullInteractionDetails(double angle, double torsion, double distance) {
			int angleI1 = (int)((angle-0.5*DENSITY_BIN_SIZE)/DENSITY_BIN_SIZE);
			int angleI2 = Math.min(DENSITY_BINS-1, angleI1+1);
			if (angleI1 == -1)
				angleI1 = 0;
			double angleF2 = ((angle / DENSITY_BIN_SIZE) + 0.5) % 1.0;
			double angleF1 = 1.0 - angleF2;

			int torsionI1 = (int)((torsion-0.5*DENSITY_BIN_SIZE)/DENSITY_BIN_SIZE);
			int torsionI2 = Math.min(DENSITY_BINS-1, torsionI1+1);
			if (torsionI1 == -1)
				torsionI1 = 0;
			double torsionF2 = ((torsion / DENSITY_BIN_SIZE) + 0.5) % 1.0;
			double torsionF1 = 1.0 - torsionF2;

			double shortDensity = angleF1 * torsionF1 * (double)(0xFF & densityGrids[SHORT][angleI1][torsionI1])
					+ angleF1 * torsionF2 * (double)(0xFF & densityGrids[SHORT][angleI1][torsionI2])
					+ angleF2 * torsionF1 * (double)(0xFF & densityGrids[SHORT][angleI2][torsionI1])
					+ angleF2 * torsionF2 * (double)(0xFF & densityGrids[SHORT][angleI2][torsionI2]);

			double longDensity =  angleF1 * torsionF1 * (double)(0xFF & densityGrids[LONG][angleI1][torsionI1])
					+ angleF1 * torsionF2 * (double)(0xFF & densityGrids[LONG][angleI1][torsionI2])
					+ angleF2 * torsionF1 * (double)(0xFF & densityGrids[LONG][angleI2][torsionI1])
					+ angleF2 * torsionF2 * (double)(0xFF & densityGrids[LONG][angleI2][torsionI2]);

			double distanceF2 = (distance < distances[SHORT]) ? 0.0 : (distance > distances[LONG]) ? 1.0
					: (distance - distances[SHORT]) / (distances[LONG] - distances[SHORT]);
			double distanceF1 = 1.0 - distanceF2;

			double f = Math.max(MIN_DENSITY, distanceF1 * shortDensity / meanDensities[SHORT]
										   + distanceF2 * longDensity / meanDensities[LONG]);

			return "f:"+DoubleFormat.toString(f, 3)
				+" ang:"+Math.round(180/Math.PI*angle)
				+" tor:"+Math.round(180/Math.PI*torsion)
				+" dis:"+DoubleFormat.toString(distance,3)
				+" pos:"+DoubleFormat.toString(distanceF2,3)
				+" denS:"+DoubleFormat.toString(shortDensity,3)
				+" denL:"+DoubleFormat.toString(longDensity,3);
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
