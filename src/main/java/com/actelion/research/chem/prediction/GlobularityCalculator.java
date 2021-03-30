package com.actelion.research.chem.prediction;

import com.actelion.research.calc.SingularValueDecomposition;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.conf.ConformerSetGenerator;
import org.openmolecules.chem.conf.gen.ConformerGenerator;

import java.util.Arrays;

public class GlobularityCalculator {
	private static final boolean MINIMIZE = false;

	/**
	 * Generates conformers from a 2D- or 3D-molecule and then calculates the globularity
	 * (flat=0, round=1) from the conformer's atom coordinates.
	 * @param mol molecule used to generate conformers
	 * @param maxConformerCount maximum number of conformers to generate
	 * @return globularity (flat=0, round=1)
	 */
	public static double assessGlobularity(StereoMolecule mol, int maxConformerCount) {
		ConformerSet cs = new ConformerSetGenerator(maxConformerCount, ConformerGenerator.STRATEGY_LIKELY_RANDOM, MINIMIZE,0).generateConformerSet(mol);
		return assessGlobularity(cs);
		}

	/**
	 * Generates conformers from a 2D- or 3D-molecule and then calculates the globularity
	 * (flat=0, round=1) from the conformer's atom coordinates.
	 * @param cs conformer set from which to calculate and average the globularity
	 * @return globularity (flat=0, round=1)
	 */
	public static double assessGlobularity(ConformerSet cs) {
		if (cs == null || cs.size() == 0)
			return Double.NaN;

		double globularity = 0.0;
		for (Conformer conf:cs)
			globularity += assessGlobularityFromSVD(conf);

		return globularity / cs.size();
		}

	/**
	 * Calculates the globularity of a conformer (flat=0, round=1) from 3D coordinates
	 * @param conf
	 * @return
	 */
	public static double assessGlobularityFromSVD(Conformer conf) {
		int atomCount = conf.getSize();
		Coordinates cog = new Coordinates();	// center of gravity
		for (int i=0; i<atomCount; i++)
			cog.add(conf.getCoordinates(i));
		cog.scale(1.0/atomCount);

		double[][] squareMatrix = new double[3][3];
		for (int i=0; i<atomCount; i++) {
			conf.getCoordinates(i).sub(cog);
			squareMatrix[0][0] += conf.getX(i) * conf.getX(i);
			squareMatrix[0][1] += conf.getX(i) * conf.getY(i);
			squareMatrix[0][2] += conf.getX(i) * conf.getZ(i);
			squareMatrix[1][0] += conf.getY(i) * conf.getX(i);
			squareMatrix[1][1] += conf.getY(i) * conf.getY(i);
			squareMatrix[1][2] += conf.getY(i) * conf.getZ(i);
			squareMatrix[2][0] += conf.getZ(i) * conf.getX(i);
			squareMatrix[2][1] += conf.getZ(i) * conf.getY(i);
			squareMatrix[2][2] += conf.getZ(i) * conf.getZ(i);
			}

		SingularValueDecomposition svd = new SingularValueDecomposition(squareMatrix, null, null);

		double[][] u = svd.getU();
		double[] temp = new double[3];
		for (int i=0; i<atomCount; i++) {
			Arrays.fill(temp, 0);
			Coordinates c = conf.getCoordinates(i);
			for (int j = 0; j<3; j++) {
				temp[j] += c.x * u[0][j];
				temp[j] += c.y * u[1][j];
				temp[j] += c.z * u[2][j];
				}
			c.set(temp[0], temp[1], temp[2]);
			}
		double xmin = Double.MAX_VALUE;
		double xmax = Double.MIN_VALUE;
		double ymin = Double.MAX_VALUE;
		double ymax = Double.MIN_VALUE;
		double zmin = Double.MAX_VALUE;
		double zmax = Double.MIN_VALUE;
		for (int i=0; i<atomCount; i++) {
			Coordinates c = conf.getCoordinates(i);
			if (xmin > c.x)
				xmin = c.x;
			if (xmax < c.x)
				xmax = c.x;
			if (ymin > c.y)
				ymin = c.y;
			if (ymax < c.y)
				ymax = c.y;
			if (zmin > c.z)
				zmin = c.z;
			if (zmax < c.z)
				zmax = c.z;
			}
		double dx = xmax - xmin;
		double dy = ymax - ymin;
		double dz = zmax - zmin;

		return dx == 0.0 ? Double.NaN : dz / dx;
		}
	}
