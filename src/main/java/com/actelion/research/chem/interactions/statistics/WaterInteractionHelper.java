package com.actelion.research.chem.interactions.statistics;

public class WaterInteractionHelper {
	private static final double WATER_OXYGEN_DISTANCE = 2.73;
	private static final double WATER_NITROGEN_DISTANCE = 2.87;
	private static final double WATER_ANY_HETERO_DISTANCE = 2.73;
	private static final double WATER_OXYGEN_DISTANCE_TOLERANCE = 2 * 0.08;	// mean distance plus 2.0 * stddev
	private static final double WATER_NITROGEN_DISTANCE_TOLERANCE = 2 * 0.11;
	private static final double WATER_ANY_HETERO_DISTANCE_TOLERANCE = 2 * 0.09;
	private static final double WATER_NEIGHBOR_ANGLE = 111.0/180.0*Math.PI;
	private static final double WATER_NEIGHBOR_ANGLE_TOLERANCE = 15.0/180.0*Math.PI;

	public static boolean qualifiesAsWaterNeighbour(int atomicNo, double distance) {
		if (atomicNo == 7) {
			if (Math.abs(distance-WATER_NITROGEN_DISTANCE) < WATER_NITROGEN_DISTANCE_TOLERANCE)
				return true;
		}
		else if (atomicNo == 8) {
			if (Math.abs(distance-WATER_OXYGEN_DISTANCE) < WATER_OXYGEN_DISTANCE_TOLERANCE)
				return true;
		}
		else {
			if (Math.abs(distance-WATER_ANY_HETERO_DISTANCE) < WATER_ANY_HETERO_DISTANCE_TOLERANCE)
				return true;
		}
		return false;
	}

	public static boolean qualifiesAsWaterAngle(double angle) {
		return (Math.abs(angle) - WATER_NEIGHBOR_ANGLE) < WATER_NEIGHBOR_ANGLE_TOLERANCE;
	}
}
