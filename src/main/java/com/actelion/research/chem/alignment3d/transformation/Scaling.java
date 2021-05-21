package com.actelion.research.chem.alignment3d.transformation;

import java.util.Base64;
import java.util.Base64.Encoder;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.phesa.EncodeFunctions;

/**
 * @author wahljo1
 *
 */
public class Scaling implements Transformation {
	
	public double getScalingFactor() {
		return scalingFactor;
	}

	public void setScalingFactor(double scalingFactor) {
		this.scalingFactor = scalingFactor;
	}

	private double scalingFactor;
	
	public Scaling(double scalingFactor) {
		this.scalingFactor = scalingFactor;
	}

	@Override
	public void apply(Coordinates coords) {
		coords.scale(scalingFactor);
		
	}

	@Override
	public void apply(double[] coords) {
		coords[0]*=scalingFactor;
		coords[1]*=scalingFactor;
		coords[2]*=scalingFactor;
		
	}

	@Override
	public String encode() {
		Encoder encoder = Base64.getEncoder();
		StringBuilder sb = new StringBuilder();
		sb.append(Transformation.SCALING);
		sb.append(Transformation.DELIMITER);
		sb.append(encoder.encodeToString(EncodeFunctions.doubleToByteArray(scalingFactor)));
		return sb.toString();
	}
	
	public static Scaling decode(String s) {
		double scaling = EncodeFunctions.byteArrayToDouble(Base64.getDecoder().decode(s));
		return new Scaling(scaling);
		
	}
	
	

	
}
