package com.actelion.research.chem.alignment3d.transformation;

import java.util.Base64;
import java.util.Base64.Encoder;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.phesa.EncodeFunctions;

public class Translation implements Transformation {
	
	private static final long serialVersionUID = 20211905L;
	
	
	private double[] translation;
	
	public Translation(double[] translation) {
		this.translation = translation;
	}
	
	public Translation(Coordinates translation) {
		this.translation = new double[3];
		this.translation[0] = translation.x;
		this.translation[1] = translation.y;
		this.translation[2] = translation.z;
	}
	
	public Translation(double x, double y, double z) {
		this.translation = new double[3];
		this.translation[0] = x;
		this.translation[1] = y;
		this.translation[2] = z;
	}

	@Override
	public void apply(Coordinates coords) {
		coords.add(translation[0], translation[1], translation[2]);
		
	}

	@Override
	public void apply(double[] coords) {
		coords[0] += translation[0];
		coords[1] += translation[1];
		coords[2] += translation[2];
	}
	
	public static Translation decode(String s) {
		double[] t = EncodeFunctions.byteArrayToDoubleArray(Base64.getDecoder().decode(s));
		return new Translation(t);
		
	}

	@Override
	public String encode() {
		Encoder encoder = Base64.getEncoder();
		StringBuilder sb = new StringBuilder();
		sb.append(Transformation.TRANSLATION);
		sb.append(Transformation.DELIMITER);
		sb.append(encoder.encodeToString(EncodeFunctions.doubleArrayToByteArray(translation)));
		return sb.toString();
	}

}
