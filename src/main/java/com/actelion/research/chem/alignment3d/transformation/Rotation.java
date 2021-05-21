package com.actelion.research.chem.alignment3d.transformation;

import java.util.Base64;
import java.util.Base64.Encoder;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.phesa.EncodeFunctions;

public class Rotation implements Transformation {
	
	private static final String DELIMITER = ";";

	public double[][] getRotation() {
		return m;
	}

	private double[][] m;
	
	public Rotation(double[][] m) {
		this.m = m;
	}
	
	public Rotation() {
		this(new double[3][3]);
	}

	public void setRotation(double[][] m) {
		this.m = m;
	}
	
	public Rotation getInvert() {
		double[][] invert = new double[3][3];
        for (int ii = 0; ii < 3; ii++) {
            for (int jj = 0; jj < 3; jj++) {
                invert[jj][ii] = m[ii][jj];
            }
        }
		return new Rotation(invert);
	}

	@Override
	public void apply(Coordinates coords) {
		coords.rotate(m);
		
	}

	@Override
	public void apply(double[] coords) {
		assert(coords.length==3);
		Coordinates c = new Coordinates(coords[0],coords[1],coords[2]);
		apply(c);
		coords[0] = c.x;
		coords[1] = c.y;
		coords[2] = c.z;
	}

	@Override
	public String encode() {
		Encoder encoder = Base64.getEncoder();
		StringBuilder sb = new StringBuilder();
		sb.append(Transformation.ROTATION);
		sb.append(Transformation.DELIMITER);
		for(double[] d : m) {
			sb.append(encoder.encodeToString(EncodeFunctions.doubleArrayToByteArray(d)));
			sb.append(DELIMITER);
		}
		return sb.toString();
	}
	
	public static Rotation decode(String s) {
		String[] str = s.split(DELIMITER);
		double[] r0 = EncodeFunctions.byteArrayToDoubleArray(Base64.getDecoder().decode(str[0]));
		double[] r1 = EncodeFunctions.byteArrayToDoubleArray(Base64.getDecoder().decode(str[1]));
		double[] r2 = EncodeFunctions.byteArrayToDoubleArray(Base64.getDecoder().decode(str[2]));
		double[][] rotation = new double[3][3];
		rotation[0] = r0;
		rotation[1] = r1;
		rotation[2] = r2;
		return new Rotation(rotation);
		
	}

}
