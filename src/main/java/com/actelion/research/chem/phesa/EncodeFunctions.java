package com.actelion.research.chem.phesa;

import java.nio.ByteBuffer;

/**
 * February 2018
 * Author: J. Wahl
 * helper functions for encoding shape properties using a Base64 encoder
 **/

public class EncodeFunctions {
	
	
	public static byte[] doubleToByteArray(double value) {
		byte[] bytes = new byte[8];
		ByteBuffer buffer = ByteBuffer.allocate(bytes.length);
		buffer.putDouble(value);
		return buffer.array();
		
	}
	
	public static byte[] doubleArrayToByteArray(double[] doubleArray) {
		ByteBuffer buffer = ByteBuffer.allocate(Double.SIZE / Byte.SIZE * doubleArray.length);
		buffer.asDoubleBuffer().put(doubleArray);
		return buffer.array();
		
	}
	
	public static byte[] intArrayToByteArray(int[] intArray) {
		ByteBuffer buffer = ByteBuffer.allocate(Integer.SIZE / Byte.SIZE * intArray.length);
		buffer.asIntBuffer().put(intArray);
		return buffer.array();
		
	}
	
	public static byte[] intToByteArray(int value) {
		byte[] bytes = new byte[4];
		ByteBuffer buffer = ByteBuffer.allocate(bytes.length);
		buffer.putInt(value);
		return buffer.array();
	}
	
	public static double byteArrayToDouble(byte[] bytes) {
		return ByteBuffer.wrap(bytes).getDouble();
	}
	
	public static double[] byteArrayToDoubleArray(byte[] bytes) {
		int times = Double.SIZE / Byte.SIZE;
		double[] doubles = new double[bytes.length / times];
		for(int i=0; i<doubles.length;i++) {
			doubles[i] = ByteBuffer.wrap(bytes, i*times, times).getDouble();
		}
		return doubles;
	}
	
	public static int[] byteArrayToIntArray(byte[] bytes) {
		int times = Integer.SIZE / Byte.SIZE;
		int[] ints = new int[bytes.length / times];
		for(int i=0; i<ints.length;i++) {
			ints[i] = ByteBuffer.wrap(bytes, i*times, times).getInt();
		}
		return ints;
	}
	
	public static int byteArrayToInt(byte[] bytes) {
		return ByteBuffer.wrap(bytes).getInt();
	}

}
