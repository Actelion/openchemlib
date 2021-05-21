package com.actelion.research.chem.alignment3d.transformation;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

import java.util.Base64;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;



public interface Transformation {
	
	public static final String ROTATION = "r";
	public static final String TRANSLATION = "t";
	public static final String SCALING = "s";
	public static final String DELIMITER = " ";
	
	public void apply(Coordinates coords);
	
	public void apply(double[] coords);
	
	public default void apply(StereoMolecule mol) {
		for(int a=0;a<mol.getAllAtoms();a++) {
			apply(mol.getCoordinates(a));
		}
	}
	
	public default void apply(Conformer conf) {
		for(int a=0;a<conf.getMolecule().getAllAtoms();a++) {
			apply(conf.getCoordinates(a));
		}
	}
	
	public String encode();
	
	public static Transformation decode(String s) {
		if(s.startsWith(Transformation.ROTATION)) {
			return Rotation.decode(s);
		}
		else if(s.startsWith(Transformation.TRANSLATION)) {
			return Translation.decode(s);
		}
		else if(s.startsWith(Transformation.SCALING)) {
			return Scaling.decode(s);
		}
		else if(s.startsWith(String.valueOf(TransformationSequence.DELIMITER_START))) {
			return TransformationSequence.decode(s);
		}
		else
			return null;
}
}
