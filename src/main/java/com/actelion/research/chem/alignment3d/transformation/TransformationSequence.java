package com.actelion.research.chem.alignment3d.transformation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.stream.Collectors;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.phesa.EncodeFunctions;

public class TransformationSequence implements Transformation {

	
	public static final char DELIMITER_START = '{';
	public static final char DELIMITER_END = '}';
	private static final char DELIMITER = ',';
	
	List<Transformation> transformationSequence;
	
	public TransformationSequence(Quaternion rotor) {
		this();
		Rotation rot = new Rotation(rotor.getRotMatrix().getArray());
		double normFactor = 1/rotor.normSquared();
		Scaling scaling = new Scaling(normFactor);
		transformationSequence.add(rot);
		transformationSequence.add(scaling);
	}
	
	public TransformationSequence() {
		transformationSequence = new ArrayList<>();
	}
	
	public void addTransformation(Transformation transformation) {
		transformationSequence.add(transformation);
	}
	
	public List<Transformation> getTransformations() {
		return transformationSequence;
	}
	
	public void apply(Coordinates coords) {
		for(Transformation transformation : transformationSequence) {
			transformation.apply(coords);
		}
	}
	
	public void apply(double[] coords) {
		for(Transformation transformation : transformationSequence) {
			transformation.apply(coords);
		}
	}
	
	public void apply(StereoMolecule mol) {
		for(Transformation transformation : transformationSequence) {
			transformation.apply(mol);
		}
	}
	
	private static List<Transformation> processSubstring(List<Character> chars) {
		List<Transformation> transforms = new ArrayList<>();
		String s = chars.stream()              // Stream<Character>
                .map(String::valueOf)   // Stream<String>
                .collect(Collectors.joining());
		String[] st = s.split(String.valueOf(DELIMITER));
		for(String substring : st) {
			String[] substrings = substring.split(Transformation.DELIMITER);
			String identifier = substrings[0];
			if(identifier.equals(Transformation.ROTATION)) {
				transforms.add(Rotation.decode(substrings[1]));
			}
			else if(identifier.equals(Transformation.SCALING)) {
				transforms.add(Scaling.decode(substrings[1]));
			}
			else if(identifier.equals(Transformation.TRANSLATION)) {
				transforms.add(Translation.decode(substrings[1]));
			}

			}
		
		return transforms;
	}
	public static TransformationSequence decode(String s) {
		char[] chars = s.toCharArray();
		List<Character> charsList = new ArrayList<>();
		for(char c : chars) {
			charsList.add(c);
		}
		Iterator<Character> iterator = charsList.iterator();
		iterator.next();
		return decode(iterator);
	}
	
	private static TransformationSequence decode(Iterator<Character> s) {

		TransformationSequence sequence = new TransformationSequence();
		List<Character> chars = new ArrayList<>();
		while(s.hasNext()) {
			char c = s.next();
			if(c==DELIMITER_END) {
				List<Transformation> transforms = processSubstring(chars);
				transforms.stream().forEach(e -> sequence.addTransformation(e));
				chars.clear();
				return sequence;
			}
			else if(c==DELIMITER_START) {
				List<Transformation> transforms = processSubstring(chars);
				transforms.stream().forEach(e -> sequence.addTransformation(e));
				chars.clear();
				sequence.addTransformation(decode(s));
			}
			else {
				chars.add(c);
			}

		}
		return sequence;
		
	}

	@Override
	public String encode() {
		StringBuilder sb = new StringBuilder();
		sb.append(DELIMITER_START);
		for(Transformation t : transformationSequence) {
			sb.append(t.encode());
			sb.append(DELIMITER);
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append(DELIMITER_END);
		return sb.toString();
		
	}
	

}
