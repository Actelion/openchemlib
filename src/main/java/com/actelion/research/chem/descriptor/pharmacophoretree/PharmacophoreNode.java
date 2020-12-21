package com.actelion.research.chem.descriptor.pharmacophoretree;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.actelion.research.chem.phesa.EncodeFunctions;


/**
 * A PharmacophoreNode represents a node in a PharmacophoreTree. A PharmacophoreNode is associated with
 * a list of features (hbond donor, hbond acceptor, +charge,-charge,aromatic, lipophilic) which define its chemical
 * properties. The similarity between two pharmacophore nodes is calculated by using a combination of steric
 * and chemical similarity.
 * @author joel
 *
 */

public class PharmacophoreNode {

	
	public static final int[] FUNCTIONALITY_WEIGHTS = new int[] {3,3,3,3,1,1};
	public static final double CHEM_SIM_WEIGHT = 0.7;
	public static final int ZERO_NODE = 1;
	public static final int LINK_NODE = 6;
	

	private int[] functionalities;
	private List<Integer> atoms;
	private List<Double> weights;
	private List<Double> volumes;
	private double size;
	private double vol;
	private int role;
	private boolean isRing;
	private boolean isAromatic;
	
	public PharmacophoreNode(List<Integer> atoms, int[][] atomFunctionalities, double[] atomVolumes, int role, boolean isRing, boolean isAromatic) {
		this.functionalities = new int[FUNCTIONALITY_WEIGHTS.length];
		this.atoms = atoms;
		this.isRing = isRing;
		this.isAromatic = isAromatic;
		weights = new ArrayList<Double>();
		volumes = new ArrayList<Double>();
		this.role = role;
		if((role  & ZERO_NODE)==0) {
			this.atoms.stream().forEach(a -> {
				weights.add(1.0);
				volumes.add(atomVolumes[a]);
			});
			this.atoms.stream().forEach(a -> {
				int[] f = atomFunctionalities[a];
				IntStream.range(0, f.length).forEach(i -> functionalities[i]+=f[i]);
			});
			calculate();
		}


	}
	
	public PharmacophoreNode(List<Integer> atoms, int[] functionalities, List<Double> volumes,
			List<Double> weights, int role) {
		this.atoms = atoms;
		this.functionalities = functionalities;
		this.volumes = volumes;
		this.weights = weights;
		this.role = role;
		calculate();
		
	}
	
	public PharmacophoreNode(List<Integer> atoms, int[][] atomFunctionalities, double[] atomVolumes, boolean isRing, boolean isAromatic) {
		this(atoms, atomFunctionalities, atomVolumes, 0, isRing, isAromatic);
	}
	
	public void updateWeights(Map<Integer,List<Integer>> atomToNodes) {
		if((role  & ZERO_NODE)==0) {
			IntStream.range(0, atoms.size()).forEach((e) -> {
				int a = atoms.get(e);
				double weight = 1.0/atomToNodes.get(a).size();
				weights.set(e, weight);
			});
			calculate();
		}
			
	}
	
	/**
	 * calculates the total volume and size of the node from its elements
	 */
	public void calculate() {
		size = 0.0;
		vol = 0.0;
		for(double w : weights) {
			size+=w;
		}
		vol = IntStream.range(0, volumes.size()).mapToDouble(e -> 
		volumes.get(e)*weights.get(e)).sum();
	}
	

	
	private static double calcStericSim ( double sum1, double sum2) {
		double sim = 0.0;
		if(sum1+sum2<0.001) // zero-nodes
			sim = 1.0;
		else {
			sim = 2*Math.min(sum1, sum2)/(sum1+sum2);
		}
		return sim;
		
	}

	
	public static double calcFeatureSim(int[] features1, int[] features2) {
		double sim = 0.0;
		double nom = 0.0;
		double denom = 0.0;
		for(int i=0;i<features1.length;i++) {
			int weight = FUNCTIONALITY_WEIGHTS[i];
			nom+=weight*Math.min(features1[i], features2[i]);
			denom+=weight*(features1[i] + features2[i]);
		}
		if(denom!=0) {
			sim = 2*nom/denom;
		}
		return sim;	
	}
	
	
	public static double getSimilarity(Collection<Integer> nodes1,Collection<Integer> nodes2, List<PharmacophoreNode> allNodes1,
			List<PharmacophoreNode> allNodes2) {
		double size1 = 0.0;
		double vol1 = 0.0;
		double sterSim = 0.0;
		double chemSim = 0.0;
		int[] functionalities1 = new int[FUNCTIONALITY_WEIGHTS.length];

		for(int n : nodes1) {
			PharmacophoreNode node = allNodes1.get(n);
			size1+=node.size;
			vol1+=node.vol;
			for(int i=0;i<functionalities1.length;i++) {
				functionalities1[i]+=node.functionalities[i];
			}
				
		}
		double size2 = 0.0;
		double vol2 = 0.0;
		int[] functionalities2 = new int[FUNCTIONALITY_WEIGHTS.length];
		for(int n : nodes2) {
			PharmacophoreNode node = allNodes2.get(n);
			size2+=node.size;
			vol2+=node.vol;
			for(int i=0;i<functionalities2.length;i++) {
				functionalities2[i]+=node.functionalities[i];
			}
		}
		if(size1/size2 > TreeMatcher.SIZE_RATIO || size1/size2 < (1.0/TreeMatcher.SIZE_RATIO )) {
			sterSim = 0.0;
			chemSim = 0.0;
		}
		else{
			sterSim = 0.5*calcStericSim(size1,size2)+0.5*calcStericSim(vol1,vol2);
			chemSim = calcFeatureSim(functionalities1,functionalities2);
		}
		if(nodes1.size() == 1 && nodes2.size()==1) { //matching link Nodes have similarity of 1;
			PharmacophoreNode n1 = allNodes1.get((Integer)nodes1.toArray()[0]);
			PharmacophoreNode n2 = allNodes2.get((Integer)nodes2.toArray()[0]);
			if(n1.isLinkNode() && n2.isLinkNode()) {
				sterSim = 1.0;
				chemSim = 1.0;
			}
		}

		return (1.0-CHEM_SIM_WEIGHT)*sterSim+CHEM_SIM_WEIGHT*chemSim;	
		
	}
	
	
	
	
	public static double getSimilarity(Collection<PharmacophoreNode> nodes1,Collection<PharmacophoreNode> nodes2) {
		double size1 = 0.0;
		double vol1 = 0.0;
		double sterSim = 0.0;
		double chemSim = 0.0;
		int[] functionalities1 = new int[FUNCTIONALITY_WEIGHTS.length];

		for(PharmacophoreNode node : nodes1) {
			size1+=node.size;
			vol1+=node.vol;
			for(int i=0;i<functionalities1.length;i++) {
				functionalities1[i]+=node.functionalities[i];
			}
				
		}
		double size2 = 0.0;
		double vol2 = 0.0;
		int[] functionalities2 = new int[FUNCTIONALITY_WEIGHTS.length];
		for(PharmacophoreNode node : nodes2) {
			size2+=node.size;
			vol2+=node.vol;
			for(int i=0;i<functionalities2.length;i++) {
				functionalities2[i]+=node.functionalities[i];
			}
		}
		if(size1/size2 > TreeMatcher.SIZE_RATIO || size1/size2 < (1.0/TreeMatcher.SIZE_RATIO )) {
			sterSim = 0.0;
			chemSim = 0.0;
		}
		else{
			sterSim = 0.5*calcStericSim(size1,size2)+0.5*calcStericSim(vol1,vol2);
			chemSim = calcFeatureSim(functionalities1,functionalities2);
		}
		if(nodes1.size() == 1 && nodes2.size()==1) { //matching link Nodes have similarity of 1;
			if(((PharmacophoreNode)nodes1.toArray()[0]).isLinkNode() && ((PharmacophoreNode)nodes2.toArray()[0]).isLinkNode()) {
				sterSim = 1.0;
				chemSim = 1.0;
			}
		}
	
		return (1.0-CHEM_SIM_WEIGHT)*sterSim+CHEM_SIM_WEIGHT*chemSim;	
		
	}
	
	public List<Integer> getAtoms() {
		return atoms;
	}
	
	public List<Double> getWeights() {
		return weights;
	}
	
	public int[] getFunctionalities() {
		return functionalities;
	}
	
	public List<Double> getVolumes() {
		return weights;
	}
	
	public double getSize() {
		return size;
	}
	
	public void setFunctionalities(int[] functionalities) {
		this.functionalities = functionalities;
	}
	
	public void setVolumes(List<Double> volumes) {
		this.volumes = volumes;
	}
	
	public void setWeights(List<Double> weights) {
		this.weights = weights;
	}
	
	public void setAtoms(List<Integer> atoms) {
		this.atoms = atoms;
	}
	
	public void setRole(int role) {
		this.role = role;
	}
	
	
	
	public boolean isLinkNode() {
		 boolean b = (role & LINK_NODE)==0 ? false : true;
		 return b;
	}
		
	public static String encode(PharmacophoreNode node) {
		String s1 = Base64.getEncoder().encodeToString(EncodeFunctions.intArrayToByteArray(node.functionalities));
		String s2 =  Base64.getEncoder().encodeToString(EncodeFunctions.intArrayToByteArray(node.atoms.stream().
				mapToInt(i -> i.intValue()).toArray()));
		String s3 =  Base64.getEncoder().encodeToString(EncodeFunctions.doubleArrayToByteArray(node.volumes.stream().
				mapToDouble(i -> i.doubleValue()).toArray()));
		String s4 =  Base64.getEncoder().encodeToString(EncodeFunctions.doubleArrayToByteArray(node.weights.stream().
				mapToDouble(i -> i.doubleValue()).toArray()));
		String s5 = Integer.toString(node.role);
		StringBuilder sb = new StringBuilder();
		sb.append(s1);
		sb.append(","); //comma is not part of base64 encoding chars
		sb.append(s2);
		sb.append(",");
		sb.append(s3);
		sb.append(",");
		sb.append(s4);
		sb.append(",");
		sb.append(s5);
		sb.append(",");
		return sb.toString();
	}
	
	public static PharmacophoreNode decode (String s) {
		String[] strings = s.split(",");
		int[] functionalities = EncodeFunctions.byteArrayToIntArray(Base64.getDecoder().decode(strings[0]));
		List<Integer> atoms = Arrays.stream(EncodeFunctions.byteArrayToIntArray(Base64.getDecoder().decode(strings[1])))
				.boxed().collect(Collectors.toList());
		List<Double> volumes = Arrays.stream(EncodeFunctions.byteArrayToDoubleArray(Base64.getDecoder().decode(strings[2])))
				.boxed().collect(Collectors.toList());
		List<Double> weights = Arrays.stream(EncodeFunctions.byteArrayToDoubleArray(Base64.getDecoder().decode(strings[3])))
				.boxed().collect(Collectors.toList());
		int role = Integer.parseInt(strings[4]);
		return new PharmacophoreNode(atoms,functionalities,volumes,weights,role);
	}

	public boolean isRing() {
		return isRing;
	}

	public void setRing(boolean isRing) {
		this.isRing = isRing;
	}

	public boolean isAromatic() {
		return isAromatic;
	}

	public void setAromatic(boolean isAromatic) {
		this.isAromatic = isAromatic;
	}
	
	
	
}
