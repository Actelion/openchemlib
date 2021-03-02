package com.actelion.research.chem.descriptor.pharmacophoretree;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorConstants;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorInfo;
import com.actelion.research.chem.descriptor.pharmacophoretree.PharmacophoreNode;
import com.actelion.research.chem.descriptor.pharmacophoretree.PharmacophoreTree;
import com.actelion.research.chem.descriptor.pharmacophoretree.PharmacophoreTreeGenerator;
import com.actelion.research.chem.descriptor.pharmacophoretree.TreeMatcher;

public class DescriptorHandlerPTree implements DescriptorHandler<PharmacophoreTree,StereoMolecule> {
	
	private static final int MAX_NODE_SIZE = 40;

	private static DescriptorHandlerPTree INSTANCE;
	
	public static final PharmacophoreTree FAILED_OBJECT = new PharmacophoreTree(new ArrayList<PharmacophoreNode>(), 
			new ArrayList<int[]>());

	
	@Override
	public float getSimilarity(PharmacophoreTree pt1, PharmacophoreTree pt2) {
		TreeMatcher matcher = new TreeMatcher(pt1,pt2);
		TreeMatcher.TreeMatching matching = matcher.matchSearch();
		return (float)matching.getSim();
	}
	
	@Override
	public String getVersion() {
		return DescriptorConstants.DESCRIPTOR_PTREE.version;
	}
	
	@Override
	public DescriptorInfo getInfo() {
		return DescriptorConstants.DESCRIPTOR_PTREE;
	}

	@Override
	public String encode(PharmacophoreTree pTree) {
		StringBuilder sb = new StringBuilder();
		for(PharmacophoreNode ppNode : pTree.getNodes()) {
			sb.append(PharmacophoreNode.encode(ppNode));
			sb.append(";");
		}
		sb.deleteCharAt(sb.length()-1); // remove trailing ;
		sb.append(":");
		for(int[] edge: pTree.getEdges()) {
			sb.append(edge[0]);
			sb.append(" ");
			sb.append(edge[1]);
			sb.append(":");
		}
		sb.deleteCharAt(sb.length()-1);
		return sb.toString();
	}

	private PharmacophoreTree getDecodedObject(String s) {
		String[] splitString = s.split(":");
		String[] nodeStrings = splitString[0].split(";");
		List<PharmacophoreNode> nodes = new ArrayList<PharmacophoreNode>();
		List<int[]> edges = new ArrayList<int[]>();
		for(String nodeString : nodeStrings)
			nodes.add(PharmacophoreNode.decode(nodeString));
		for(int i=1;i<splitString.length;i++) {
			String[] edgeStrings = splitString[i].split(" ");
			edges.add(new int[] {Integer.valueOf(edgeStrings[0]),Integer.valueOf(edgeStrings[1])});
		}
		
		
			
		
		return new PharmacophoreTree(nodes,edges);
	}
	
	
	public PharmacophoreTree decode(byte[] arr) {

		return decode(new String(arr));
		
	}
	
	@Override
	public PharmacophoreTree decode(String s) {
		try {
			return s == null || s.length() == 0 ? null
					: s.equals(FAILED_STRING) ? FAILED_OBJECT
					:                           getDecodedObject(s);
		} catch (RuntimeException e1) {
			return FAILED_OBJECT;
		}
	}

	@Override
	public PharmacophoreTree createDescriptor(StereoMolecule mol) {
		mol.stripSmallFragments();
		
		PharmacophoreTree pharmTree =  PharmacophoreTreeGenerator.generate(mol);
		if(pharmTree.getNodes().size()>MAX_NODE_SIZE)
			pharmTree = FAILED_OBJECT;
		return pharmTree;

	}

	@Override
	public boolean calculationFailed(PharmacophoreTree o) {
		boolean failed = o.getNodes().size() == 0?  true : false;
		return failed;
	}

	@Override
	public DescriptorHandler<PharmacophoreTree, StereoMolecule> getThreadSafeCopy() {
		return new DescriptorHandlerPTree();
	}


	public static DescriptorHandlerPTree getDefaultInstance(){

		if(INSTANCE==null){
			INSTANCE = new DescriptorHandlerPTree();
		}

		return INSTANCE;
	}

	

}
