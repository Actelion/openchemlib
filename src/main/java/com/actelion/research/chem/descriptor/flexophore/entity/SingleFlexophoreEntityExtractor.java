package com.actelion.research.chem.descriptor.flexophore.entity;

import com.actelion.research.calc.ArrayUtilsCalc;
import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.ConstantsDWAR;
import com.actelion.research.chem.descriptor.flexophore.Molecule3DExtractor;
import com.actelion.research.chem.descriptor.flexophore.Molecule3DFunctions;
import com.actelion.research.chem.descriptor.flexophore.MolDistHistViz;
import com.actelion.research.chem.descriptor.flexophore.PPNodeViz;
import com.actelion.research.chem.descriptor.flexophore.UnparametrizedAtomTypeException;
import com.actelion.research.chem.descriptor.flexophore.generator.FlexophoreCreateFunctions;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;


/**
 * SingleFlexophoreEntityExtractor
 * Extract pharmacophore points and linker.
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Jan 14, 2013 MvK Start implementation
 */
public class SingleFlexophoreEntityExtractor {

	public static final int OFFSET_LABEL_RGROUP_FUSION = 10;
	
	private String idcode; 
	
	private String coordinates;
	
	private MolDistHistViz mdhv;
	
	private Molecule3D ff;
		
	private List<FlexophorePoint> liFlexophorePoint;
	
	private List<Linker> liLinker;
	
	private int linkerIdCounter;
	
	
	public SingleFlexophoreEntityExtractor(String idcode, String coordinates) {
		
		this.idcode = idcode;
		
		this.coordinates = coordinates;
		
		liFlexophorePoint = new ArrayList<FlexophorePoint>();
	}
		
	public void process() throws UnparametrizedAtomTypeException {
		
		mdhv = FlexophoreCreateFunctions.create(idcode, coordinates);

		ff = mdhv.getMolecule();
		
		int sizeFlexophore = mdhv.getNumPPNodes();
		
		for (int i = 0; i < sizeFlexophore; i++) {
			PPNodeViz node = mdhv.getNode(i);
			
			FlexophorePoint fp =  getFlexophorePoint(node, ff, i);
							
			liFlexophorePoint.add(fp);
			
		}
		
		extractLinker();
		
		for(FlexophorePoint fp2ReArrange : liFlexophorePoint){
			// addIdCode2FlexophorePointWithRGroups(fp2ReArrange, ff);
			addIdCode2FlexophorePoint(fp2ReArrange, ff);
		}
		
	}
	
	
	public List<FlexophorePoint> getPharmacophorePointList(){
		return liFlexophorePoint;
	}
	
	public List<Linker> getLinkerList(){
		return liLinker;
	}
	
	private void extractLinker(){
		
		liLinker = new ArrayList<Linker>();
		
		for (int i = 0; i < liFlexophorePoint.size(); i++) {
			
			for (int j = i+1; j < liFlexophorePoint.size(); j++) {
				
				Linker linker = extractLinker(liFlexophorePoint.get(i), liFlexophorePoint.get(j));
				
				if(linker!=null){
					liLinker.add(linker);
				}
			}
		}
	}
	
	private Linker extractLinker(FlexophorePoint ia1, FlexophorePoint ia2){
		
		int [] arrShortestIndexPath = getShortestPathWithoutPharmacophorePointBetween(ia1, ia2);
		
		if(arrShortestIndexPath==null){
			return null;
		}
		
		
		Linker linker = new Linker(linkerIdCounter++);
		
		linker.addOriginalAtomIndex(arrShortestIndexPath);
		
		linker.setDistanceHistogram(mdhv.getDistHist(ia1.getId(), ia2.getId()));
		
		linker.setIdFlexophorePoint1(ia1.getId());
		
		linker.setIdFlexophorePoint2(ia2.getId());
		
		addIdCode2Linker(linker, ff);
		
		
		AtomIndexLinkerId atomIndexLinkerId1 = new AtomIndexLinkerId(arrShortestIndexPath[0], arrShortestIndexPath[1], linker.getId());
		
		ia1.addArrOriginalAtomIndexRGroups(atomIndexLinkerId1);
		
		AtomIndexLinkerId atomIndexLinkerId2 = new AtomIndexLinkerId(arrShortestIndexPath[arrShortestIndexPath.length-1], arrShortestIndexPath[arrShortestIndexPath.length-2], linker.getId());
		
		ia2.addArrOriginalAtomIndexRGroups(atomIndexLinkerId2);
		
		return linker;
	}
	
	
	private FlexophorePoint getFlexophorePoint(PPNodeViz node, Molecule3D ff, int id){
		
		FlexophorePoint flexophorePoint = new FlexophorePoint(id);
		
		flexophorePoint.setInteractionType(node.get());
		
		List<Integer> liIndexAtomOrig = node.getListIndexOriginalAtoms();
				
		//
		// Fill gaps between atoms belonging to the pharmacophore point.
		//
		HashSet<Integer> hsAtomIndex = new HashSet<Integer>(liIndexAtomOrig);
		for (int i = 0; i < liIndexAtomOrig.size(); i++) {
			
			for (int j = i+1; j < liIndexAtomOrig.size(); j++) {
				int [] arrIndex = Molecule3DFunctions.getPath(ff, liIndexAtomOrig.get(i), liIndexAtomOrig.get(j));
				
				for (int k = 0; k < arrIndex.length; k++) {
					hsAtomIndex.add(arrIndex[k]);
				}
			}
		}
				
		List<Integer> liIndexAtomExpanded = new ArrayList<Integer>(hsAtomIndex);
		for (int atomIndex : liIndexAtomExpanded) {
			flexophorePoint.addOriginalAtomIndex(atomIndex);
		}
		
		return flexophorePoint;
	}
	

	/**
	 * Fused Flexophore points are indicated with a R label > <code>OFFSET_LABEL_RGROUP_FUSION</code>. 
	 * @param fp
	 * @param ffOriginal
	 */
	private static void addIdCode2FlexophorePointWithRGroups(FlexophorePoint fp, Molecule3D ffOriginal){
		
		List<Integer> liIndexAtomWithRGroups = new ArrayList<Integer>(fp.getOriginalAtomIndex());
		
		List<AtomIndexLinkerId> liAtomIndexOriginalRGroups = fp.getAtomIndexLinkerId();
		
		Molecule3D molRGroups = new Molecule3D(ffOriginal);
		
		int arrayIndexAtomLabelRGroup = 0;
		
		int arrayIndexAtomLabelFusionRGroup = OFFSET_LABEL_RGROUP_FUSION;
		
		for(AtomIndexLinkerId atomIndexLinkerId : liAtomIndexOriginalRGroups){
			
			int indexAtomFlexPoint = atomIndexLinkerId.getAtomIndexConnectionFlexophorePoint();
			
			int indexAtomRGroup = atomIndexLinkerId.getAtomIndexFirstLinkerAtom();
						
			int rGroupLabel = -1;
			
			int bondOrder = -1;
			
			if(indexAtomFlexPoint == indexAtomRGroup){ // This were two fused pharmacophore points
				
				rGroupLabel = ConstantsDWAR.arrAtomLabelRGroups[arrayIndexAtomLabelFusionRGroup++];
				
				bondOrder = Molecule.cBondTypeSingle;
								
			} else {
				rGroupLabel = ConstantsDWAR.arrAtomLabelRGroups[arrayIndexAtomLabelRGroup++];
				
				int indexBond = molRGroups.getBond(indexAtomFlexPoint, indexAtomRGroup);
				
				bondOrder = molRGroups.getBondOrder(indexBond);
			}
						
			int indexAtomRGroupNew = molRGroups.addAtom(rGroupLabel);
			
			int connected = molRGroups.getConnAtoms(indexAtomFlexPoint);
			
			if(connected==8){
				System.out.println("SingleFlexophoreEntityExtractor addIdCode2FlexophorePoint(...) atom " + indexAtomFlexPoint + " has already 8 connections:");
				
				StereoMolecule mol = new Molecule3D(molRGroups);
				
				mol.ensureHelperArrays(Molecule.cHelperRings);
				
				Canonizer can = new Canonizer(mol);
				
				System.out.println(can.getIDCode());
				
				System.out.println();
			} else {
				molRGroups.addBond(indexAtomFlexPoint, indexAtomRGroupNew, bondOrder);
				
				liIndexAtomWithRGroups.add(indexAtomRGroupNew);
			}
			
			
			
			// Reset, start with R1 again.
			if(arrayIndexAtomLabelRGroup == OFFSET_LABEL_RGROUP_FUSION){
				arrayIndexAtomLabelRGroup = 0;
			}
			// Reset, start with R OFFSET_LABEL_RGROUP_FUSION again.
			if(arrayIndexAtomLabelFusionRGroup == ConstantsDWAR.arrAtomLabelRGroups.length){
				arrayIndexAtomLabelFusionRGroup = 0;
			}
		}

		Molecule3D mol = Molecule3DExtractor.extract(molRGroups, liIndexAtomWithRGroups);

		mol.ensureHelperArrays(Molecule.cHelperRings);
		
		Canonizer can = new Canonizer(mol);
		
		fp.setIdCode(can.getIDCode());

	}
	
	private static void addIdCode2FlexophorePoint(FlexophorePoint fp, Molecule3D ffOriginal){
		
		List<Integer> liIndexAtomWithRGroups = new ArrayList<Integer>(fp.getOriginalAtomIndex());
		
		Molecule3D ffRGroups = new Molecule3D(ffOriginal);
		
		Molecule3D ffSub = Molecule3DExtractor.extract(ffRGroups, liIndexAtomWithRGroups);
		
		StereoMolecule mol = new Molecule3D(ffSub);
		
		mol.setFragment(true);
		
		mol.ensureHelperArrays(Molecule.cHelperRings);
		
		Canonizer can = new Canonizer(mol);
		
		fp.setIdCode(can.getIDCode());

	}
	
	/**
	 * If the linker links two fused Flexophore points than no idcode is set, because there are no linking atoms.
	 * @param linker
	 * @param ff
	 */
	private static void addIdCode2Linker(Linker linker, Molecule3D ff){
		
		int [] arrPath = linker.getOriginalAtomIndex();
		
		//
		// Fusion?
		//
		if(arrPath.length==2){
			if(arrPath[0]==arrPath[1]){
				return;
			}
		}
		
		List<Integer> liIndexAtom = ArrayUtilsCalc.toList(arrPath);

		Molecule3D ffSub = Molecule3DExtractor.extract(ff, liIndexAtom);
		
		StereoMolecule mol = new Molecule3D(ffSub);
		
		mol.ensureHelperArrays(Molecule.cHelperRings);
		
		mol.setAtomicNo(0, ConstantsDWAR.ATOM_LABEL_R1);
		
		mol.setAtomicNo(liIndexAtom.size()-1, ConstantsDWAR.ATOM_LABEL_R2);
				
		Canonizer can = new Canonizer(mol);

		linker.setIdCode(can.getIDCode());
	}

	
	
	/**
	 * Start index is taken from ia1 and the end index is taken from ia2. 
	 * If both Flexophore points were fused the path has a length of two and contains two times the same atom index.
	 * @param ia1
	 * @param ia2
	 * @return
	 */
	
	private int [] getShortestPathWithoutPharmacophorePointBetween(FlexophorePoint ia1, FlexophorePoint ia2){
		
		HashSet<Integer> hsAtomIndexPPPoints = new HashSet<Integer>();
		
		for (int i = 0; i < liFlexophorePoint.size(); i++) {
			
			FlexophorePoint fp = liFlexophorePoint.get(i);
			
			if((fp.getId() != ia1.getId()) && (fp.getId() != ia2.getId())){
				
				List<Integer> liAtomIndexOriginal = fp.getOriginalAtomIndex();
				
				
				hsAtomIndexPPPoints.addAll(liAtomIndexOriginal);
				
			}
		}
		
		List<Integer> liOrigIndex1 = ia1.getOriginalAtomIndex();
		
		List<Integer> liOrigIndex2 = ia2.getOriginalAtomIndex();
		
		List<int []> liValidPathList = new ArrayList<int[]>();
		
		for (int i = 0; i < liOrigIndex1.size(); i++) {
			
			int atomIndexFlexPoint1 = liOrigIndex1.get(i);
			
			for (int j = 0; j < liOrigIndex2.size(); j++) {
				
				int atomIndexFlexPoint2 = liOrigIndex2.get(j);
				
				int [] arrIndexPath = null;
				
				if(atomIndexFlexPoint1 == atomIndexFlexPoint2){ // Two fused Flexophore points?
					arrIndexPath = new int [2];
					arrIndexPath[0]=atomIndexFlexPoint1;
					arrIndexPath[1]=atomIndexFlexPoint2;
				} else {
					arrIndexPath = Molecule3DFunctions.getPath(ff, atomIndexFlexPoint1, atomIndexFlexPoint2);
				}
				
				boolean pathContainsPPPoint = false;
								
				for (int k = 0; k < arrIndexPath.length; k++) {
					if(hsAtomIndexPPPoints.contains(arrIndexPath[k])){
						pathContainsPPPoint = true;
						break;
					}
				}
								
				if(!pathContainsPPPoint){
					
					liValidPathList.add(arrIndexPath);
				}
			}
		}
		
		if(liValidPathList.size()==0){
			return null;
		}
		
		
		int [] arrShortestIndexPath = liValidPathList.get(0);
		
		for (int i = 0; i < liValidPathList.size(); i++) {
			if(liValidPathList.get(i).length<arrShortestIndexPath.length){
				arrShortestIndexPath = liValidPathList.get(i); 
			}
		}
		
		return arrShortestIndexPath;
	}
	


}
