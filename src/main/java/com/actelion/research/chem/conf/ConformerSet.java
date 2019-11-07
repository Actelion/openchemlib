package com.actelion.research.chem.conf;


import java.util.Iterator;
import java.util.TreeSet;
import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParserWithoutCoordinateInvention;
import com.actelion.research.chem.StereoMolecule;


public class ConformerSet extends TreeSet<Conformer> {
	
	public ConformerSet() {
		this(null);
	}
	
	public ConformerSet(String s) {
		super();
		if (s != null) {
			String[] confString = s.split(" ");
			if (confString.length >= 2) {
				IDCodeParserWithoutCoordinateInvention parser = new IDCodeParserWithoutCoordinateInvention();
				StereoMolecule mol = parser.getCompactMolecule(confString[0], confString[1]);
				Conformer firstConformer = new Conformer(mol);
				add(firstConformer);

				int[] atomAndBondCount = parser.getAtomAndBondCounts(confString[0], null);
				for(int i=2; i<confString.length;i++) {
					Conformer conf = new Conformer(firstConformer);
					try {
						parser.parseCoordinates(confString[i].getBytes(), 0, conf, atomAndBondCount[0], atomAndBondCount[1]);
						add(conf);
					} catch (Exception e) {}
				}
			}
		}
	}
	
	
	public ConformerSet getSubset(int size) {
		ConformerSet treeSet = new ConformerSet();
		this.stream().limit(size).forEach(c->treeSet.add(c));
		return treeSet;
		
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		if (!isEmpty()) {
			Iterator<Conformer> iterator = iterator();
			StereoMolecule mol = iterator.next().toMolecule();
			Canonizer can = new Canonizer(mol, Canonizer.COORDS_ARE_3D);
			sb.append(can.getIDCode());
			sb.append(" ");
			sb.append(can.getEncodedCoordinates(true));
			while (iterator.hasNext()) {
				can.invalidateCoordinates();
				iterator.next().toMolecule(mol);
				sb.append(" ");
				sb.append(can.getEncodedCoordinates(true));
			}
		}
		return sb.toString();
	}
}
