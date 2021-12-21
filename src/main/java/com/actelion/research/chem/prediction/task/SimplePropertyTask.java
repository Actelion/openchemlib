package com.actelion.research.chem.prediction.task;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.prediction.MolecularPropertyHelper;

public class SimplePropertyTask implements IPredictorTask {

	private String propertyCode;
	
	public SimplePropertyTask(String propertyCode) {
		this.propertyCode = propertyCode;
	}

	@Override
	public String[] run(String idcode) {
		IDCodeParser parser = new IDCodeParser();
		StereoMolecule mol = new StereoMolecule();
		parser.parse(mol, idcode);;
		int type = MolecularPropertyHelper.getTypeFromCode(propertyCode);
		float value = MolecularPropertyHelper.calculateProperty(mol, type);
		String[] propertyArray = new String[1];
		propertyArray[0] = Float.toString(value);
		return propertyArray;
	}
	
	
	
	
}
