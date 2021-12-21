package com.actelion.research.chem.prediction.task;

public class SimplePropertyTaskCreator implements IPredictorTaskCreator {
	
	public static final String TASK_CODE = "simpleProperty";

	@Override
	public IPredictorTask resolve(String propertyCode) {
		return new SimplePropertyTask(propertyCode);
	}

}
