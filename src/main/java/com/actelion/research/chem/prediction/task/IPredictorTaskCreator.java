package com.actelion.research.chem.prediction.task;

public interface IPredictorTaskCreator {
	
	public IPredictorTask resolve(String taskParameters);

}
