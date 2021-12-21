package com.actelion.research.chem.prediction.task;

public interface IPredictorTask {
	/**
	 * takes an idcode as an input and calculates a property (return an array since certain tasks 
	 * may return several values (e.g. property plus confidence interval for a prediction task, or
	 * docking score plus the pose)
	 * @param idcode
	 * @return
	 */
	public String[] run(String idcode); 
	

}
