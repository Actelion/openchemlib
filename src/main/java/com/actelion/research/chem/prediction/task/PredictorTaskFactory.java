package com.actelion.research.chem.prediction.task;

import java.util.HashMap;
import java.util.Map;

public class PredictorTaskFactory {
	
	public Map<String,IPredictorTaskCreator> taskCreators;
	
	public PredictorTaskFactory() {
		this.taskCreators = new HashMap<String,IPredictorTaskCreator>();
		taskCreators.put(SimplePropertyTaskCreator.TASK_CODE, new SimplePropertyTaskCreator()); //add task creator for simple property tasks
		
	}
	
	public void addTaskCreator(String taskName, IPredictorTaskCreator taskCreator) {
		taskCreators.put(taskName, taskCreator);
	}
	
	public IPredictorTask resolve(String taskName, String taskParameters) throws TaskResolveException {
		IPredictorTaskCreator taskCreator = taskCreators.get(taskName);
		IPredictorTask task = null;
		if(taskCreator == null)
			throw new TaskResolveException();
		else {
			try {
				task = taskCreator.resolve(taskParameters);
			}
			catch(Exception e) {
				throw new TaskResolveException();
			}
		}
		return task;
	}
	
	public static class TaskResolveException extends Exception {
		private static final long serialVersionUID = 1L;
		private static final String ERROR_MESSAGE = "Could not resolve task";
		
		public TaskResolveException() {
			super(ERROR_MESSAGE);
		}
	}

}
