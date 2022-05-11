package com.actelion.research.gui.generic;

public interface GenericCanvas {
	void repaint();
	double getCanvasWidth();
	double getCanvasHeight();
	int getBackgroundRGB();
	GenericDrawContext getDrawContext();
}
