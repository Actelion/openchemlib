package com.actelion.research.gui.generic;

public interface GenericCanvas {
	void repaint();
	int getWidth();
	int getHeight();
	GenericDrawContext getDrawContext();
}
