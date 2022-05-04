package com.actelion.research.gui.generic;

public interface GenericImage {
	void scale(int width, int height);
	Object get();
	int getWidth();
	int getHeight();
	void setRGB(int x, int y, int argb);
	}
