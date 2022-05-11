package com.actelion.research.gui.generic;

public interface GenericImage {
	void scale(int width, int height);
	Object get();
	int getWidth();
	int getHeight();
	int getRGB(int x, int y);
	void setRGB(int x, int y, int argb);
	}
