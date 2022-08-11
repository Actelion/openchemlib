package com.actelion.research.gui.generic;

public interface GenericDrawContext {
	float getLineWidth();
	void setLineWidth(float lineWidth);
	int getRGB();
	void setRGB(int rgb);
	int getFontSize();
	void setFont(int size, boolean isBold, boolean isItalic);
	void drawLine(double x1, double y1, double x2, double y2);
	void drawDottedLine(double x1, double y1, double x2, double y2);
	void drawRectangle(double x, double y, double w, double h);
	void fillRectangle(double x, double y, double w, double h);
	void drawCircle(double x, double y, double d);
	void fillCircle(double x, double y, double d);
	void drawPolygon(GenericPolygon p);
	void fillPolygon(GenericPolygon p);
	void drawString(double x, double y, String s);
	void drawCenteredString(double x, double y, String s);
	void drawImage(GenericImage image, double x, double y);
	void drawImage(GenericImage image, double sx, double sy, double dx, double dy, double w, double h);
	void drawImage(GenericImage image, double sx, double sy, double sw, double sh, double dx, double dy, double dw, double dh);
//	void setClip(double x, double y, double w, double h);
	GenericRectangle getBounds(String s);
	boolean isDarkBackground();
	int getForegroundRGB();
	int getBackgroundRGB();
	int getSelectionBackgroundRGB();
	GenericImage createARGBImage(int width, int height);
}
