package com.actelion.research.gui.generic;

import java.awt.*;
import java.awt.geom.Rectangle2D;

public interface GenericDrawContext {
	// for the time being we use the awt COlor object. For better compatibility with JS and Android, we might use ARGB int instead
	float getLineWidth();
	void setLineWidth(float lineWidth);
	Color getColor();
	void setColor(Color color);
	int getFontSize();
	void setFont(int size, boolean isBold, boolean isItalic);
	void drawLine(double x1, double y1, double x2, double y2);
	void drawDottedLine(double x1, double y1, double x2, double y2);
	void drawRectangle(double x1, double y1, double x2, double y2);
	void fillRectangle(double x1, double y1, double x2, double y2);
	void drawCircle(double x, double y, double d);
	void fillCircle(double x, double y, double d);
	void drawPolygon(GenericPolygon p);
	void fillPolygon(GenericPolygon p);
	void drawString(double x, double y, String s);
	void drawCenteredString(double x, double y, String s);
	void drawImage(double x, double y, Image image);
	void setClip(double x, double y, double w, double h);
	Rectangle2D getBounds(String s);
}
