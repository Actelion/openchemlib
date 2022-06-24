package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.GenericDrawContext;
import com.actelion.research.gui.generic.GenericImage;
import com.actelion.research.gui.generic.GenericPolygon;
import com.actelion.research.gui.generic.GenericRectangle;
import javafx.geometry.Bounds;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.image.Image;
import javafx.scene.paint.Color;
import javafx.scene.shape.StrokeLineCap;
import javafx.scene.text.Font;
import javafx.scene.text.FontPosture;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;

import javax.swing.*;

public class FXDrawContext implements GenericDrawContext {
	private GraphicsContext mG;
	private int mRGB;

	public FXDrawContext(GraphicsContext g) {
		mG= g;
		mG.setLineCap(StrokeLineCap.ROUND);
	}

	@Override
	public int getFontSize() {
		return (int)mG.getFont().getSize();
	}

	@Override
	public void setFont(int size, boolean isBold, boolean isItalic) {
		FontWeight fw = isBold ? FontWeight.BOLD : FontWeight.NORMAL;
		FontPosture fp = isItalic ? FontPosture.ITALIC : FontPosture.REGULAR;
		mG.setFont(Font.font(null, fw, fp, size));
	}

	@Override
	public void drawLine(double x1, double y1, double x2, double y2) {
		mG.strokeLine(x1, y1, x2, y2);
	}

	@Override
	public void drawDottedLine(double x1, double y1, double x2, double y2) {
		double[] orig = mG.getLineDashes();
		mG.setLineDashes(3*mG.getLineWidth());
		mG.strokeLine(x1, y1, x2, y2);
		mG.setLineDashes(orig);
	}

	@Override
	public void drawRectangle(double x, double y, double w, double h) {
		mG.strokeRect(x, y, w, h);
	}

	@Override
	public void fillRectangle(double x, double y, double w, double h) {
		mG.fillRect(x, y, w, h);
	}

	@Override
	public void drawCircle(double x, double y, double d) {
		mG.strokeOval(x, y, d, d);
	}

	@Override
	public void fillCircle(double x, double y, double d) {
		mG.fillOval(x, y, d, d);
	}

	@Override
	public void drawPolygon(GenericPolygon p) {
		mG.strokePolygon(p.getX(), p.getY(), p.getSize());
	}

	@Override
	public void fillPolygon(GenericPolygon p) {
		mG.fillPolygon(p.getX(), p.getY(), p.getSize());
	}

	@Override
	public float getLineWidth() {
		return (float)mG.getLineWidth();
	}

	@Override
	public void setLineWidth(float lineWidth) {
		mG.setLineWidth(lineWidth);
	}

	@Override
	public int getRGB() {
		return mRGB;
	}

	@Override
	public void setRGB(int rgb) {
		mRGB = rgb;
		float r = (float)((rgb >> 16) & 0xFF) / 255f;
		float g = (float)((rgb >> 8) & 0xFF) / 255f;
		float b = (float)(rgb & 0xFF) / 255f;
		Color color = new Color(r, g, b, 1f);
		mG.setStroke(color);
		mG.setFill(color);
	}

	@Override
	public void drawString(double x, double y, String s) {
		mG.fillText(s, x, y);
	}

	@Override
	public void drawCenteredString(double x, double y, String s) {
		double w = getBounds(s).getWidth();
		mG.fillText(s, x-w/2.0, y+mG.getFont().getSize()/3.0);
	}

	@Override
	public GenericRectangle getBounds(String s) {
		Text t = new Text(s);
		t.setFont(mG.getFont());
		Bounds b = t.getLayoutBounds();
		return new GenericRectangle(b.getMinX(), b.getMinY(), b.getWidth(), b.getHeight());
	}

	@Override
	public void drawImage(GenericImage image, double x, double y) {
		mG.drawImage((Image)image.get(), x, y);
	}

	@Override
	public void drawImage(GenericImage image, double sx, double sy, double dx, double dy, double w, double h) {
		mG.drawImage((Image)image.get(), sx, sy, w, h, dx, dy, w, h);
	}

	@Override
	public void drawImage(GenericImage image, double sx, double sy, double sw, double sh, double dx, double dy, double dw, double dh) {
		mG.drawImage((Image)image.get(), sx, sy, sw, sh, dx, dy, dw, dh);
	}

//	@Override
//	public void setClip(double x, double y, double w, double h) {
//		mParentRegion.setClip(new Rectangle(x, y, w, h));
//	}

	public boolean isDarkBackground() {
		return false;
	}

	@Override
	public int getForegroundRGB() {
		return UIManager.getColor("TextArea.foreground").getRGB();
	}

	@Override
	public int getBackgroundRGB() {
		return UIManager.getColor("TextArea.background").getRGB();
	}

	@Override
	public int getSelectionBackgroundRGB() {
		return UIManager.getColor("TextArea.selectionBackground").getRGB();
	}

	@Override
	public GenericImage createARGBImage(int width, int height) {
		return new FXImage(width, height);
	}
}
