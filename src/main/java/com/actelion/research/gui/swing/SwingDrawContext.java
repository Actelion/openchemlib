package com.actelion.research.gui.swing;

import com.actelion.research.gui.LookAndFeelHelper;
import com.actelion.research.gui.generic.GenericDrawContext;
import com.actelion.research.gui.generic.GenericImage;
import com.actelion.research.gui.generic.GenericPolygon;
import com.actelion.research.gui.generic.GenericRectangle;

import javax.swing.*;
import java.awt.*;
import java.awt.font.GlyphVector;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;

public class SwingDrawContext implements GenericDrawContext {
	private static boolean isMac = (System.getProperty("os.name").toLowerCase().indexOf("mac") >= 0);

	private Graphics2D mG;

	public SwingDrawContext(Graphics2D g) {
		mG= g;
	}

	@Override
	public int getFontSize() {
		return mG.getFont().getSize();
		}

	@Override
	public void setFont(int size, boolean isBold, boolean isItalic) {
		int style = (isBold ? Font.BOLD : 0) + (isItalic ? Font.ITALIC : 0);
		mG.setFont(mG.getFont().deriveFont(style, (float)size));
		}

	@Override
	public void drawLine(double x1, double y1, double x2, double y2) {
		// Lines on OSX are shifted left and down when drawn in Graphics2D by 0.5. We must compensate.
		if (isMac)
			mG.draw(new Line2D.Double(x1-0.5, y1-0.5, x2-0.5, y2-0.5));
		else
			mG.draw(new Line2D.Double(x1, y1, x2, y2));
		}

	@Override
	public void drawDottedLine(double x1, double y1, double x2, double y2) {
		Stroke stroke = mG.getStroke();
		float width = ((BasicStroke)stroke).getLineWidth();
		mG.setStroke(new BasicStroke(width, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, width, new float[] {3.0f*width},	0f));
		drawLine(x1, y1, x2, y2);
		mG.setStroke(stroke);
		}

	@Override
	public void drawRectangle(double x, double y, double w, double h) {
		mG.draw(new Rectangle2D.Double(x, y, w, h));
		}

	@Override
	public void fillRectangle(double x, double y, double w, double h) {
		mG.fill(new Rectangle2D.Double(x, y, w, h));
		}

	@Override
	public void drawCircle(double x, double y, double d) {
		if (isMac)
			mG.draw(new Ellipse2D.Double(x-0.5f, y-0.5f, d, d));
		else
			mG.draw(new Ellipse2D.Double(x, y, d, d));
		}

	@Override
	public void fillCircle(double x, double y, double d) {
		if (isMac)
			mG.fill(new Ellipse2D.Double(x-0.5f, y-0.5f, d, d));
		else
			mG.fill(new Ellipse2D.Double(x, y, d, d));
		}

	@Override
	public void drawPolygon(GenericPolygon p) {
		GeneralPath polygon = new GeneralPath(GeneralPath.WIND_NON_ZERO, p.getSize());

		if (isMac) {
			polygon.moveTo((float)p.getX(0)-0.5f, (float)p.getY(0)-0.5f);
			for (int i = 1; i<p.getSize(); i++)
				polygon.lineTo((float)p.getX(i)-0.5f, (float)p.getY(i)-0.5f);
			polygon.closePath();
			}
		else {
			polygon.moveTo((float)p.getX(0), (float)p.getY(0));
			for (int i = 1; i<p.getSize(); i++)
				polygon.lineTo((float)p.getX(i), (float)p.getY(i));
			polygon.closePath();
		}

		mG.draw(polygon);
		}

	@Override
	public void fillPolygon(GenericPolygon p) {
		GeneralPath polygon = new GeneralPath(GeneralPath.WIND_NON_ZERO, p.getSize());
		polygon.moveTo((float)p.getX(0), (float)p.getY(0));
		for (int i = 1; i<p.getSize(); i++)
			polygon.lineTo((float)p.getX(i), (float)p.getY(i));
		polygon.closePath();

		mG.fill(polygon);

		if (isMac) {
			polygon = new GeneralPath(GeneralPath.WIND_NON_ZERO, p.getSize());
			polygon.moveTo((float)p.getX(0)-0.5f, (float)p.getY(0)-0.5f);
			for (int i = 1; i<p.getSize(); i++)
				polygon.lineTo((float)p.getX(i)-0.5f, (float)p.getY(i)-0.5f);
			polygon.closePath();
			}

		mG.draw(polygon);
		}

	@Override
	public float getLineWidth() {
		return ((BasicStroke)mG.getStroke()).getLineWidth();
		}

	@Override
	public void setLineWidth(float lineWidth) {
		mG.setStroke(new BasicStroke(lineWidth, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
		}

	@Override
	public int getRGB() {
		return mG.getColor().getRGB();
		}

	@Override
	public void setRGB(int rgb) {
		Color color = new Color(rgb);
		mG.setColor(color);
		mG.setPaint(color);
		}

	@Override
	public void drawString(double x, double y, String s) {
		mG.drawString(s, Math.round(x), Math.round(y));
		}

	@Override
	public void drawCenteredString(double x, double y, String s) {
		GlyphVector gv = mG.getFont().createGlyphVector(mG.getFontRenderContext(), s);
		double width = gv.getLogicalBounds().getWidth();
		mG.drawGlyphVector(gv, (float)(x-width/2.0), (float)(y+(double)mG.getFont().getSize()/3.0));
		}

	@Override
	public GenericRectangle getBounds(String s) {
		Rectangle2D r = mG.getFontMetrics().getStringBounds(s, mG);
		return new GenericRectangle(r.getX(), r.getY(), r.getWidth(), r.getHeight());
		}

	@Override
	public void drawImage(GenericImage image, double x, double y) {
		mG.drawImage((Image)image.get(), Math.round((float)x), Math.round((float)y), null);
		}

	@Override
	public void drawImage(GenericImage image, double sx, double sy, double dx, double dy, double w, double h) {
		mG.drawImage((Image)image.get(),
				Math.round((float)dx), Math.round((float)dy), Math.round((float)(dx+w)), Math.round((float)(dy+h)),
				Math.round((float)sx), Math.round((float)sy), Math.round((float)(sx+w)), Math.round((float)(sy+h)), null);
		}

	@Override
	public void drawImage(GenericImage image, double sx, double sy, double sw, double sh, double dx, double dy, double dw, double dh) {
		mG.drawImage((Image)image.get(),
				Math.round((float)dx), Math.round((float)dy), Math.round((float)(dx+dw)), Math.round((float)(dy+dh)),
				Math.round((float)sx), Math.round((float)sy), Math.round((float)(sx+sw)), Math.round((float)(sy+sh)), null);
		}

//	@Override
//	public void setClip(double x, double y, double w, double h) {
//		mG.setClip(Math.round((float)x), Math.round((float)y), Math.round((float)w), Math.round((float)h));
//		}

	@Override
	public boolean isDarkBackground() {
		return LookAndFeelHelper.isDarkLookAndFeel();
		}

	@Override
	public int getForegroundRGB() {
		Color color = UIManager.getColor("TextArea.foreground");
		return color == null ? 0 : color.getRGB();
		}

	@Override
	public int getBackgroundRGB() {
		Color color = UIManager.getColor("TextArea.background");
		return color == null ? 0 : color.getRGB();
		}

	@Override
	public int getSelectionBackgroundRGB() {
		Color color = UIManager.getColor("TextArea.selectionBackground");
		return color == null ? 0 : color.getRGB();
		}

	@Override
	public GenericImage createARGBImage(int width, int height) {
		return new SwingImage(width, height);
		}
	}
