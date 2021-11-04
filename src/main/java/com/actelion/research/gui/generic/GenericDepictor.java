package com.actelion.research.gui.generic;

import com.actelion.research.chem.AbstractDepictor;
import com.actelion.research.chem.StereoMolecule;

import java.awt.*;

public class GenericDepictor extends AbstractDepictor<GenericDrawContext> {
	private int		mpTextSize;
	private float	mLineWidth;

	public GenericDepictor(StereoMolecule mol) {
		super(mol);
	}

	public GenericDepictor(StereoMolecule mol, int displayMode) {
		super(mol, displayMode);
	}

	protected void init() {
		super.init();
		mLineWidth = 1.0f;
	}

	protected void drawBlackLine(DepictorLine theLine) {
		mContext.drawLine(theLine.x1, theLine.y1, theLine.x2, theLine.y2);
	}

	protected void drawDottedLine(DepictorLine theLine) {
		mContext.drawDottedLine(theLine.x1, theLine.y1, theLine.x2, theLine.y2);
	}

	protected void drawString(String theString, double x, double y) {
		mContext.drawCenteredString(x, y, theString);
	}

	protected void drawPolygon(GenericPolygon p) {
		mContext.fillPolygon(p);
	}

	protected void fillCircle(double x, double y, double d) {
		mContext.fillCircle(x, y, d);
	}

	protected double getStringWidth(String theString) {
		return mContext.getBounds(theString).getWidth();
	}

	protected void setTextSize(int theSize) {
		mpTextSize = theSize;
		if (mContext != null)
			mContext.setFont(theSize, false, false);
	}

	public int getTextSize() {
		return mpTextSize;
	}

	protected double getLineWidth() {
		return mLineWidth;
	}

	protected void setLineWidth(double lineWidth) {
		mLineWidth = (float)lineWidth;
		mContext.setLineWidth(mLineWidth);
	}

	protected void setColor(Color theColor) {
		mContext.setColor(theColor);
	}
}