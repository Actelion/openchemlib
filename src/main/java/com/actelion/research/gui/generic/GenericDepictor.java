package com.actelion.research.gui.generic;

import com.actelion.research.chem.AbstractDepictor;
import com.actelion.research.chem.StereoMolecule;

public class GenericDepictor extends AbstractDepictor<GenericDrawContext> {
	private int     mTextSize;
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
		mTextSize = theSize;
		if (mContext != null)
			mContext.setFont(theSize, false, false);
	}

	public int getTextSize() {
		return mTextSize;
	}

	@Override
	/**
	 * Draws a circular atom background that fades from center towards circle border
	 * @param argb if alpha < 1 then the background is mixed in accordingly
	 * @param radius <= 1.0; if null, then a default of 0.5 of the average bond length is used
	 */
	public void hiliteAtomBackgrounds(int[] atomARGB, float[] radius) {
		double avbl = getTransformation().getScaling() * getMolecule().getAverageBondLength();
		double maxRadius = (radius == null) ? 0.6 * avbl : 0.75 * avbl;
		GenericRectangle rect = simpleCalculateBounds();

		int imageX = (int)Math.floor(rect.x - maxRadius);
		int imageY = (int)Math.floor(rect.y - maxRadius);
		int imageW = (int)Math.ceil(rect.x + rect.width + maxRadius) - imageX + 1;
		int imageH = (int)Math.ceil(rect.y + rect.height + maxRadius) - imageY + 1;

		float[] imageARGB = new float[4 * imageW * imageH];
		float[] pixelARGB = new float[4];   // pixel color buffer with changing alpha values

		float[] background = getBackgroundColor().getRGBColorComponents(null);

		for (int atom=0; atom<getMolecule().getAtoms(); atom++) {
			if ((atomARGB[atom] & 0xFF000000) != 0) {
				int x0 = (int)Math.round(getAtomX(atom)) - imageX;
				int y0 = (int)Math.round(getAtomY(atom)) - imageY;
				float alpha = (float)((atomARGB[atom] & 0xFF000000) >>> 24) / 255f;
				pixelARGB[1] = (float)((atomARGB[atom] & 0x00FF0000) >> 16) / 255f;
				pixelARGB[2] = (float)((atomARGB[atom] & 0x0000FF00) >> 8) / 255f;
				pixelARGB[3] = (float)(atomARGB[atom] & 0x000000FF) / 255f;
				if (alpha != 1f) {
					pixelARGB[1] = background[0]+alpha*(pixelARGB[1]-background[0]);
					pixelARGB[2] = background[1]+alpha*(pixelARGB[2]-background[1]);
					pixelARGB[3] = background[2]+alpha*(pixelARGB[3]-background[2]);
				}
				int atomRadius = (radius == null) ? (int)maxRadius : (int)(maxRadius * radius[atom]);
				int atomSquare = atomRadius * atomRadius;

				for (int rx=0; rx<=atomRadius; rx++) {
					for (int ry=0; ry<=atomRadius; ry++) {
						float squareDistance = rx*rx + ry*ry;
						if (squareDistance <= atomSquare) {
							pixelARGB[0] = 1f - (float)Math.sqrt(squareDistance) / atomRadius;

							if (rx != 0 && ry != 0)
								mixInColor(imageARGB, 4 * ((y0 - ry) * imageW + x0 - rx), pixelARGB);

							if (ry != 0)
								mixInColor(imageARGB, 4 * ((y0 - ry) * imageW + x0 + rx), pixelARGB);

							if (rx != 0)
								mixInColor(imageARGB, 4 * ((y0 + ry) * imageW + x0 - rx), pixelARGB);

							mixInColor(imageARGB, 4 * ((y0 + ry) * imageW + x0 + rx), pixelARGB);
						}
					}
				}
			}
		}

		GenericImage image = mContext.createARGBImage(imageW, imageH);
		int index = 0;
		for (int y = 0; y<imageH; y++) {
			for (int x=0; x<imageW; x++) {
				int c = 0;
				if (imageARGB[index] != 0f) {
					if (imageARGB[index] > 1f) {
						float intensity = imageARGB[index++];
						c = 0xFF000000;
						c += ((int)(255f * imageARGB[index++] / intensity)) << 16;
						c += ((int)(255f * imageARGB[index++] / intensity)) << 8;
						c += ((int)(255f * imageARGB[index++] / intensity));
					}
					else {
						c = ((int)(255f * imageARGB[index++])) << 24;
						c += ((int)(255f * imageARGB[index++])) << 16;
						c += ((int)(255f * imageARGB[index++])) << 8;
						c += (int)(255f * imageARGB[index++]);
					}
				}
				else {
					index += 4;
				}
				image.setRGB(x, y, c);
			}
		}

		mContext.drawImage(image, imageX, imageY);
	}

	private void mixInColor(float[] imageARGB, int index, float[] pixelARGB) {
		if (imageARGB[index] == 0) {    // if there is no previous color assigned
			imageARGB[index++] = pixelARGB[0];
			imageARGB[index++] = pixelARGB[1];
			imageARGB[index++] = pixelARGB[2];
			imageARGB[index++] = pixelARGB[3];
		}
		else {  // we have to mix in the proper ratio
			float oldIntensity = imageARGB[index];
			float newIntensity = imageARGB[index] + pixelARGB[0];
			imageARGB[index] = newIntensity;
			index++;
			imageARGB[index] = (imageARGB[index] * oldIntensity + pixelARGB[1] * pixelARGB[0]) / newIntensity;
			index++;
			imageARGB[index] = (imageARGB[index] * oldIntensity + pixelARGB[2] * pixelARGB[0]) / newIntensity;
			index++;
			imageARGB[index] = (imageARGB[index] * oldIntensity + pixelARGB[3] * pixelARGB[0]) / newIntensity;
		}
	}

	protected double getLineWidth() {
		return mLineWidth;
	}

	protected void setLineWidth(double lineWidth) {
		mLineWidth = (float)lineWidth;
		mContext.setLineWidth(mLineWidth);
	}

	protected void setRGB(int rgb) {
		mContext.setRGB(rgb);
	}
}