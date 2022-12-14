/*
 * Copyright (c) 1997 - 2016
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the the copyright holder nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Thomas Sander
 */

package com.actelion.research.chem;

import com.actelion.research.gui.editor.AtomQueryFeatureDialogBuilder;
import com.actelion.research.gui.generic.GenericPoint;
import com.actelion.research.gui.generic.GenericPolygon;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.util.ColorHelper;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

public abstract class AbstractDepictor<T> {
	/*
	 * If displayMode includes cDModeColorizeAtomLabels, then all atom labels are drawn
	 * in a default color, which may be overridden by explicit atom colors or atom error colors.
	 * The default color affects the atom label and potentially shown implicit hydrogens,
	 * while explicit atom colors and atom error colors affect half of the connected bonds
	 * in addition.
	 * Color values taken from jmol:
	 * http://jmol.sourceforge.net/jscolors/#Atoms%20%28%27CPK%27%20colors,%20default%20element%20colors%29
	 */
	private static final int[] ATOM_LABEL_COLOR = {     0x00000000,
		0x00FFFFFF, 0x00D9FFFF, 0x00CC80FF, 0x00C2FF00, 0x00FFB5B5, 0x00909090, 0x003050F8, 0x00FF0D0D,
		0x0090E050, 0x00B3E3F5, 0x00AB5CF2, 0x008AFF00, 0x00BFA6A6, 0x00F0C8A0, 0x00FF8000, 0x00FFFF30,
		0x001FF01F, 0x0080D1E3, 0x008F40D4, 0x003DFF00, 0x00E6E6E6, 0x00BFC2C7, 0x00A6A6AB, 0x008A99C7,
		0x009C7AC7, 0x00E06633, 0x00F090A0, 0x0050D050, 0x00C88033, 0x007D80B0, 0x00C28F8F, 0x00668F8F,
		0x00BD80E3, 0x00FFA100, 0x00A62929, 0x005CB8D1, 0x00702EB0, 0x0000FF00, 0x0094FFFF, 0x0094E0E0,
		0x0073C2C9, 0x0054B5B5, 0x003B9E9E, 0x00248F8F, 0x000A7D8C, 0x00006985, 0x00C0C0C0, 0x00FFD98F,
		0x00A67573, 0x00668080, 0x009E63B5, 0x00D47A00, 0x00940094, 0x00429EB0, 0x0057178F, 0x0000C900,
		0x0070D4FF, 0x00FFFFC7, 0x00D9FFC7, 0x00C7FFC7, 0x00A3FFC7, 0x008FFFC7, 0x0061FFC7, 0x0045FFC7,
		0x0030FFC7, 0x001FFFC7, 0x0000FF9C, 0x0000E675, 0x0000D452, 0x0000BF38, 0x0000AB24, 0x004DC2FF,
		0x004DA6FF, 0x002194D6, 0x00267DAB, 0x00266696, 0x00175487, 0x00D0D0E0, 0x00FFD123, 0x00B8B8D0,
		0x00A6544D, 0x00575961, 0x009E4FB5, 0x00AB5C00, 0x00754F45, 0x00428296, 0x00420066, 0x00007D00,
		0x0070ABFA, 0x0000BAFF, 0x0000A1FF, 0x00008FFF, 0x000080FF, 0x00006BFF, 0x00545CF2, 0x00785CE3,
		0x008A4FE3, 0x00A136D4, 0x00B31FD4, 0x00B31FBA, 0x00B30DA6, 0x00BD0D87, 0x00C70066, 0x00CC0059,
		0x00D1004F, 0x00D90045, 0x00E00038, 0x00E6002E, 0x00EB0026, 0x00000000, 0x00000000, 0x00000000,
		0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
		0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
		0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
		0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
		0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
		0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
		0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
		0x00000000, 0x00000000,
		0x00C8C8C8, 0x00145AFF, 0x0000DCDC, 0x00E60A0A, 0x00E6E600, 0x0000DCDC, 0x00E60A0A, 0x00EBEBEB,
		0x008282D2, 0x000F820F, 0x000F820F, 0x00145AFF, 0x00E6E600, 0x003232AA, 0x00DC9682, 0x00FA9600,
		0x00FA9600, 0x00B45AB4, 0x003232AA, 0x000F820F
		};

    private static final int BOND_FG_HILITE_COLOR = 0xFFFF8000;
	private static final int BOND_BG_HILITE_COLOR = 0xFF5CA0FF;

	private static final int FG_EXCLUDE_GROUP_COLOR = 0xFFA00040;
	private static final int BG_EXCLUDE_GROUP_COLOR = 0xFFFFA0FF;

	private static final int COLOR_SELECTED = Molecule.cAtomColorRed;
	private static final int COLOR_CIP_LETTER = Molecule.cAtomColorDarkRed;
	private static final int COLOR_CHIRALITY_TEXT = Molecule.cAtomColorDarkRed;
	private static final int COLOR_UNDEFINED = -1;
    private static final int COLOR_HILITE_BOND_BG = -2;
    private static final int COLOR_HILITE_BOND_FG = -3;
    private static final int COLOR_OVERRULED = -4;
    private static final int COLOR_RGB = -5;
	private static final int COLOR_CUSTOM_FOREGROUND = -6;
	private static final int COLOR_EXCLUDE_GROUP_BG = -7;
	private static final int COLOR_EXCLUDE_GROUP_FG = -8;
	private static final int COLOR_RESTORE_PREVIOUS = -9;
	private static final int COLOR_INITIALIZE = -10;

	public static final int COLOR_BLUE = 0xFF2060FF;
	public static final int COLOR_RED = 0xFFFF0000;
	public static final int COLOR_GREEN = 0xFF00FF00;
	public static final int COLOR_MAGENTA = 0xFFC000FF;
	public static final int COLOR_ORANGE = 0xFFFFA000;
	public static final int COLOR_DARK_GREEN = 0xFF008000;
	public static final int COLOR_DARK_RED = 0xFFA00000;

	public static final int cOptAvBondLen = (int)Molecule.cDefaultAVBL;
	public static final int cColorGray = 1;	// avoid the Molecule.cAtomFlagsColor range

	protected static final int cModeMaxBondLength			= 0x0FFFF;
	protected static final int cModeInflateToAVBL			= 0x30000;
	private static final int cModeChiralTextLocation		= 0xC0000;

	// options for validateView() and updateCoords()
	public static final int cModeInflateToMaxAVBL			= 0x10000;
	public static final int cModeInflateToHighResAVBL		= 0x20000;	// like cModeInflateToMaxAVBL, but avbl is multiplied by 256 before encoding
	public static final int cModeChiralTextBelowMolecule	= 0x00000;
	public static final int cModeChiralTextAboveMolecule	= 0x40000;
	public static final int cModeChiralTextOnFrameTop		= 0x80000;
	public static final int cModeChiralTextOnFrameBottom	= 0xC0000;

	public static final int	cDModeNoTabus = 0x0001;	// these are the display mode options
	public static final int	cDModeAtomNo = 0x0002;
	public static final int	cDModeBondNo = 0x0004;
	public static final int	cDModeHiliteAllQueryFeatures = 0x0008;	// Hilite also atoms and bonds with obvious query features
	public static final int	cDModeShowMapping = 0x0010;
	public static final int	cDModeSuppressChiralText = 0x0020;
	public static final int	cDModeSuppressCIPParity = 0x0040;
	public static final int	cDModeSuppressESR = 0x0080;

	private static final int cDModeShowSymmetryAny = 0x0300;
	public static final int cDModeShowSymmetrySimple = 0x0100;
    public static final int cDModeShowSymmetryStereoHeterotopicity = 0x0200;
	public static final int	cDModeNoImplicitAtomLabelColors = 0x0400;
	public static final int	cDModeNoStereoProblem = 0x0800;
	public static final int	cDModeNoColorOnESRAndCIP = 0x1000;
	public static final int cDModeNoImplicitHydrogen = 0x2000;
	public static final int cDModeDrawBondsInGray = 0x4000;

	private static final double cFactorTextSize = 0.6;
	private static final double cFactorChiralTextSize = 0.5;
	private static final double cFactorBondSpacing = 0.15;
	private static final double cFactorBondHiliteRadius = 0.38;
	private static final double cFactorExcludeGroupRadius = 0.47;
	private static final double cFactorDotDiameter = 0.12;
	private static final double cFactorQFDiameter = 0.40;
	private static final double cFactorLineWidth = 0.06;

	private boolean[]				mAtomIsConnected;
	private boolean[]				mAtomLabelDisplayed;
	private double					mpBondSpacing,mpDotDiameter,mpLineWidth,mpQFDiameter,mpBondHiliteRadius,
									mFactorTextSize,mpExcludeGroupRadius,mChiralTextSize;
	private int						mpLabelSize,mStandardForegroundColor,mDisplayMode,mCurrentColor,mPreviousColor;
	private boolean                 mIsValidatingView;
	private ArrayList<GenericRectangle> mpTabuZone;
    private ArrayList<DepictorDot>  mpDot;
	private StereoMolecule     		mMol;
	private GenericRectangle		mBoundingRect = new GenericRectangle();
	private DepictorTransformation	mTransformation;
	private GenericPoint			mChiralTextLocation;
	private int[]					mAtomColor,mAtomHiliteColor;
	private float[]                 mAtomHiliteRadius;
	private String[]				mAtomText;
	private GenericPoint[]  		mAlternativeCoords;
	private int					    mOverruleForeground,mOverruleBackground,mBondBGHiliteColor,mBondFGHiliteColor,
									mExcludeGroupFGColor,mExcludeGroupBGColor,mCustomForeground,mCustomBackground,
									mRGBColor;
	protected T                     mContext;

	public AbstractDepictor(StereoMolecule mol) {
		this(mol, 0);
		}


	public AbstractDepictor(StereoMolecule mol, int displayMode) {
		mMol = mol;
		mDisplayMode = displayMode;
		init();
		}


	public void setDisplayMode(int displayMode) {
		mDisplayMode = displayMode;
		}


	/**
	 * Defines additional atom text to be displayed in top right
	 * position of some/all atom label. If the atom is charged, then
	 * the atom text follows the charge information.
	 * @param atomText null or String array matching atom indexes (may contain null entries)
	 */
	public void setAtomText(String[] atomText) {
		mAtomText = atomText;
		}


	/**
	 * If the foreground color is set (!=0), then molecule is drawn in the foreground
	 * color except for non carbon atoms, which are drawn in atomicNo specific
	 * colors. If a background color is given (!=0), then atom coloring and highlighting
	 * is optimized to achieve good contrasts.
	 * @param foreground null (black) or color to be used for molecule drawing
	 * @param background null (white) or color on which the molecule is drawn
	 */
	@Deprecated
	public void setForegroundColor(Color foreground, Color background) {
		setForegroundColor(foreground == null ? 0 : foreground.getRGB(),
						   background == null ? 0 : background.getRGB());
		}

	/**
	 * If the foreground color is set (!=0), then molecule is drawn in the foreground
	 * color except for non carbon atoms, which are drawn in atomicNo specific
	 * colors. If a background color is given (!=0), then atom coloring and highlighting
	 * is optimized to achieve good contrasts.
	 * @param foreground 0 or ARGB values of color to be used for molecule drawing
	 * @param background 0 or ARGB values of color on which the molecule is drawn
	 */
	public void setForegroundColor(int foreground, int background) {
		mStandardForegroundColor = COLOR_CUSTOM_FOREGROUND;
		mCustomForeground = foreground;
		mCustomBackground = background;
		updateBondHiliteColor();
		}


	/**
	 * If the overrule color is set, the entire molecule is drawn in the foreground
	 * color neglecting any atom color information. The background color is used
	 * to construct a proper bond hilite color, if bond hiliting is used.
	 * @param foreground null or color to be used for molecule drawing
	 * @param background may be null
	 */
	@Deprecated
	public void setOverruleColor(Color foreground, Color background) {
		setOverruleColor(foreground == null ? 0 : foreground.getRGB(),
						 background == null ? 0 : background.getRGB());
		}

	/**
	 * If the overrule color is set (!=0), then entire molecule is drawn in the foreground
	 * color neglecting any atom color information. The background color is used
	 * to construct a proper bond hilite color, if bond hiliting is used.
	 * @param foreground 0 or ARGB values of color to be used for molecule drawing
	 * @param background 0 or ARGB values of background color
	 */
	public void setOverruleColor(int foreground, int background) {
		mOverruleForeground = foreground;
		mOverruleBackground = background;
		updateBondHiliteColor();
		}

	/**
	 * If you want the Depictor to draw an atom background with specific colors for every atom,
	 * then you need to call this method before calling paint().
	 * @param argb values with a==0 are not considered
	 * @param radius <= 1.0; if null, then a default of 0.5 of the average bond length is used
	 */
	public void setAtomHighlightColors(int[] argb, float[] radius) {
		mAtomHiliteColor = argb;
		mAtomHiliteRadius = radius;
		}


	public void setTransformation(DepictorTransformation t) {
		mTransformation = t;
		}


	/**
	 * Sets a multiplication factor to the text size of all labels. The default is 1.0.
	 * @param factor text size factor
	 */
	public void setFactorTextSize(double factor) {
		mFactorTextSize = factor;
		}


	public DepictorTransformation getTransformation() {
		return mTransformation;
		}


	public void applyTransformation(DepictorTransformation t) {
		t.applyTo(mTransformation);
		t.applyTo(mBoundingRect);
		t.applyTo(mChiralTextLocation);
		}



	/**
	 * Scales and translates the depictor's molecule's coordinates such that the molecule
	 * fits into the supplied rectangle.
	 * This method considers atom labels when generating the bounding box, which includes
	 * atom labels of drawn implicit hydrogen atoms.
	 * @param context
	 * @param viewRect
	 * @param mode is typically (cModeInflateToMaxAVBL | maximum_desired_bond_length)
	 * @return transformation that moved/scaled the original molecule into viewRect
	 */
	public DepictorTransformation updateCoords(T context, GenericRectangle viewRect, int mode) {
		validateView(context, viewRect, mode);
		if (mTransformation.isVoidTransformation()) {
			return null;
			}
		else {
			DepictorTransformation t = mTransformation;
			mTransformation.applyTo(mMol);
			mTransformation = new DepictorTransformation();
			return t;
			}
		}


	/**
	 * Scales and translates the depictor's molecule's coordinates such that the molecule
	 * fits into the supplied rectangle.
	 * This simple method creates a transformation that places all atom coordinates in
	 * the viewRect but does not know about nor considers atom labels.
	 * @param viewRect
	 * @param mode
	 * @return transformation that moved/scaled the original molecule into viewRect
	 */
	public DepictorTransformation simpleUpdateCoords(GenericRectangle viewRect, int mode) {
		simpleValidateView(viewRect, mode);
		if (mTransformation.isVoidTransformation()) {
			return null;
			}
		else {
			DepictorTransformation t = mTransformation;
			mTransformation.applyTo(mMol);
			mTransformation = new DepictorTransformation();
			return t;
			}
		}


	public StereoMolecule getMolecule() {
		return mMol;
		}


	/**
	 * A depictor maintains a DepictorTransformation object, which defines translation and scaling
	 * needed to project all of the molecule's coordinates into the viewRect.
	 * This method updates the depictor's transformation
	 * such that the molecule is centered in viewRect and reduced in size if it doesn't fit.
	 * If mode contains cModeInflateToMaxAVBL, then it is scaled to reach the desired mean bond length,
	 * or to fill viewRect, whatever is reached first.
	 * @param context
	 * @param viewRect
	 * @param mode (<0> or cModeInflateToMaxAVBL) + <desired mean bond length>
	 * @return incremental transformation being applied to depictor's current transformation 
	 */
	public DepictorTransformation validateView(T context, GenericRectangle viewRect, int mode) {
		if (mMol.getAllAtoms() == 0)
			return null;

		DepictorTransformation t1 = simpleValidateView(viewRect, mode);

		mMol.ensureHelperArrays(requiredHelperArrays());
		markIsolatedAtoms();

		mpDot.clear();
		mpTabuZone.clear();

	    mContext = context;
		calculateParameters();

		mpSetNormalLabelSize();

		mIsValidatingView = true;
		for (int i=0; i<mMol.getAllAtoms(); i++)
	    	mpDrawAtom(i, null);
		mIsValidatingView = false;

		double avbl = mTransformation.getScaling() * mMol.getAverageBondLength();
		expandBoundsByTabuZones(avbl);

		setChiralTextLocation(viewRect, avbl, mode);

		if (viewRect == null || viewRect.contains(mBoundingRect))
			return t1;

		DepictorTransformation t2 = new DepictorTransformation(mBoundingRect, viewRect, avbl, mode);
		t2.applyTo(mTransformation);
		t2.applyTo(mBoundingRect);
		t2.applyTo(mChiralTextLocation);

		if (t1 == null)
			return t2;

		t2.applyTo(t1);
		return t1;
		}


	public DepictorTransformation simpleValidateView(GenericRectangle viewRect, int mode) {
	// returns incremental transformation that moves/scales already transformed molecule into viewRect
		if (mMol.getAllAtoms() == 0)
			return null;

		mBoundingRect = simpleCalculateBounds();

		double avbl = mTransformation.getScaling() * mMol.getAverageBondLength();
		DepictorTransformation t = new DepictorTransformation(mBoundingRect, viewRect, avbl, mode);

		if (t.isVoidTransformation()) {
			t = null;
			}
		else {
			t.applyTo(mTransformation);
			t.applyTo(mBoundingRect);
			}

		setChiralTextLocation(viewRect, avbl, mode);

		return t;
		}

    // This might be overridden by subclasses (e.g. SVG Depictor)
    protected void onDrawBond(int bond, double x1, double y1, double x2, double y2)
    {
        // NOOP
    }

    // This might be overridden by subclasses (e.g. SVG Depictor)
    protected void onDrawAtom(int atom, String symbol, double x, double y)
    {
        // NOOP
    }


    public GenericRectangle simpleCalculateBounds() {
	    double minx = getAtomX(0);	// determine size of molecule
	    double maxx = getAtomX(0);
	    double miny = getAtomY(0);
	    double maxy = getAtomY(0);

		for (int i=0; i<mMol.getAllAtoms(); i++) {
			if (minx > getAtomX(i)) minx = getAtomX(i);
			if (maxx < getAtomX(i)) maxx = getAtomX(i);
			if (miny > getAtomY(i)) miny = getAtomY(i);
			if (maxy < getAtomY(i)) maxy = getAtomY(i);
			}

		return new GenericRectangle(minx, miny, maxx-minx, maxy-miny);
		}


	private void expandBoundsByTabuZones(double avbl) {
		for (int i=0; i<mpTabuZone.size(); i++)
			mBoundingRect = mBoundingRect.union(mpTabuZone.get(i));

		expandByHiliteBackgrounds(avbl);

		double border = 0.1 * avbl;
		mBoundingRect.x -= border;
		mBoundingRect.y -= border;
		mBoundingRect.width += 2.0*border;
		mBoundingRect.height += 2.0*border;
		}


	private void expandByHiliteBackgrounds(double avbl) {
		boolean[] isAtomHilited = new boolean[mMol.getAllAtoms()];
		for (int i=0; i<mMol.getAllBonds(); i++) {
			if (mMol.isBondBackgroundHilited(i)) {
				isAtomHilited[mMol.getBondAtom(0, i)] = true;
				isAtomHilited[mMol.getBondAtom(1, i)] = true;
				}
			}

		GenericRectangle rect = new GenericRectangle();
		for (int i=0; i<mMol.getAllAtoms(); i++) {
			double radius = (mMol.getAtomQueryFeatures(i) & Molecule.cAtomQFExcludeGroup) != 0 ? avbl*cFactorExcludeGroupRadius
					: isAtomHilited[i] ? avbl*cFactorBondHiliteRadius : 0;

			if (radius != 0) {
				double x = mTransformation.transformX(mMol.getAtomX(i));
				double y = mTransformation.transformY(mMol.getAtomY(i));
				rect.set(x-radius, y-radius, radius*2, radius*2);
				mBoundingRect = mBoundingRect.union(rect);
				}
			}
		}

	private void setChiralTextLocation(GenericRectangle viewRect, double avbl, int mode) {
		double spacing = avbl / 2.0;
		switch (mode & cModeChiralTextLocation) {
		case cModeChiralTextOnFrameBottom:
			if (viewRect != null) {
				mChiralTextLocation.x = viewRect.x + viewRect.width/2.0;
				mChiralTextLocation.y = viewRect.y + viewRect.height - spacing;
				break;
				}
		case cModeChiralTextBelowMolecule:
			mChiralTextLocation.x = mBoundingRect.x + mBoundingRect.width/2.0;
			mChiralTextLocation.y = mBoundingRect.y + mBoundingRect.height + spacing;
			if (viewRect != null && mChiralTextLocation.y > viewRect.y + viewRect.height - spacing)
				mChiralTextLocation.y = viewRect.y + viewRect.height - spacing;
			break;
		case cModeChiralTextOnFrameTop:
			if (viewRect != null) {
				mChiralTextLocation.x = viewRect.x + viewRect.width/2.0;
				mChiralTextLocation.y = viewRect.y + spacing;
				break;
				}
		case cModeChiralTextAboveMolecule:
			mChiralTextLocation.x = mBoundingRect.x + mBoundingRect.width/2.0;
			mChiralTextLocation.y = mBoundingRect.y - spacing;
			if (viewRect != null && mChiralTextLocation.y < viewRect.y + spacing)
				mChiralTextLocation.y = viewRect.y + spacing;
			break;
			}
		}


	/**
	 * @param atom
	 * @return atom x-coordinate in Depictor's display space
	 */
	public double getAtomX(int atom) {
		return mTransformation.transformX(mMol.getAtomX(atom));
		}


	/**
	 * @param atom
	 * @return atom y-coordinate in Depictor's display space
	 */
	public double getAtomY(int atom) {
		return mTransformation.transformY(mMol.getAtomY(atom));
		}


	/**
	 * Requires a call of updateCoords() or validateView() before calling this method.
	 * @return the bounding rectangle in device coordinates (of the moved/scaled molecule)
	 */
	public final GenericRectangle getBoundingRect() {
		return mBoundingRect;
		}


	protected void init() {
		mFactorTextSize = 1.0;
		mTransformation = new DepictorTransformation();
		mpTabuZone = new ArrayList<>();
		mpDot = new ArrayList<>();
		mAtomLabelDisplayed = new boolean[mMol.getAllAtoms()];
		mChiralTextLocation = new GenericPoint();
		mStandardForegroundColor = Molecule.cAtomColorNone;
		mCurrentColor = COLOR_UNDEFINED;
		updateBondHiliteColor();
		}


    private void updateBondHiliteColor() {
   		int background = (mOverruleBackground != 0) ? mOverruleBackground
			             : (mCustomBackground != 0) ? mCustomBackground : 0xFFFFFFFF;

		mBondBGHiliteColor = ColorHelper.intermediateColor(background, BOND_BG_HILITE_COLOR, 0.3f);
	    mBondFGHiliteColor = ColorHelper.getContrastColor(BOND_FG_HILITE_COLOR, background);
		mExcludeGroupBGColor = BG_EXCLUDE_GROUP_COLOR;
		mExcludeGroupFGColor = FG_EXCLUDE_GROUP_COLOR;
    	}


    private void calculateParameters() {
	    double averageBondLength = mTransformation.getScaling() * mMol.getAverageBondLength();
		mpLineWidth = averageBondLength * cFactorLineWidth;
		mpBondSpacing = averageBondLength * cFactorBondSpacing;
		mpBondHiliteRadius = averageBondLength * cFactorBondHiliteRadius;
		mpExcludeGroupRadius = averageBondLength * cFactorExcludeGroupRadius;
		mpLabelSize    = (int)(averageBondLength * mFactorTextSize * cFactorTextSize + 0.5);
		mpDotDiameter = averageBondLength * cFactorDotDiameter;
		mpQFDiameter = averageBondLength * cFactorQFDiameter;
		mChiralTextSize = averageBondLength * cFactorChiralTextSize + 0.5f;
		}


	public synchronized void paint(T context) {
		if (mMol.getAllAtoms() == 0)
		    return;

		mMol.ensureHelperArrays(requiredHelperArrays());

		mContext = context;
		calculateParameters();

		int[][] esrGroupMemberCount = mMol.getESRGroupMemberCounts();

		boolean explicitAtomColors = false;
		mAtomColor = new int[mMol.getAllAtoms()];
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			mAtomColor[atom] = mMol.getAtomColor(atom);
			if (mAtomColor[atom] != Molecule.cAtomColorNone)
				explicitAtomColors = true;
			if (mMol.isSelectedAtom(atom))
				mAtomColor[atom] = COLOR_SELECTED;
			if (mMol.getStereoProblem(atom) && (mDisplayMode & cDModeNoStereoProblem) == 0)
				mAtomColor[atom] = Molecule.cAtomColorMagenta;
			}

		setColor_(COLOR_INITIALIZE);	// to initialize the color tracking mechanism

		if (mAtomHiliteColor != null && (mAtomHiliteColor.length >= mMol.getAtoms()))
			hiliteAtomBackgrounds(mAtomHiliteColor, mAtomHiliteRadius);

		hiliteExcludeGroups();
		hiliteBondBackgrounds();
		indicateQueryFeatures();
		addChiralInfo();

		mpSetNormalLabelSize();
		setLineWidth(mpLineWidth);

		setColor_(mStandardForegroundColor);
		markIsolatedAtoms();

		mpDot.clear();
		mpTabuZone.clear();

		if ((mDisplayMode & cDModeNoTabus) != 0) {
			// we draw bonds first
			mpDrawAllBonds(esrGroupMemberCount);
			mpDrawAllDots();
			mpDrawBondQueryFeatures();
			}

		for (int i=0; i<mMol.getAllAtoms(); i++) {
			if (isHighlightedAtom(i)) {
				setColor_(COLOR_HILITE_BOND_FG);
	    		mpDrawAtom(i, esrGroupMemberCount);
				setColor_(mStandardForegroundColor);
				}
			else if (mAtomColor[i] != 0) {
				setColor_(mAtomColor[i]);
	    		mpDrawAtom(i, esrGroupMemberCount);
				setColor_(mStandardForegroundColor);
				}
			else if (!explicitAtomColors
				  && mMol.getMoleculeColor() != Molecule.cMoleculeColorNeutral
				  && mMol.getAtomicNo(i) != 1
				  && mMol.getAtomicNo(i) != 6
				  && ((mDisplayMode & cDModeNoImplicitAtomLabelColors) == 0)
				  && mMol.getAtomList(i) == null
				  && mMol.getAtomicNo(i) < ATOM_LABEL_COLOR.length) {
				setRGBColor(getContrastColor(ATOM_LABEL_COLOR[mMol.getAtomicNo(i)], i));
	    		mpDrawAtom(i, esrGroupMemberCount);
				setColor_(mStandardForegroundColor);
				}
			else {
	    		mpDrawAtom(i, esrGroupMemberCount);
				}
			}

		if ((mDisplayMode & cDModeNoTabus) == 0) {
			mpDrawAllDots();
			mpDrawBondQueryFeatures();
			mpDrawAllBonds(esrGroupMemberCount);
			}
		}


	public Color getBackgroundColor() {
		return new Color(getBackgroundRGB());
		}


	public int getBackgroundRGB() {
		return mOverruleBackground != 0 ? mOverruleBackground : mCustomBackground != 0 ? mCustomBackground : 0xFFFFFFFF;
		}


	private int getContrastColor(int rgb, int atom) {
		int bg = (mOverruleBackground != 0) ? mOverruleBackground
				: (mAtomHiliteColor != null && atom < mAtomHiliteColor.length && (mAtomHiliteColor[atom] & 0xFF000000) != 0) ? mAtomHiliteColor[atom]
				: (mCustomBackground != 0) ? mCustomBackground : 0xFFFFFFFF;
		return ColorHelper.getContrastColor(rgb, bg);
		}


	private boolean isHighlightedAtom(int atom) {
		if (mMol.getAllConnAtoms(atom) == 0)
			return false;

		for (int i=0; i<mMol.getAllConnAtoms(atom); i++)
			if (!mMol.isBondForegroundHilited(mMol.getConnBond(atom, i)))
				return false;

		return true;
		}


	private int requiredHelperArrays() {
	    return ((mDisplayMode & cDModeShowSymmetrySimple) != 0) ? Molecule.cHelperSymmetrySimple
	         : ((mDisplayMode & cDModeShowSymmetryStereoHeterotopicity) != 0) ? Molecule.cHelperSymmetryStereoHeterotopicity
             : Molecule.cHelperCIP;
	    }


	private void markIsolatedAtoms() {
		mAtomIsConnected = new boolean[mMol.getAllAtoms()];
		for (int bnd=0; bnd<mMol.getAllBonds(); bnd++) {
			mAtomIsConnected[mMol.getBondAtom(0,bnd)] = true;
			mAtomIsConnected[mMol.getBondAtom(1,bnd)] = true;
			}
		}


	private void addChiralInfo() {
		if ((mDisplayMode & cDModeSuppressChiralText) != 0)
			return;

		String chiralText = mMol.getChiralText();

		if (chiralText != null) {
			if (mChiralTextLocation.x == 0.0 && mChiralTextLocation.y == 0.0) {
				double avbl = mTransformation.getScaling() * mMol.getAverageBondLength();
				mBoundingRect = simpleCalculateBounds();
				expandBoundsByTabuZones(avbl);
				setChiralTextLocation(null, avbl, 0);
				}

			setTextSize((int)mChiralTextSize);
			if (mMol.getMoleculeColor() != Molecule.cMoleculeColorNeutral)
				setColor_(COLOR_CHIRALITY_TEXT);
			drawString(chiralText, mChiralTextLocation.x, mChiralTextLocation.y+0.3f*mChiralTextSize);
			}
		}


	/**
	 * May be overridden to get a more appealing highlighting than just round circles
	 * @param argb if alpha < 1 then the background is mixed in accordingly
	 * @param radius <= 1.0; if null, then a default of 0.5 of the average bond length is used
	 */
	public void hiliteAtomBackgrounds(int[] argb, float[] radius) {
		int background = (mOverruleBackground != 0) ? mOverruleBackground
				: (mCustomBackground != 0) ? mCustomBackground : 0xFFFFFFFF;

		double avbl = mTransformation.getScaling() * mMol.getAverageBondLength();
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			int alpha = (argb[atom] & 0xFF000000) >>> 24;
			if (alpha != 0) {
				int rgb = argb[atom];
				if (alpha != 255)
					rgb = ColorHelper.intermediateColor(background, argb[atom], (float)alpha/255f);
				double r = (radius == null) ? 0.5 * avbl : 0.6 * radius[atom] * avbl;
				setRGB(rgb);
	            fillCircle(getAtomX(atom)-r, getAtomY(atom)-r, 2*r);
				}
			}
		}


	private void hiliteBondBackgrounds() {
        setLineWidth(2*mpBondHiliteRadius);
        DepictorLine line = new DepictorLine();
        for (int bond=0; bond<mMol.getAllBonds(); bond++) {
			int atom1 = mMol.getBondAtom(0, bond);
			int atom2 = mMol.getBondAtom(1, bond);
        	if (mMol.isBondBackgroundHilited(bond)) {
	        	line.x1 = getAtomX(atom1);
	        	line.y1 = getAtomY(atom1);
	            line.x2 = getAtomX(atom2);
	            line.y2 = getAtomY(atom2);
				setColor_(COLOR_HILITE_BOND_BG);
	            drawBlackLine(line);
        		}
        	}
		}


	private void hiliteExcludeGroups() {
		if (mMol.isFragment()) {
			double radius = mpExcludeGroupRadius;
			setColor_(COLOR_EXCLUDE_GROUP_BG);
			for (int atom = 0; atom < mMol.getAtoms(); atom++)
				if ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFExcludeGroup) != 0)
					fillCircle(getAtomX(atom)-radius, getAtomY(atom)-radius, 2*radius);

			setLineWidth(2f*mpExcludeGroupRadius);
			DepictorLine line = new DepictorLine();
			for (int bond = 0; bond < mMol.getAllBonds(); bond++) {
				int atom1 = mMol.getBondAtom(0, bond);
				int atom2 = mMol.getBondAtom(1, bond);
				if ((mMol.getAtomQueryFeatures(atom1)
						& mMol.getAtomQueryFeatures(atom2)
						& Molecule.cAtomQFExcludeGroup) != 0) {
					line.x1 = getAtomX(atom1);
					line.y1 = getAtomY(atom1);
					line.x2 = getAtomX(atom2);
					line.y2 = getAtomY(atom2);
					drawBlackLine(line);
					}
				}
			}
		}

	private void indicateQueryFeatures() {
		if (mMol.isFragment()) {
			setColor_(Molecule.cAtomColorOrange);
			if (((mDisplayMode & cDModeHiliteAllQueryFeatures) != 0))
				for (int atom=0; atom<mMol.getAtoms(); atom++)
					if ((mMol.getAtomQueryFeatures(atom) & ~Molecule.cAtomQFExcludeGroup) != 0)
						fillCircle(getAtomX(atom)-mpQFDiameter/2,
								   getAtomY(atom)-mpQFDiameter/2,
								   mpQFDiameter);
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				if (mMol.getBondQueryFeatures(bond) != 0) {
					int atom1 = mMol.getBondAtom(0, bond);
					int atom2 = mMol.getBondAtom(1, bond);
					fillCircle((getAtomX(atom1)+getAtomX(atom2)-mpQFDiameter)/2,
							   (getAtomY(atom1)+getAtomY(atom2)-mpQFDiameter)/2,
							   mpQFDiameter);
					}
				}
			}
		}


	private void mpDrawAllBonds(int[][] esrGroupMemberCount) {
		int origColor = mStandardForegroundColor;
		int origForeground = mCustomForeground;
		if ((mDisplayMode & cDModeDrawBondsInGray) != 0) {
			mStandardForegroundColor = COLOR_CUSTOM_FOREGROUND;
			mCustomForeground = 0xFF808080;
			setColor_(cColorGray);
			}

		mAlternativeCoords = new GenericPoint[mMol.getAllAtoms()];

    		// add all double bonds first because they may set alternative coords for single bonds
		for (int i=0; i<mMol.getAllBonds(); i++)
			if (mMol.getBondType(i) == Molecule.cBondTypeDouble
			 || mMol.getBondType(i) == Molecule.cBondTypeCross
			 || mMol.getBondType(i) == Molecule.cBondTypeDelocalized)
				mpDrawBond(i);

		for (int i=0; i<mMol.getAllBonds(); i++)
			if (mMol.getBondType(i) != Molecule.cBondTypeDouble
			 && mMol.getBondType(i) != Molecule.cBondTypeCross
			 && mMol.getBondType(i) != Molecule.cBondTypeDelocalized)
				mpDrawBond(i);

		if ((mDisplayMode & cDModeSuppressCIPParity) == 0) {
			for (int i=0; i<mMol.getAllBonds(); i++) {
				if (mMol.getBondCIPParity(i) != 0) {
					String cipStr = null;

					if (mMol.getBondCIPParity(i) == Molecule.cBondCIPParityEorP
					 || mMol.getBondCIPParity(i) == Molecule.cBondCIPParityZorM) {
						if (mMol.getBondOrder(i) == 2
						 || mMol.getBondESRType(i) == Molecule.cESRTypeAbs
						 || esrGroupMemberCount == null
						 || esrGroupMemberCount[mMol.getBondESRType(i)][mMol.getBondESRGroup(i)] > 1) {
							if (mMol.getBondCIPParity(i) == Molecule.cBondCIPParityEorP)
								cipStr = (mMol.getBondOrder(i) == 2) ? "E" : mMol.isBondParityPseudo(i) ? "p" : "P";
							else
								cipStr = (mMol.getBondOrder(i) == 2) ? "Z" : mMol.isBondParityPseudo(i) ? "m" : "M";
							}
						}
					else {
						cipStr = "?";
						}

					if (cipStr != null) {
						mpSetSmallLabelSize();
						setColor_(mMol.isBondForegroundHilited(i) ? COLOR_HILITE_BOND_FG :
								(mMol.getMoleculeColor() == Molecule.cMoleculeColorNeutral || (mDisplayMode & cDModeNoColorOnESRAndCIP) != 0) ?
										mStandardForegroundColor : COLOR_CIP_LETTER);
						int atom1 = mMol.getBondAtom(0,i);
						int atom2 = mMol.getBondAtom(1,i);
						double x = (getAtomX(atom1) + getAtomX(atom2)) / 2;
						double y = (getAtomY(atom1) + getAtomY(atom2)) / 2;
						double dx = (getAtomX(atom1) - getAtomX(atom2)) / 3;
						double dy = (getAtomY(atom1) - getAtomY(atom2)) / 3;
						mpDrawString(x+dy,y-dx, cipStr,true);
						setColor_(mStandardForegroundColor);
						mpSetNormalLabelSize();
						}
					}
				}
			}

		if ((mDisplayMode & cDModeBondNo) != 0) {
			mpSetSmallLabelSize();
			setColor_(Molecule.cAtomColorDarkGreen);
			for (int i=0; i<mMol.getAllBonds(); i++) {
				int atom1 = mMol.getBondAtom(0,i);
				int atom2 = mMol.getBondAtom(1,i);
				String type = mMol.isDelocalizedBond(i) ? "d"
							: mMol.isAromaticBond(i) ? "a" : "";
				double x = (getAtomX(atom1) + getAtomX(atom2)) / 2;
				double y = (getAtomY(atom1) + getAtomY(atom2)) / 2;
				mpDrawString(x,y,type+String.valueOf(i),true);
				}
			setColor_(mStandardForegroundColor);
			mpSetNormalLabelSize();
			}

		if ((mDisplayMode & cDModeDrawBondsInGray) != 0) {
			mStandardForegroundColor = origColor;
			mCustomForeground = origForeground;
			}
		}


	private void mpDrawBond(int bnd) {
		DepictorLine theLine = new DepictorLine();
		DepictorLine aLine = new DepictorLine();
		DepictorLine bLine = new DepictorLine();
		GenericPoint piBondOffset = new GenericPoint();
		GenericPoint nextBondOffset = new GenericPoint();

		int atom1 = mMol.getBondAtom(0,bnd);
		int atom2 = mMol.getBondAtom(1,bnd);

        onDrawBond(bnd,getAtomX(atom1),getAtomY(atom1),getAtomX(atom2),getAtomY(atom2));

		// if one of the bond atoms is part of an exclude group
		if (!mMol.isSelectedAtom(atom1)
		 && !mMol.isSelectedAtom(atom2)
		 && ((mMol.getAtomQueryFeatures(atom1)
			| mMol.getAtomQueryFeatures(atom2)) & Molecule.cAtomQFExcludeGroup) != 0)
			setColor_(COLOR_EXCLUDE_GROUP_FG);

		if (mAlternativeCoords[atom1] == null) {
			theLine.x1 = getAtomX(atom1);
			theLine.y1 = getAtomY(atom1);
			}
		else {
			theLine.x1 = mAlternativeCoords[atom1].x;
			theLine.y1 = mAlternativeCoords[atom1].y;
			}

		if (mAlternativeCoords[atom2] == null) {
			theLine.x2 = getAtomX(atom2);
			theLine.y2 = getAtomY(atom2);
			}
		else {
			theLine.x2 = mAlternativeCoords[atom2].x;
			theLine.y2 = mAlternativeCoords[atom2].y;
			}

        if ((mMol.getBondQueryFeatures(bnd) & Molecule.cBondQFBridge) != 0) {
            mpHandleDottedLine(theLine, atom1, atom2);
			setColor_(COLOR_RESTORE_PREVIOUS);
            return;
            }

        int bondOrder = (mMol.getBondType(bnd) == Molecule.cBondTypeDelocalized) ? 0
                      : (mMol.getBondType(bnd) == Molecule.cBondTypeMetalLigand) ? 1
					  : mMol.getBondOrder(bnd);

		switch (bondOrder) {
		case 1:
			int bondType = mMol.getBondType(bnd);
			if ((mDisplayMode & cDModeSuppressESR) != 0
			 && (bondType == Molecule.cBondTypeUp || bondType == Molecule.cBondTypeDown)) {
				int stereoCenter = mMol.getBondAtom(0, bnd);
				int esrType = mMol.getAtomESRType(stereoCenter);
				if (esrType != Molecule.cESRTypeAbs) {
					int esrGroup = mMol.getAtomESRGroup(stereoCenter);
					int count = 0;
					for (int atom=0; atom<mMol.getAtoms(); atom++)
						if (mMol.getAtomESRType(atom) == esrType
						 && mMol.getAtomESRGroup(atom) == esrGroup)
							count++;
					if (count == 1)
						bondType = Molecule.cBondTypeSingle;
					}
				}
            switch (bondType) {
			case Molecule.cBondTypeSingle:
				mpHandleLine(theLine, atom1, atom2);
				break;
			case Molecule.cBondTypeUp:
				mpHandleWedge(theLine, atom1, atom2);
				break;
			case Molecule.cBondTypeDown:
				double xdiff = theLine.x2 - theLine.x1;
				double ydiff = theLine.y2 - theLine.y1;

                int color1,color2;
        		if (mMol.isBondForegroundHilited(mMol.getBond(atom1, atom2))) {
        			color1 = COLOR_HILITE_BOND_FG;
        			color2 = COLOR_HILITE_BOND_FG;
        			}
        		else {
        			color1 = mAtomColor[atom1];
        	        color2 = getESRColor(atom1);
        	        if (color1 == mMol.getAtomColor(atom1))	// if it is not selected or stereo error
        	            color1 = color2;
        			}

                for (int i=2; i<17; i+=2) {
					aLine.x1 = theLine.x1 + i*xdiff/17 - i*ydiff/128;
					aLine.y1 = theLine.y1 + i*ydiff/17 + i*xdiff/128;
					aLine.x2 = theLine.x1 + i*xdiff/17 + i*ydiff/128;
					aLine.y2 = theLine.y1 + i*ydiff/17 - i*xdiff/128;
					if (mpProperLine(aLine)) {
						setColor_((i<9) ? color1 : color2);
						drawBlackLine(aLine);
						setColor_(mStandardForegroundColor);
						}
					}
				break;
            case Molecule.cBondTypeMetalLigand:
                mpHandleShortDashedLine(theLine, atom1, atom2);
                break;
			}
			break;
		case 0:	// bonds defined to be aromatic
		case 2:
			if ((mAtomLabelDisplayed[atom1] || mMol.getAtomPi(atom1) == 2)
		 	 && (mAtomLabelDisplayed[atom2] || mMol.getAtomPi(atom2) == 2)
		 	 && !mMol.isRingBond(bnd)
		 	 && bondOrder == 2) {
											// double bond connecting two atoms both being either
											// a hetero atom or the central atom of an allene
				if (!mpProperLine(theLine))
					break;

				mpCalcPiBondOffset(theLine.x2 - theLine.x1,
								   theLine.y2 - theLine.y1,piBondOffset);
				double xdiff = piBondOffset.x / 2;
				double ydiff = piBondOffset.y / 2;

				aLine.x1 = theLine.x1 + xdiff;
				aLine.y1 = theLine.y1 + ydiff;
				aLine.x2 = theLine.x2 + xdiff;
				aLine.y2 = theLine.y2 + ydiff;
				bLine.x1 = theLine.x1 - xdiff;
				bLine.y1 = theLine.y1 - ydiff;
				bLine.x2 = theLine.x2 - xdiff;
				bLine.y2 = theLine.y2 - ydiff;

				if (mMol.getBondType(bnd) == Molecule.cBondTypeCross)
					mpMakeCrossBond(aLine,bLine);

				drawLine(aLine, atom1, atom2);
				if (bondOrder == 2)
					drawLine(bLine, atom1, atom2);
				else
					drawDashedLine(bLine, atom1, atom2);
				}
			else if ((mAtomLabelDisplayed[atom2] || mMol.getAtomPi(atom2) == 2)
			      && bondOrder == 2) {
											// C=X double bond with atm1 is the carbon
											// or R2C=C=CR2 with atm2 is central atom
				mpDBFromNonLabelToLabel(theLine, bnd, false);
				}
			else if ((mAtomLabelDisplayed[atom1] || mMol.getAtomPi(atom1) == 2)
		          && bondOrder == 2) {
											// C=X double bond with atm2 is the carbon
											// or R2C=C=CR2 with atm1 is central atom
				mpDBFromNonLabelToLabel(theLine, bnd, true);
				}
			else {
								// standard carbon-carbon double bond. Thus,one bond
								// connects the atom centers and the other lies aside.
				int side = mpPreferredSide(bnd);
				if (side == 0) side = 1;
				aLine.x1 = theLine.x1;
				aLine.y1 = theLine.y1;
				aLine.x2 = theLine.x2;
				aLine.y2 = theLine.y2;

				mpCalcPiBondOffset(theLine.x2 - theLine.x1,
								   theLine.y2 - theLine.y1,piBondOffset);

				if (side > 0) {
					bLine.x1 = theLine.x1 + piBondOffset.x;
					bLine.y1 = theLine.y1 + piBondOffset.y;
					bLine.x2 = theLine.x2 + piBondOffset.x;
					bLine.y2 = theLine.y2 + piBondOffset.y;

					if (mpCalcNextBondOffset(atom1, atom2, 1, nextBondOffset)
					 || (mMol.getConnAtoms(atom1) > 1)) {
						bLine.x1 += nextBondOffset.x + piBondOffset.y;
						bLine.y1 += nextBondOffset.y - piBondOffset.x;
						}

					if (mpCalcNextBondOffset(atom2, atom1, -1, nextBondOffset)
					 || (mMol.getConnAtoms(atom2) > 1)) {
						bLine.x2 += nextBondOffset.x - piBondOffset.y;
						bLine.y2 += nextBondOffset.y + piBondOffset.x;
						}
					}
				else {
					bLine.x1 = theLine.x1 - piBondOffset.x;
					bLine.y1 = theLine.y1 - piBondOffset.y;
					bLine.x2 = theLine.x2 - piBondOffset.x;
					bLine.y2 = theLine.y2 - piBondOffset.y;

					if (mpCalcNextBondOffset(atom1, atom2, -1, nextBondOffset)
					 || (mMol.getConnAtoms(atom1) > 1)) {
						bLine.x1 += nextBondOffset.x + piBondOffset.y;
						bLine.y1 += nextBondOffset.y - piBondOffset.x;
						}

					if (mpCalcNextBondOffset(atom2, atom1, 1, nextBondOffset)
					 || (mMol.getConnAtoms(atom2) > 1)) {
						bLine.x2 += nextBondOffset.x - piBondOffset.y;
						bLine.y2 += nextBondOffset.y + piBondOffset.x;
						}
					}

				if (mMol.getBondType(bnd) == Molecule.cBondTypeCross)
					mpMakeCrossBond(aLine,bLine);

				mpHandleLine(aLine, atom1, atom2);
				if (bondOrder == 2)
					mpHandleLine(bLine, atom1, atom2);
				else
					mpHandleDashedLine(bLine, atom1, atom2);
				}
			break;
		case 3:
			if (mpProperLine(theLine)) {
				drawLine(theLine, atom1, atom2);
				mpCalcPiBondOffset(theLine.x2 - theLine.x1,
								   theLine.y2 - theLine.y1, piBondOffset);
				drawOffsetLine(theLine, atom1, atom2, piBondOffset.x, piBondOffset.y, aLine);
				drawOffsetLine(theLine, atom1, atom2, -piBondOffset.x, -piBondOffset.y, aLine);
				}
			break;
		case 4:
			if (mpProperLine(theLine)) {
				mpCalcPiBondOffset(theLine.x2 - theLine.x1,
						theLine.y2 - theLine.y1, piBondOffset);
				drawOffsetLine(theLine, atom1, atom2, 1.5*piBondOffset.x, 1.5*piBondOffset.y, aLine);
				drawOffsetLine(theLine, atom1, atom2, 0.5*piBondOffset.x, 0.5*piBondOffset.y, aLine);
				drawOffsetLine(theLine, atom1, atom2, -0.5*piBondOffset.x, -0.5*piBondOffset.y, aLine);
				drawOffsetLine(theLine, atom1, atom2, -1.5*piBondOffset.x, -1.5*piBondOffset.y, aLine);
				}
			break;
		case 5:
			if (mpProperLine(theLine)) {
				drawLine(theLine, atom1, atom2);
				mpCalcPiBondOffset(theLine.x2 - theLine.x1,
						theLine.y2 - theLine.y1, piBondOffset);
				drawOffsetLine(theLine, atom1, atom2, 2*piBondOffset.x, 2*piBondOffset.y, aLine);
				drawOffsetLine(theLine, atom1, atom2, piBondOffset.x, piBondOffset.y, aLine);
				drawOffsetLine(theLine, atom1, atom2, -piBondOffset.x, -piBondOffset.y, aLine);
				drawOffsetLine(theLine, atom1, atom2, -2*piBondOffset.x, -2*piBondOffset.y, aLine);
				}
			break;
            }

		if (mCurrentColor == COLOR_EXCLUDE_GROUP_FG)
			setColor_(COLOR_RESTORE_PREVIOUS);
		}


	private void drawOffsetLine(DepictorLine theLine, int atom1, int atom2, double dx, double dy, DepictorLine aLine) {
		aLine.x1 = theLine.x1 + dx;
		aLine.y1 = theLine.y1 + dy;
		aLine.x2 = theLine.x2 + dx;
		aLine.y2 = theLine.y2 + dy;
		drawLine(aLine, atom1, atom2);
		}


	private void mpDBFromNonLabelToLabel(DepictorLine theLine, int bnd, boolean inverted) {
		DepictorLine aLine = new DepictorLine();
		DepictorLine bLine = new DepictorLine();
		GenericPoint piBondOffset = new GenericPoint();
		GenericPoint nextBondOffset = new GenericPoint();

		int atm1 = mMol.getBondAtom(0,bnd);
		int atm2 = mMol.getBondAtom(1,bnd);

		if (inverted) {
			double td = theLine.x1;
			theLine.x1 = theLine.x2;
			theLine.x2 = td;
			td = theLine.y1;
			theLine.y1 = theLine.y2;
			theLine.y2 = td;
			int ti = atm1;
			atm1 = atm2;
			atm2 = ti;
			}

		if (!mpProperLine(theLine))
			return;

		if (mMol.isRingBond(bnd)) {	// don't draw a centered double bond when
									// bond is in a ring
			aLine.x1 = theLine.x1;
			aLine.y1 = theLine.y1;
			aLine.x2 = theLine.x2;
			aLine.y2 = theLine.y2;

			int side = (inverted) ? -mpPreferredSide(bnd) : mpPreferredSide(bnd);
			if (side == 0) side = 1;

			mpCalcPiBondOffset(theLine.x2 - theLine.x1,
							   theLine.y2 - theLine.y1,piBondOffset);

			if (side > 0) {
				bLine.x1 = theLine.x1 + piBondOffset.x;
				bLine.y1 = theLine.y1 + piBondOffset.y;
				bLine.x2 = theLine.x2 + piBondOffset.x;
				bLine.y2 = theLine.y2 + piBondOffset.y;

				if (mpCalcNextBondOffset(atm1,atm2,1,nextBondOffset)
				 || (mMol.getConnAtoms(atm1) > 1)) {
					bLine.x1 += nextBondOffset.x + piBondOffset.y;
					bLine.y1 += nextBondOffset.y - piBondOffset.x;
					}
				}
			else {
				bLine.x1 = theLine.x1 - piBondOffset.x;
				bLine.y1 = theLine.y1 - piBondOffset.y;
				bLine.x2 = theLine.x2 - piBondOffset.x;
				bLine.y2 = theLine.y2 - piBondOffset.y;

				if (mpCalcNextBondOffset(atm1,atm2,-1,nextBondOffset)
				 || (mMol.getConnAtoms(atm1) > 1)) {
					bLine.x1 += nextBondOffset.x + piBondOffset.y;
					bLine.y1 += nextBondOffset.y - piBondOffset.x;
					}
				}

			if (mMol.getBondType(bnd) == Molecule.cBondTypeCross)
				mpMakeCrossBond(aLine,bLine);

			mpHandleLine(aLine,atm1,atm2);// the central line
			if (mMol.getBondType(bnd) == Molecule.cBondTypeDelocalized)
				mpHandleDashedLine(bLine,atm1,atm2);
			else
				mpHandleLine(bLine,atm1,atm2);
			}
		else {
			mpCalcPiBondOffset(theLine.x2 - theLine.x1,
							   theLine.y2 - theLine.y1,piBondOffset);
			double xdiff = piBondOffset.x / 2;
			double ydiff = piBondOffset.y / 2;

			boolean aLineIsInnerLine = false;
//			boolean bLineIsInnerLine = false;

			aLine.x1 = theLine.x1 + xdiff;
			aLine.y1 = theLine.y1 + ydiff;
			aLine.x2 = theLine.x2 + xdiff;
			aLine.y2 = theLine.y2 + ydiff;

			if (mMol.getConnAtoms(atm1) > 1) {
				if (!mpCalcNextBondOffset(atm1,atm2,1,nextBondOffset)) {
					mAlternativeCoords[atm1] = new GenericPoint(aLine.x1, aLine.y1);
//					bLineIsInnerLine = true;
					}
				else {
					aLine.x1 += nextBondOffset.x;
					aLine.y1 += nextBondOffset.y;
					if (mMol.getConnAtoms(atm1) == 2) {
						if (nextBondOffset.x != 0 || nextBondOffset.y != 0) {
							aLine.x1 += piBondOffset.y;
							aLine.y1 -= piBondOffset.x;
							}
						}
					}
				}

			bLine.x1 = theLine.x1 - xdiff;
			bLine.y1 = theLine.y1 - ydiff;
			bLine.x2 = theLine.x2 - xdiff;
			bLine.y2 = theLine.y2 - ydiff;

			if (mMol.getConnAtoms(atm1) > 1) {
				if (!mpCalcNextBondOffset(atm1,atm2,0,nextBondOffset)) {
					mAlternativeCoords[atm1] = new GenericPoint(bLine.x1, bLine.y1);
					aLineIsInnerLine = true;
					}
				else {
					bLine.x1 += nextBondOffset.x;
					bLine.y1 += nextBondOffset.y;
					if (mMol.getConnAtoms(atm1) == 2) {
						if (nextBondOffset.x != 0 || nextBondOffset.y != 0) {
							bLine.x1 += piBondOffset.y;
							bLine.y1 -= piBondOffset.x;
							}
						}
					}
				}

			if (mMol.getBondType(bnd) == Molecule.cBondTypeCross)
				mpMakeCrossBond(aLine,bLine);

			if (mMol.getBondType(bnd) == Molecule.cBondTypeDelocalized) {
				if (aLineIsInnerLine) {
					drawDashedLine(aLine,atm1,atm2);
					drawLine(bLine,atm1,atm2);
					}
				else {
					drawLine(aLine,atm1,atm2);
					drawDashedLine(bLine,atm1,atm2);
					}
				}
			else {
				drawLine(aLine,atm1,atm2);
				drawLine(bLine,atm1,atm2);
				}
			}
		}


	private void mpMakeCrossBond(DepictorLine aLine, DepictorLine bLine) {
		double temp;
		temp = aLine.x2;
		aLine.x2 = bLine.x2;
		bLine.x2 = temp;
		temp = aLine.y2;
		aLine.y2 = bLine.y2;
		bLine.y2 = temp;
		}


	private void mpCalcPiBondOffset(double dx, double dy, GenericPoint piBondOffset) {
		if (dx == 0) {
			if (dy < 0)
				piBondOffset.x =   mpBondSpacing;
			else
				piBondOffset.x = - mpBondSpacing;
			piBondOffset.y = 0;
			return;
			}

		double alpha = Math.atan(dy / dx);
		if (dx < 0)
			alpha += Math.PI;
		piBondOffset.x = - (mpBondSpacing * Math.sin(alpha));
		piBondOffset.y =   (mpBondSpacing * Math.cos(alpha));
		}


	private boolean mpProperLine(DepictorLine theLine) {
		// cuts line ends according to needs of involved atoms and returns
		// 'false' if line lies entirely in tabuZones,otherwise it returns 'true'
		boolean endsExchanged,retval;

		if (theLine.x1 == theLine.x2 && theLine.y1 == theLine.y2) {
			for (GenericRectangle tabuZone:mpTabuZone)
				if (tabuZone.contains(theLine.x1, theLine.y1))
					return false;

			return true;
			}

		GenericRectangle theFrame = mpGetFrame(theLine);

		endsExchanged = false;
		if (theLine.x1 > theLine.x2) {  // first point is the one with smaller x
			mpExchangeLineEnds(theLine);
			endsExchanged = true;
			}

		for (int i=0; i<mpTabuZone.size(); i++) {
			GenericRectangle tabuZone = mpTabuZone.get(i);
				// cannot use tabuZone.intersects(theFrame) because overlapping zero-width rects would not intersect
			if (tabuZone.x > theFrame.x + theFrame.width
			 || tabuZone.y > theFrame.y + theFrame.height
			 || theFrame.x > tabuZone.x + tabuZone.width
			 || theFrame.y > tabuZone.y + tabuZone.height)
				continue;

			if (mpInTabuZone(theLine.x1,theLine.y1,i)) {
				if (mpInTabuZone(theLine.x2,theLine.y2,i)) {
					if (endsExchanged)
						mpExchangeLineEnds(theLine);
					return false;		  // entire Line lies within tabuZone boundaries
					}
				mpShortenLine(theLine,0,i);
				retval = mpProperLine(theLine);
				if (endsExchanged)
					mpExchangeLineEnds(theLine);
				return retval;
				}
			if (mpInTabuZone(theLine.x2,theLine.y2,i)) {
				mpShortenLine(theLine,1,i);
				retval = mpProperLine(theLine);
				if (endsExchanged)
					mpExchangeLineEnds(theLine);
				return retval;
				}
			}
		if (endsExchanged)
			mpExchangeLineEnds(theLine);
		return true;
		}

	private boolean mpCalcNextBondOffset(int atm1,int atm2,int side,GenericPoint nextBondOffset) {
		final double RO_LIMIT = 2.617993878;	// right outer angle limit = 150 degrees
		final double LO_LIMIT = 3.665191429;	// left  outer angle limit = 210 degrees
		final double RI_LIMIT = 0.523598776;	// right inner angle limit =  30 degrees
		final double LI_LIMIT = 5.759586531;	// left  inner angle limit = 330 degrees

		boolean retval;
		int i,remoteAtm,bnd;
		double bondAngle,theBondAngle,testAngle;
		double angleDiff,currentAngleDiff,distance;

		retval = false;
		nextBondOffset.x = 0;		// default offset if no bond within angle limit
		nextBondOffset.y = 0;

		if (side > 0)
			angleDiff = RO_LIMIT;
		else
			angleDiff = LO_LIMIT;

		theBondAngle = mMol.getBondAngle(atm1,atm2);

		for (i=0; i<mMol.getConnAtoms(atm1); i++) {
			bnd = mMol.getConnBond(atm1,i);
			bondAngle = theBondAngle;

			if (mMol.getBondAtom(0,bnd) == atm1)
				remoteAtm = mMol.getBondAtom(1,bnd);
			else
				remoteAtm = mMol.getBondAtom(0,bnd);

			if (remoteAtm == atm2)
				continue;

			testAngle = mMol.getBondAngle(atm1,remoteAtm);
			if (bondAngle < testAngle)
				bondAngle += Math.PI*2;

			currentAngleDiff = bondAngle - testAngle;

			if (side > 0) {
				if (currentAngleDiff < Math.PI)
					retval = true;
								// a bond is leading away from double bond's right side

				if (currentAngleDiff > RO_LIMIT)
					currentAngleDiff = RO_LIMIT;

				if (currentAngleDiff < RI_LIMIT)
					currentAngleDiff = RI_LIMIT;

				if (currentAngleDiff <= angleDiff) {
					angleDiff = currentAngleDiff;
					distance = mpBondSpacing * Math.tan(angleDiff - Math.PI/2) / 2;
					nextBondOffset.x = - (distance * Math.sin(bondAngle));
					nextBondOffset.y = - (distance * Math.cos(bondAngle));
					}
				}
			else {
				if (currentAngleDiff >= Math.PI)
					retval = true;
								// a bond is leading away from double bond's left side

				if (currentAngleDiff < LO_LIMIT)
					currentAngleDiff = LO_LIMIT;

				if (currentAngleDiff > LI_LIMIT)
					currentAngleDiff = LI_LIMIT;

				if (currentAngleDiff >= angleDiff) {
					angleDiff = currentAngleDiff;
					distance = mpBondSpacing * Math.tan(4.712388981 - angleDiff) / 2;
					nextBondOffset.x = - (distance * Math.sin(bondAngle));
					nextBondOffset.y = - (distance * Math.cos(bondAngle));
					}
				}
			}
		return retval;
		}


	private void mpExchangeLineEnds(DepictorLine theLine) {
		double temp;
		temp = theLine.x1;
		theLine.x1 = theLine.x2;
		theLine.x2 = temp;
		temp = theLine.y1;
		theLine.y1 = theLine.y2;
		theLine.y2 = temp;
		}


	private void mpHandleLine(DepictorLine theLine,int atm1,int atm2) {
		if (mpProperLine(theLine))
			drawLine(theLine,atm1,atm2);
		}


	private void mpHandleDashedLine(DepictorLine theLine,int atm1,int atm2) {
		if (mpProperLine(theLine))
			drawDashedLine(theLine,atm1,atm2);
		}


    private void mpHandleShortDashedLine(DepictorLine theLine,int atm1,int atm2) {
        if (mpProperLine(theLine))
            drawShortDashedLine(theLine, atm1, atm2);
        }


	private void mpHandleDottedLine(DepictorLine theLine,int atm1,int atm2) {
		if (mpProperLine(theLine))
			drawDottedLine(theLine);
		}


	private void mpHandleWedge(DepictorLine origWedge,int atm1,int atm2) {
		DepictorLine theWedge = new DepictorLine();

		if (origWedge.x1 == origWedge.x2
		 && origWedge.y1 == origWedge.y2)
			return;

		theWedge.x1 = origWedge.x1;	// use copy of data for recursive processing
		theWedge.y1 = origWedge.y1;
		theWedge.x2 = origWedge.x2;
		theWedge.y2 = origWedge.y2;

		GenericRectangle theFrame = mpGetFrame(theWedge);

		for (int i=0; i<mpTabuZone.size(); i++) {
			GenericRectangle tabuZone = mpTabuZone.get(i);
				// cannot use tabuZone.intersects(theFrame) because overlapping zero-width rects would not intersect
			if (tabuZone.x > theFrame.x + theFrame.width
			 || tabuZone.y > theFrame.y + theFrame.height
			 || theFrame.x > tabuZone.x + tabuZone.width
			 || theFrame.y > tabuZone.y + tabuZone.height)
				continue;

			if (mpInTabuZone(theWedge.x1,theWedge.y1,i)) {
				if (mpInTabuZone(theWedge.x2,theWedge.y2,i))
					return;		// entire Wedge lies within tabuZone boundaries
				mpShortenLine(theWedge,0,i);
				mpHandleWedge(theWedge,atm1,atm2);
				return;
				}
			if (mpInTabuZone(theWedge.x2,theWedge.y2,i)) {
				mpShortenLine(theWedge,1,i);
				mpHandleWedge(theWedge,atm1,atm2);
				return;
				}
			}

		drawWedge(theWedge,atm1,atm2);
		}


	private GenericRectangle mpGetFrame(DepictorLine theLine) {
		GenericRectangle theFrame = new GenericRectangle();
		if (theLine.x1 <= theLine.x2) {
			theFrame.x = theLine.x1;
			theFrame.width = theLine.x2 - theLine.x1;
			}
		else {
			theFrame.x  = theLine.x2;
			theFrame.width = theLine.x1 - theLine.x2;
			}
		if (theLine.y1 <= theLine.y2) {
			theFrame.y	 = theLine.y1;
			theFrame.height = theLine.y2 - theLine.y1;
			}
		else {
			theFrame.y	 = theLine.y2;
			theFrame.height = theLine.y1 - theLine.y2;
			}
		return theFrame;
		}


	private boolean mpInTabuZone(double x, double y, int tabuZoneNo) {
		if ((mDisplayMode & cDModeNoTabus) != 0)
			return false;

		GenericRectangle tabuZone = mpTabuZone.get(tabuZoneNo);

			// cannot use tabuZone.contains() because points on edge would be considered to be within the react
		return (x > tabuZone.x
			 && x < tabuZone.x + tabuZone.width
			 && y > tabuZone.y
			 && y < tabuZone.y + tabuZone.height);
		}


	private void mpShortenLine(DepictorLine theLine,int pointNo,int tabuZoneNo) {
		double x1,y1,x2,y2,dx,dy,tabuX,tabuY,sx,sy;

		if (pointNo == 0) {
			x1 = theLine.x1;
			y1 = theLine.y1;
			x2 = theLine.x2;
			y2 = theLine.y2;
			}
		else {
			x1 = theLine.x2;
			y1 = theLine.y2;
			x2 = theLine.x1;
			y2 = theLine.y1;
			}

		GenericRectangle tabuZone = mpTabuZone.get(tabuZoneNo);
		tabuX = (x2 > x1) ? tabuZone.x+tabuZone.width : tabuZone.x;
		tabuY = (y2 > y1) ? tabuZone.y+tabuZone.height : tabuZone.y;

		dx = x2 - x1;
		dy = y2 - y1;
		if (Math.abs(dx) > Math.abs(dy)) {
			if (y1 == y2) {
				sx = tabuX;
				sy = y1;
				}
			else {
				sx = x1 + dx*(tabuY-y1)/dy;
				if ((x2>x1)==(tabuX>sx)) {
					sy = tabuY;
					}
				else {
					sx = tabuX;
					sy = y1 + dy*(tabuX-x1)/dx;
					}
				}
			}
		else {
			if (x1 == x2) {
				sx = x1;
				sy = tabuY;
				}
			else {
				sy = y1 + dy*(tabuX-x1)/dx;
				if ((y2>y1)==(tabuY>sy)) {
					sx = tabuX;
					}
				else {
					sx = x1 + dx*(tabuY-y1)/dy;
					sy = tabuY;
					}
				}
			}

		if (pointNo == 0) {
			theLine.x1 = sx;
			theLine.y1 = sy;
			}
		else {
			theLine.x2  = sx;
			theLine.y2 = sy;
			}
		}


	private int mpPreferredSide(int bnd) {
		boolean[] isAromatic = new boolean[ExtendedMolecule.cMaxConnAtoms];
		boolean[] isInRing = new boolean[ExtendedMolecule.cMaxConnAtoms];
		double[] angle = new double[ExtendedMolecule.cMaxConnAtoms];
		double[] bondAngle = new double[2];

		int angles = 0;
		for (int i=0; i<2; i++) {
			int atm = mMol.getBondAtom(i,bnd);

			for (int j=0; j<mMol.getConnAtoms(atm); j++) {
				int connBond = mMol.getConnBond(atm,j);
				if (connBond == bnd)
					continue;

				if (angles == 4)
					return 0;

				isAromatic[angles] = mMol.isAromaticBond(connBond);
				isInRing[angles] = mMol.isRingBond(connBond);
				angle[angles++] = mMol.getBondAngle(atm, mMol.getConnAtom(atm,j));
				}
			}

		boolean changed;
		bondAngle[0] = mMol.getBondAngle(mMol.getBondAtom(0,bnd),mMol.getBondAtom(1,bnd));
		if (bondAngle[0] < 0) {
			bondAngle[1] = bondAngle[0] + Math.PI;
			changed = false;
			}
		else {
			bondAngle[1] = bondAngle[0];
			bondAngle[0] = bondAngle[1] - Math.PI;
			changed = true;
			}

		int side = 0;
		for (int i=0; i<angles; i++) {
			int value;
			if (isAromatic[i])
				value = 20;
			else if (isInRing[i])
				value = 17;
			else
				value = 16;

			if ((angle[i] > bondAngle[0]) && (angle[i] < bondAngle[1]))
				side -= value;
			else
				side += value;
			}

		return (changed) ? -side : side;
		}


	private void mpDrawAtom(int atom, int[][] esrGroupMemberCount) {
		double chax,chay,xdiff,ydiff,x,y;

        if (!mIsValidatingView)
            onDrawAtom(atom,mMol.getAtomLabel(atom), getAtomX(atom), getAtomY(atom));

		String propStr = null;
		if (mMol.getAtomCharge(atom) != 0) {
			String valStr = (Math.abs(mMol.getAtomCharge(atom)) == 1) ? ""
			                    : String.valueOf(Math.abs(mMol.getAtomCharge(atom)));
			propStr = (mMol.getAtomCharge(atom) < 0) ? valStr + "-" : valStr + "+";
			}
		if (mAtomText != null && (atom < mAtomText.length) && mAtomText[atom] != null && mAtomText[atom].length() > 0)
			propStr = append(propStr, mAtomText[atom]);

		String isoStr = null;
		long queryFeatures = mMol.getAtomQueryFeatures(atom);
		if (queryFeatures != 0) {
			if ((queryFeatures & Molecule.cAtomQFIsStereo) != 0)
				isoStr = append(isoStr, "*");
			if ((queryFeatures & Molecule.cAtomQFIsNotStereo) != 0)
				isoStr = append(isoStr, "!*");
			if ((queryFeatures & Molecule.cAtomQFHeteroAromatic) != 0)
				isoStr = append(isoStr, "ha");
			else if ((queryFeatures & Molecule.cAtomQFAromatic) != 0)
				isoStr = append(isoStr, "a");
			else if ((queryFeatures & Molecule.cAtomQFNotAromatic) != 0)
				isoStr = append(isoStr, "!a");
			if ((queryFeatures & Molecule.cAtomQFMoreNeighbours) != 0)
				isoStr = append(isoStr, "s");
            if ((queryFeatures & Molecule.cAtomQFHydrogen) != 0) {
                long hydrogens = (queryFeatures & Molecule.cAtomQFHydrogen);
    			if (hydrogens == Molecule.cAtomQFNot1Hydrogen+Molecule.cAtomQFNot2Hydrogen+Molecule.cAtomQFNot3Hydrogen)
    				isoStr = append(isoStr, "h0");
    			else if (hydrogens == Molecule.cAtomQFNot0Hydrogen+Molecule.cAtomQFNot2Hydrogen+Molecule.cAtomQFNot3Hydrogen)
    				isoStr = append(isoStr, "h1");
                else if (hydrogens == Molecule.cAtomQFNot0Hydrogen+Molecule.cAtomQFNot1Hydrogen+Molecule.cAtomQFNot3Hydrogen)
                    isoStr = append(isoStr, "h2");
    			else if (hydrogens == Molecule.cAtomQFNot0Hydrogen)
    				isoStr = append(isoStr, "h>0");
    			else if (hydrogens == Molecule.cAtomQFNot0Hydrogen+Molecule.cAtomQFNot1Hydrogen)
    				isoStr = append(isoStr, "h>1");
			    else if (hydrogens == Molecule.cAtomQFNot0Hydrogen+Molecule.cAtomQFNot1Hydrogen+Molecule.cAtomQFNot2Hydrogen)
				    isoStr = append(isoStr, "h>2");
                else if (hydrogens == Molecule.cAtomQFNot3Hydrogen)
                    isoStr = append(isoStr, "h<3");
                else if (hydrogens == Molecule.cAtomQFNot2Hydrogen+Molecule.cAtomQFNot3Hydrogen)
                    isoStr = append(isoStr, "h<2");
                }
            if ((queryFeatures & Molecule.cAtomQFCharge) != 0) {
                long charge = (queryFeatures & Molecule.cAtomQFCharge);
    			if (charge == Molecule.cAtomQFNotChargePos+Molecule.cAtomQFNotChargeNeg)
    				isoStr = append(isoStr, "c0");
    			else if (charge == Molecule.cAtomQFNotCharge0+Molecule.cAtomQFNotChargeNeg)
    				isoStr = append(isoStr, "c+");
    			else if (charge == Molecule.cAtomQFNotCharge0+Molecule.cAtomQFNotChargePos)
    				isoStr = append(isoStr, "c-");
                }
            if ((queryFeatures & Molecule.cAtomQFPiElectrons) != 0) {
                long piElectrons = (queryFeatures & Molecule.cAtomQFPiElectrons);
                if (piElectrons == Molecule.cAtomQFNot1PiElectron+Molecule.cAtomQFNot2PiElectrons)
                    isoStr = append(isoStr, "pi0");
                else if (piElectrons == Molecule.cAtomQFNot0PiElectrons+Molecule.cAtomQFNot2PiElectrons)
                    isoStr = append(isoStr, "pi1");
                else if (piElectrons == Molecule.cAtomQFNot0PiElectrons+Molecule.cAtomQFNot1PiElectron)
                    isoStr = append(isoStr, "pi2");
                else if (piElectrons == Molecule.cAtomQFNot0PiElectrons)
                    isoStr = append(isoStr, "pi>0");
                }
            if ((queryFeatures & Molecule.cAtomQFNeighbours) != 0) {
                long neighbours = (queryFeatures & Molecule.cAtomQFNeighbours);
                if (neighbours == (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot1Neighbour))
                    isoStr = append(isoStr, "n1");
                else if (neighbours == (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot2Neighbours))
                    isoStr = append(isoStr, "n2");
                else if (neighbours == (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot3Neighbours))
                    isoStr = append(isoStr, "n3");
                else if (neighbours == Molecule.cAtomQFNot3Neighbours+Molecule.cAtomQFNot4Neighbours)
                    isoStr = append(isoStr, "n<3");
                else if (neighbours == Molecule.cAtomQFNot4Neighbours)
                    isoStr = append(isoStr, "n<4");
                else if (neighbours == Molecule.cAtomQFNot0Neighbours+Molecule.cAtomQFNot1Neighbour)
                    isoStr = append(isoStr, "n>1");
                else if (neighbours == Molecule.cAtomQFNot0Neighbours+Molecule.cAtomQFNot1Neighbour+Molecule.cAtomQFNot2Neighbours)
                    isoStr = append(isoStr, "n>2");
                else if (neighbours == (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot4Neighbours))
                    isoStr = append(isoStr, "n>3");
                }
			if ((queryFeatures & Molecule.cAtomQFENeighbours) != 0) {
				long eNegNeighbours = (queryFeatures & Molecule.cAtomQFENeighbours);
				if (eNegNeighbours == (Molecule.cAtomQFENeighbours & ~Molecule.cAtomQFNot0ENeighbours))
					isoStr = append(isoStr, "e0");
				else if (eNegNeighbours == (Molecule.cAtomQFENeighbours & ~Molecule.cAtomQFNot1ENeighbour))
					isoStr = append(isoStr, "e1");
				else if (eNegNeighbours == (Molecule.cAtomQFENeighbours & ~Molecule.cAtomQFNot2ENeighbours))
					isoStr = append(isoStr, "e2");
				else if (eNegNeighbours == (Molecule.cAtomQFENeighbours & ~Molecule.cAtomQFNot3ENeighbours))
					isoStr = append(isoStr, "e3");
				else if (eNegNeighbours == (Molecule.cAtomQFNot2ENeighbours | Molecule.cAtomQFNot3ENeighbours | Molecule.cAtomQFNot4ENeighbours))
					isoStr = append(isoStr, "e<2");
				else if (eNegNeighbours == (Molecule.cAtomQFNot3ENeighbours | Molecule.cAtomQFNot4ENeighbours))
					isoStr = append(isoStr, "e<3");
				else if (eNegNeighbours == Molecule.cAtomQFNot4ENeighbours)
					isoStr = append(isoStr, "e<4");
				else if (eNegNeighbours == Molecule.cAtomQFNot0ENeighbours)
					isoStr = append(isoStr, "e>0");
				else if (eNegNeighbours == (Molecule.cAtomQFNot0ENeighbours | Molecule.cAtomQFNot1ENeighbour))
					isoStr = append(isoStr, "e>1");
				else if (eNegNeighbours == (Molecule.cAtomQFNot0ENeighbours | Molecule.cAtomQFNot1ENeighbour | Molecule.cAtomQFNot2ENeighbours))
					isoStr = append(isoStr, "e>2");
				else if (eNegNeighbours == (Molecule.cAtomQFENeighbours & ~Molecule.cAtomQFNot4ENeighbours))
					isoStr = append(isoStr, "e>3");
				else if (eNegNeighbours == (Molecule.cAtomQFNot0ENeighbours | Molecule.cAtomQFNot3ENeighbours | Molecule.cAtomQFNot3ENeighbours))
					isoStr = append(isoStr, "e1-2");
				else if (eNegNeighbours == (Molecule.cAtomQFNot0ENeighbours | Molecule.cAtomQFNot4ENeighbours))
					isoStr = append(isoStr, "e1-3");
				else if (eNegNeighbours == (Molecule.cAtomQFNot0ENeighbours | Molecule.cAtomQFNot1ENeighbour | Molecule.cAtomQFNot4ENeighbours))
					isoStr = append(isoStr, "e2-3");
				}
            if ((queryFeatures & Molecule.cAtomQFRingState) != 0) {
                long ringBonds = (queryFeatures & Molecule.cAtomQFRingState);
                if (ringBonds == Molecule.cAtomQFNot2RingBonds+Molecule.cAtomQFNot3RingBonds+Molecule.cAtomQFNot4RingBonds)
                    isoStr = append(isoStr, "!r");
                else if (ringBonds == Molecule.cAtomQFNotChain)
                    isoStr = append(isoStr, "r");
                else if (ringBonds == Molecule.cAtomQFNot3RingBonds+Molecule.cAtomQFNot4RingBonds)
	                isoStr = append(isoStr, "rb<3");
                else if (ringBonds == Molecule.cAtomQFNotChain+Molecule.cAtomQFNot3RingBonds+Molecule.cAtomQFNot4RingBonds)
                    isoStr = append(isoStr, "rb2");
                else if (ringBonds == Molecule.cAtomQFNotChain+Molecule.cAtomQFNot2RingBonds+Molecule.cAtomQFNot4RingBonds)
                    isoStr = append(isoStr, "rb3");
                else if (ringBonds == Molecule.cAtomQFNotChain+Molecule.cAtomQFNot2RingBonds+Molecule.cAtomQFNot3RingBonds)
                    isoStr = append(isoStr, "rb4");
                }
            if ((queryFeatures & Molecule.cAtomQFSmallRingSize) != 0) {
                isoStr = append(isoStr, "r"+((queryFeatures & Molecule.cAtomQFSmallRingSize)>>Molecule.cAtomQFSmallRingSizeShift));
                }
			if ((queryFeatures & Molecule.cAtomQFNewRingSize) != 0) {
				isoStr = append(isoStr, createRingSizeText(queryFeatures));
				}
            if ((queryFeatures & Molecule.cAtomQFFlatNitrogen) != 0) {
                isoStr = append(isoStr, "f");
                }
			}
		if (mMol.getAtomMass(atom) != 0) {
			isoStr = append(isoStr, String.valueOf(mMol.getAtomMass(atom)));
			}

		int unpairedElectrons = 0;
		if (mMol.getAtomRadical(atom) != 0) {
			switch (mMol.getAtomRadical(atom)) {
			case Molecule.cAtomRadicalStateS:
				propStr = append(propStr, "|");
				break;
			case Molecule.cAtomRadicalStateD:
				unpairedElectrons = 1;
				break;
			case Molecule.cAtomRadicalStateT:
				unpairedElectrons = 2;
				break;
				}
			}

		String cipStr = null;
		if ((mDisplayMode & cDModeSuppressCIPParity) == 0) {
			if (mMol.isAtomConfigurationUnknown(atom))
				cipStr = "?";
			else if (mMol.getAtomCIPParity(atom) != 0) {
				if (mMol.getAtomESRType(atom) == Molecule.cESRTypeAbs
				 || esrGroupMemberCount == null
				 || esrGroupMemberCount[mMol.getAtomESRType(atom)][mMol.getAtomESRGroup(atom)] > 1) {
					if (mMol.getConnAtoms(atom) == 2) {
						switch (mMol.getAtomCIPParity(atom)) {
						case Molecule.cAtomCIPParitySorP:
							cipStr = mMol.isAtomParityPseudo(atom) ? "p" : "P";
							break;
						case Molecule.cAtomCIPParityRorM:
							cipStr = mMol.isAtomParityPseudo(atom) ? "m" : "M";
							break;
						default:
							cipStr = "*";
							break;
							}
						}
					else {
						switch (mMol.getAtomCIPParity(atom)) {
						case Molecule.cAtomCIPParityRorM:
							cipStr = mMol.isAtomParityPseudo(atom) ? "r" : "R";
							break;
						case Molecule.cAtomCIPParitySorP:
							cipStr = mMol.isAtomParityPseudo(atom) ? "s" : "S";
							break;
						default:
							cipStr = "*";
							break;
							}
						}
					}
				}
			}
        if ((mDisplayMode & cDModeShowSymmetryAny) != 0)
            cipStr = append(cipStr, String.valueOf(mMol.getSymmetryRank(atom)));

		String mapStr = null;
		if ((mDisplayMode & cDModeShowMapping) != 0
		 && mMol.getAtomMapNo(atom) != 0)
			mapStr = "" + mMol.getAtomMapNo(atom);

        String esrStr = null;
        if (mMol.getStereoBond(atom) != -1) {
        	int esrInfo = getESRTypeToDisplayAt(atom);
            if (esrInfo != -1)
            	esrStr = (esrInfo == Molecule.cESRTypeAbs) ? "abs"
            		   : (((esrInfo & 0x00FF) == Molecule.cESRTypeAnd) ? "&" : "or") + (1 + (esrInfo >> 8));
            }

		int hydrogensToAdd = 0;
		if ((mDisplayMode & cDModeNoImplicitHydrogen) == 0) {
			if (mMol.isFragment()) {
				if ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFNoMoreNeighbours) != 0)
					hydrogensToAdd = mMol.getImplicitHydrogens(atom);
				}
			else {
				if (mMol.getAtomicNo(atom) != 6
				 || mMol.getAtomMass(atom) != 0
				 || !mAtomIsConnected[atom]
				 || mMol.getAtomRadical(atom) != 0)
					hydrogensToAdd = mMol.getImplicitHydrogens(atom);
				}
			}

		boolean largeIsoString = false;
		String atomStr = mMol.getAtomCustomLabel(atom);
		if (atomStr != null && atomStr.startsWith("]")) {
			isoStr = append(atomStr.substring(1), isoStr);
			atomStr = null;
			largeIsoString = true;
			}
		if (atomStr != null) {
		    hydrogensToAdd = 0;
		    }
		else if (mMol.getAtomList(atom) != null) {
			String atmStart = ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0) ?
								"[!" : "[";
			atomStr = atmStart+mMol.getAtomListString(atom)+"]";
			if (atomStr.length() > 5)
				atomStr = atmStart+mMol.getAtomList(atom).length+"]";
			if ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFNoMoreNeighbours) != 0)
				hydrogensToAdd = -1;
			}
		else if ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0) {
			atomStr = "?";
			if ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFNoMoreNeighbours) != 0)
				hydrogensToAdd = -1;
			}
		else if (mMol.getAtomicNo(atom) != 6
		 || propStr != null
		 || isoStr != null
		 || (hydrogensToAdd > 0)
		 || !mAtomIsConnected[atom])
			atomStr = mMol.getAtomLabel(atom);

		double labelWidth = 0.0;

		if (!mMol.isSelectedAtom(atom) & (mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFExcludeGroup) != 0)
			setColor_(COLOR_EXCLUDE_GROUP_FG);

		if (atomStr != null) {
			labelWidth = getStringWidth(atomStr);
			mpDrawString(getAtomX(atom),getAtomY(atom),atomStr,true);
			mAtomLabelDisplayed[atom] = true;
			}
		else if (mpAlleneCenter(atom))
			mpDrawDot(getAtomX(atom),getAtomY(atom),atom);

		if (propStr != null) {
			mpSetSmallLabelSize();
			x = getAtomX(atom) + ((labelWidth + getStringWidth(propStr)) / 2.0 + 1);
			y = getAtomY(atom) - ((getTextSize()*4-4)/8);
			mpDrawString(x,y,propStr,true);
			mpSetNormalLabelSize();
			}

        if ((mDisplayMode & cDModeAtomNo) != 0)
            isoStr = String.valueOf(atom);

		if (isoStr != null) {
			if (largeIsoString)
				mpSetReducedLabelSize();
			else
				mpSetSmallLabelSize();
			x = getAtomX(atom) - ((labelWidth + getStringWidth(isoStr)) / 2.0f);
			y = getAtomY(atom) - ((getTextSize()*4-4)/8);
			mpDrawString(x,y,isoStr,true);
			mpSetNormalLabelSize();
			}

		if (cipStr != null) {
			mpSetSmallLabelSize();
			x = getAtomX(atom) - ((labelWidth + getStringWidth(cipStr)) / 2.0f);
			y = getAtomY(atom) + ((getTextSize()*4+4)/8);
			int theColor = mCurrentColor;
			if (mMol.getMoleculeColor() != Molecule.cMoleculeColorNeutral && (mDisplayMode & cDModeNoColorOnESRAndCIP) == 0)
				setColor_(COLOR_CIP_LETTER);
			mpDrawString(x,y,cipStr,false);
			setColor_(theColor);
			mpSetNormalLabelSize();
			}

		if (mapStr != null) {
			mpSetSmallLabelSize();
			x = getAtomX(atom) + ((labelWidth + getStringWidth(mapStr)) / 2.0 + 1);
			y = getAtomY(atom) + ((getTextSize()*4+4)/8);
			int theColor = mCurrentColor;
			setColor_(mMol.isAutoMappedAtom(atom) ? Molecule.cAtomColorDarkGreen : Molecule.cAtomColorDarkRed);
			mpDrawString(x,y,mapStr,true);
			setColor_(theColor);
			mpSetNormalLabelSize();
			}

        if (esrStr != null) {
	        double angle = mpGetFreeSpaceAngle(atom);
			mpSetSmallLabelSize();
            x = getAtomX(atom) + 0.7*getTextSize()*Math.sin(angle);
            y = getAtomY(atom) + 0.7*getTextSize()*Math.cos(angle);
            int theColor = mCurrentColor;
			if (!mIsValidatingView && mMol.getMoleculeColor() != Molecule.cMoleculeColorNeutral)
	            setColor_(getESRColor(atom));
            mpDrawString(x,y,esrStr,false);
            setColor_(theColor);
			mpSetNormalLabelSize();
            }

        if (hydrogensToAdd == 0 && unpairedElectrons == 0) {
			if (mCurrentColor == COLOR_EXCLUDE_GROUP_FG)
				setColor_(COLOR_RESTORE_PREVIOUS);
			return;
			}

		double hindrance[] = new double[4];
		for (int i=0; i<mMol.getAllConnAtomsPlusMetalBonds(atom); i++) {
			int bnd = mMol.getConnBond(atom,i);
			for (int j=0; j<2; j++) {
				if (mMol.getBondAtom(j,bnd) == atom) {
					double theAngle = mMol.getBondAngle(mMol.getBondAtom(j,bnd),mMol.getBondAtom(1-j,bnd));
					if (theAngle < -Math.PI/2) {
						hindrance[0] -= (theAngle + Math.PI/2);
						hindrance[3] += (theAngle + Math.PI);
						}
					else if (theAngle < 0) {
						hindrance[2] += (theAngle + Math.PI/2);
						hindrance[3] -= theAngle;
						}
					else if (theAngle < Math.PI/2) {
						hindrance[1] += theAngle;
						hindrance[2] += (Math.PI/2 - theAngle);
						}
					else {
						hindrance[0] += (theAngle - Math.PI/2);
						hindrance[1] += (Math.PI - theAngle);
						}
					}
				}
			}

		if (mMol.getConnAtoms(atom) == 0) {
			if (mMol.isElectronegative(atom))	// single atoms:
				hindrance[3] -= 0.2;	// position depends on electronegativity
			else
				hindrance[1] -= 0.2;
			}
		else
			hindrance[1] -= 0.1;		// connected atoms: prefer right position slightly

		if (propStr != null || mapStr != null) hindrance[1] += 10;
								// make sure not to add H's to the right of the
								// hetero atom where they would mess up the properties
		if (isoStr != null || cipStr != null) hindrance[3] += 10;
								// make sure not to add H's to the left of the
								// hetero atom where they would mess up the atom mass

		String hNoStr = "";
		if (hydrogensToAdd != 0) {
			double hydrogenWidth = getStringWidth("H");
			double hNoWidth = 0.0;
			double hHeight = getTextSize();
			if (hydrogensToAdd == -1) {
				hNoStr = "n";
				mpSetSmallLabelSize();
				hNoWidth = getStringWidth(hNoStr);
				}
			else if (hydrogensToAdd > 1) {
				hNoStr = String.valueOf(hydrogensToAdd);
				mpSetSmallLabelSize();
				hNoWidth = getStringWidth(hNoStr);
				}

			if (hindrance[1] < 0.6 || hindrance[3] < 0.6) {
								// hindrance on the left or right are reasonably small
				chay = getAtomY(atom);
				if (hindrance[1] <= hindrance[3]) {
					hindrance[1] += 10;
					chax = getAtomX(atom) + ((labelWidth + hydrogenWidth) / 2.0f);
					}
				else {
					hindrance[3] += 10;
					chax = getAtomX(atom) - ((labelWidth + hydrogenWidth) / 2.0f) - hNoWidth;
					}
				}
			else {
				chax = getAtomX(atom);
				if (hindrance[0] < hindrance[2]) {
					hindrance[0] += 10;
					chay = getAtomY(atom) - hHeight;
					}
				else {
					hindrance[2] += 10;
					chay = getAtomY(atom) + hHeight;
					}
				}

			if (hNoWidth > 0) {
				x = chax + ((hydrogenWidth + hNoWidth) / 2.0f);
				y = chay + ((getTextSize()*4+4)/8);
				mpDrawString(x,y,hNoStr,true);
				mpSetNormalLabelSize();
				}

			mpDrawString(chax,chay,"H",true);
			}

		int bestSide = 0;
		if (unpairedElectrons != 0) {
			double minHindrance = 50.0;
			double counterHindrance = 0.0;
			for (int i=0; i<4; i++) {
				int counterSide = (i > 1) ? i - 2 : i + 2;
				if (hindrance[i] < minHindrance) {
					bestSide = i;
					minHindrance = hindrance[i];
					counterHindrance = hindrance[counterSide];
					}
				else if (hindrance[i] == minHindrance) {
					if (hindrance[counterSide] > counterHindrance) {
						bestSide = i;
						counterHindrance = hindrance[counterSide];
						}
					}
				}

			switch (bestSide) {
			case 0:
				chax = getAtomX(atom);
				chay = getAtomY(atom) - mpDotDiameter - labelWidth / 2;
				break;
			case 1:
				chax = getAtomX(atom) + mpDotDiameter + labelWidth / 2;
				chay = getAtomY(atom);
				break;
			case 2:
				chax = getAtomX(atom);
				chay = getAtomY(atom) + mpDotDiameter + labelWidth / 2;
				break;
			default: // 3
				chax = getAtomX(atom) - mpDotDiameter - labelWidth / 2;
				chay = getAtomY(atom);
				break;
				}

			if (unpairedElectrons == 1) {
				mpDrawDot(chax,chay,atom);
				}
			else {
				switch (bestSide) {
				case 0:
					xdiff = 2 * mpDotDiameter;
					ydiff = 0;
					chax -= mpDotDiameter;
					break;
				case 1:
					xdiff = 0;
					ydiff = 2 * mpDotDiameter;
					chay -= mpDotDiameter;
					break;
				case 2:
					xdiff = 2 * mpDotDiameter;
					ydiff = 0;
					chax -= mpDotDiameter;
					break;
				default: // 3
					xdiff = 0;
					ydiff = 2 * mpDotDiameter;
					chay -= mpDotDiameter;
					break;
					}

				mpDrawDot(chax, chay, atom);
				mpDrawDot(chax + xdiff, chay + ydiff, atom);
				}
			}

		if (mCurrentColor == COLOR_EXCLUDE_GROUP_FG)
			setColor_(COLOR_RESTORE_PREVIOUS);
		}

	private String createRingSizeText(long queryFeatures) {
		queryFeatures &= Molecule.cAtomQFNewRingSize;
		for (int i=0; i<AtomQueryFeatureDialogBuilder.RING_SIZE_VALUES.length; i++)
			if (queryFeatures == AtomQueryFeatureDialogBuilder.RING_SIZE_VALUES[i])
				return AtomQueryFeatureDialogBuilder.RING_SIZE_SHORT_TEXT[i];

		StringBuilder customOption = new StringBuilder("R");
		if ((queryFeatures & Molecule.cAtomQFRingSize0) != 0)
			customOption.append("0");
		if ((queryFeatures & Molecule.cAtomQFRingSize3) != 0)
			customOption.append("3");
		if ((queryFeatures & Molecule.cAtomQFRingSize4) != 0)
			customOption.append("4");
		if ((queryFeatures & Molecule.cAtomQFRingSize5) != 0)
			customOption.append("5");
		if ((queryFeatures & Molecule.cAtomQFRingSize6) != 0)
			customOption.append("6");
		if ((queryFeatures & Molecule.cAtomQFRingSize7) != 0)
			customOption.append("7");
		if ((queryFeatures & Molecule.cAtomQFRingSizeLarge) != 0)
			customOption.append("8");

		return customOption.toString();
		}

	private void mpSetNormalLabelSize() {
		setTextSize(mpLabelSize);
		}

	private void mpSetReducedLabelSize() {
		setTextSize((mpLabelSize*5+1)/6);
		}

	private void mpSetSmallLabelSize() {
		setTextSize((mpLabelSize*2+1)/3);
		}

    private double mpGetFreeSpaceAngle(int atom) {
            // returns the angle from the given atom that is furthest away
            // from any bond and stereo label
	    double[] angle = new double[mMol.getAllConnAtoms(atom)];
        for (int i=0; i<mMol.getAllConnAtoms(atom); i++)
            angle[i] = mMol.getBondAngle(atom, mMol.getConnAtom(atom, i));
        Arrays.sort(angle);
	    double maxMean = mpGetMeanAngle(angle, 0);
	    double maxVal = mpGetAngleESRScore(angle, 0, maxMean);
        for (int i=1; i<angle.length; i++) {
	        double mean = mpGetMeanAngle(angle, i);
	        double val = mpGetAngleESRScore(angle, i, mean);
            if (maxVal < val) {
                maxVal = val;
                maxMean = mean;
                }
            }
        return maxMean;
        }

    private double mpGetAngleESRScore(double[] angleList, int index, double meanAngle) {
            // initial score is the angle difference between associated bonds
	    double score = (index == 0) ?
                2.0 * Math.PI + angleList[0] - angleList[angleList.length-1]
              : angleList[index] - angleList[index-1];

            // subtract a penalty for colliding with the CIP label
        if (meanAngle > -Math.PI*2/3 && meanAngle < Math.PI/3.0)
            score -= 2*Math.cos(meanAngle + Math.PI/6.0);
        else
            // add a small score for angles on the top right
            score -= 0.5*Math.cos(meanAngle + Math.PI/6.0);

            return score;
        }
    
    private double mpGetMeanAngle(double[] angle, int index) {
        if (index > 0)
            return (angle[index] + angle[index-1]) / 2.0;

	    double mean = Math.PI + (angle[0] + angle[angle.length-1]) / 2.0;
        return (mean > Math.PI) ? mean - 2.0f * Math.PI : mean;
        }

    private String append(String a, String b) {
        return (a == null) ? b : (b == null) ? a : a+","+b;
        }

	private void mpDrawString(double x, double y, String str, boolean withTabu) {
		if (withTabu) {
			double strWidth,xdiff,ydiff;

			strWidth = getStringWidth(str);
			xdiff = strWidth / 2 + getTextSize() / 8;
			ydiff = getTextSize() / 2;
			if (str == "+" || str == "-")
				ydiff = ydiff * 2 / 3;

			mpTabuZone.add(new GenericRectangle(x-xdiff, y-ydiff, 2*xdiff, 2*ydiff));
			}

		if (!mIsValidatingView)
			drawString(str,x,y);
		}


	private void mpDrawDot(double x, double y, int atm) {
		mpTabuZone.add(new GenericRectangle(x-mpDotDiameter, y-mpDotDiameter,
											  2*mpDotDiameter, 2*mpDotDiameter));

		if (!mIsValidatingView) {
			mpDot.add(new DepictorDot(x, y, isHighlightedAtom(atm) ? COLOR_HILITE_BOND_FG : mAtomColor[atm]));
			}
		}


	private void mpDrawAllDots() {
		for (DepictorDot dot:mpDot) {
			setColor_(dot.color);
			drawDot(dot.x, dot.y);
			}
		setColor_(mStandardForegroundColor);
		}


	private boolean mpAlleneCenter(int atm) {
		if (mMol.getConnAtoms(atm) != 2) return false;
		for (int i=0; i<2; i++)
			if (mMol.getConnBondOrder(atm,i) != 2)
				return false;

		return true;
		}


    private void mpDrawBondQueryFeatures() {
        boolean textSizeChanged = false;
        for (int bond=0; bond<mMol.getBonds(); bond++) {
        	String label = null;
            if (mMol.isBondBridge(bond)) {
                int minAtoms = mMol.getBondBridgeMinSize(bond);
                int maxAtoms = mMol.getBondBridgeMaxSize(bond);
                label = (minAtoms == maxAtoms) ? "["+minAtoms+"]" : "["+minAtoms+":"+maxAtoms+"]";
            	}
            else if ((mMol.getBondQueryFeatures(bond) & Molecule.cBondQFAromState) != 0) {
            	label = ((mMol.getBondQueryFeatures(bond) & Molecule.cBondQFAromState) == Molecule.cBondQFAromatic) ? "a"
            		  : ((mMol.getBondQueryFeatures(bond) & Molecule.cBondQFRingState) == Molecule.cBondQFRing) ? "r!a" : "!a";
            	}
            else if ((mMol.getBondQueryFeatures(bond) & Molecule.cBondQFRingState) != 0) {
            	label = ((mMol.getBondQueryFeatures(bond) & Molecule.cBondQFRingState) == Molecule.cBondQFRing) ? "r" : "!r";
            	}

            int ringSize = (mMol.getBondQueryFeatures(bond) & Molecule.cBondQFRingSize) >> Molecule.cBondQFRingSizeShift;
            if (ringSize != 0)
            	label = ((label == null) ? "" : label) + ringSize;

            if (label != null) {
	            int atom1 = mMol.getBondAtom(0, bond);
	            int atom2 = mMol.getBondAtom(1, bond);
	            if (!textSizeChanged) {
					mpSetSmallLabelSize();
	                textSizeChanged = true;
	                }
	            double x = (getAtomX(atom1) + getAtomX(atom2)) / 2;
	            double y = (getAtomY(atom1) + getAtomY(atom2)) / 2;
	            double dx = getAtomX(atom2) - getAtomX(atom1);
	            double dy = getAtomY(atom2) - getAtomY(atom1);
	            double d = Math.sqrt(dx*dx+dy*dy);
	            double hw = 0.60 * getStringWidth(label);	// slightly larger than 0.5f to increase label distance from bond
	            double hh = 0.55 * getTextSize();
	            if (d != 0) {
	            	if (dx > 0)
	            		mpDrawString(x+hw*dy/d, y-hh*dx/d, label, true);
	            	else
	            		mpDrawString(x-hw*dy/d, y+hh*dx/d, label, true);
	            	}
            	}
            }

        if (textSizeChanged)
			mpSetNormalLabelSize();
        }


    private void drawLine(DepictorLine theLine, int atom1, int atom2) {
    	if (mMol.isBondForegroundHilited(mMol.getBond(atom1, atom2))) {
			setColor_(COLOR_HILITE_BOND_FG);
			drawBlackLine(theLine);
			setColor_(mStandardForegroundColor);
    		}
    	else if (mAtomColor[atom1] != mAtomColor[atom2]) {
    		drawColorLine(theLine, atom1, atom2);
    		}
    	else if (mAtomColor[atom1] != Molecule.cAtomColorNone) {
			setColor_(mAtomColor[atom1]);
			drawBlackLine(theLine);
			setColor_(mStandardForegroundColor);
    		}
		else {
			drawBlackLine(theLine);
			}
		}


	private void drawColorLine(DepictorLine theLine,int atm1,int atm2) {
		DepictorLine line1 = new DepictorLine();
		DepictorLine line2 = new DepictorLine();
		line1.x1 = theLine.x1;
		line1.y1 = theLine.y1;
		line1.x2 = (theLine.x1 + theLine.x2) / 2;
		line1.y2 = (theLine.y1 + theLine.y2) / 2;
		line2.x1 = line1.x2;
		line2.y1 = line1.y2;
		line2.x2 = theLine.x2;
		line2.y2 = theLine.y2;
		if (mpProperLine(line1)) {
			setColor_(mAtomColor[atm1]);
			drawBlackLine(line1);
			}
		if (mpProperLine(line2)) {
			setColor_(mAtomColor[atm2]);
			drawBlackLine(line2);
			}
		setColor_(mStandardForegroundColor);
		}


	private void drawDashedLine(DepictorLine theLine, int atom1, int atom2) {
		double xinc = (theLine.x2 - theLine.x1) / 10;
		double yinc = (theLine.y2 - theLine.y1) / 10;

		DepictorLine aLine = new DepictorLine();

		int color1,color2;
		if (mMol.isBondForegroundHilited(mMol.getBond(atom1, atom2))) {
			color1 = COLOR_HILITE_BOND_FG;
			color2 = COLOR_HILITE_BOND_FG;
			}
		else {
			color1 = mAtomColor[atom1];
			color2 = mAtomColor[atom2];
			}

		setColor_(color1);
		aLine.x1 = theLine.x1;
		aLine.y1 = theLine.y1;
		aLine.x2 = theLine.x1 + xinc * 2;
		aLine.y2 = theLine.y1 + yinc * 2;
		drawBlackLine(aLine);
		aLine.x1 = theLine.x1 + xinc * 4;
		aLine.y1 = theLine.y1 + yinc * 4;
		aLine.x2 = theLine.x1 + xinc * 5;
		aLine.y2 = theLine.y1 + yinc * 5;
		drawBlackLine(aLine);

		setColor_(color2);
		aLine.x1 = theLine.x1 + xinc * 5;
		aLine.y1 = theLine.y1 + yinc * 5;
		aLine.x2 = theLine.x1 + xinc * 6;
		aLine.y2 = theLine.y1 + yinc * 6;
		drawBlackLine(aLine);
		aLine.x1 = theLine.x1 + xinc * 8;
		aLine.y1 = theLine.y1 + yinc * 8;
		aLine.x2 = theLine.x2;
		aLine.y2 = theLine.y2;
		drawBlackLine(aLine);

		setColor_(mStandardForegroundColor);
		}


	private void drawShortDashedLine(DepictorLine theLine, int atom1, int atom2) {
		double xdif = theLine.x2 - theLine.x1;
		double ydif = theLine.y2 - theLine.y1;
		double length = Math.sqrt(xdif * xdif + ydif * ydif);
		int points = 2 * (int)Math.round(length / (4 * mpLineWidth));

		double xinc = xdif / (points-1);
		double yinc = ydif / (points-1);

		int color1,color2;
		if (mMol.isBondForegroundHilited(mMol.getBond(atom1, atom2))) {
			color1 = COLOR_HILITE_BOND_FG;
			color2 = COLOR_HILITE_BOND_FG;
			}
		else {
			color1 = mAtomColor[atom1];
			color2 = mAtomColor[atom2];
			}

		double x = theLine.x1 - mpLineWidth/2;
		double y = theLine.y1 - mpLineWidth/2;

		setColor_(color1);
		for (int i=0; i<points/2; i++) {
			fillCircle(x, y, mpLineWidth);
			x += xinc;
			y += yinc;
			}

		setColor_(color2);
		for (int i=0; i<points/2; i++) {
			fillCircle(x, y, mpLineWidth);
			x += xinc;
			y += yinc;
			}

		setColor_(mStandardForegroundColor);
		}


	private void drawWedge(DepictorLine theWedge,int atom1, int atom2) {
		double xdiff = (theWedge.y1 - theWedge.y2) / 9;
		double ydiff = (theWedge.x2 - theWedge.x1) / 9;
		double xe1 = (theWedge.x2 + xdiff);
		double ye1 = (theWedge.y2 + ydiff);
		double xe2 = (theWedge.x2 - xdiff);
		double ye2 = (theWedge.y2 - ydiff);
		double xm1 = (theWedge.x1 + xe1) / 2;
		double ym1 = (theWedge.y1 + ye1) / 2;
		double xm2 = (theWedge.x1 + xe2) / 2;
		double ym2 = (theWedge.y1 + ye2) / 2;

		GenericPolygon p1 = new GenericPolygon(3);
		GenericPolygon p2 = new GenericPolygon(4);
		p1.addPoint(theWedge.x1, theWedge.y1);
		p1.addPoint(xm1, ym1);
		p1.addPoint(xm2, ym2);
		p2.addPoint(xm2, ym2);
		p2.addPoint(xm1, ym1);
		p2.addPoint(xe1, ye1);
		p2.addPoint(xe2, ye2);

		int color1,color2;
		if (mMol.isBondForegroundHilited(mMol.getBond(atom1, atom2))) {
			color1 = COLOR_HILITE_BOND_FG;
			color2 = COLOR_HILITE_BOND_FG;
			}
		else {
			color1 = mAtomColor[atom1];
			color2 = mAtomColor[atom2];
			if (mMol.getMoleculeColor() != Molecule.cMoleculeColorNeutral) {
				color2 = getESRColor(atom1);
				if (color1 == mMol.getAtomColor(atom1))    // if it is not selected or stereo error
					color1 = color2;
				}
			}

        setColor_(color1);
		drawPolygon(p1);
		setColor_(color2);
		drawPolygon(p2);
		setColor_(mStandardForegroundColor);
		}


	protected void drawDot(double x, double y) {
		fillCircle(x-mpDotDiameter/2, y-mpDotDiameter/2, mpDotDiameter);
		}


	private void setRGBColor(int rgb) {
		if (mOverruleForeground != 0) {
			if (mCurrentColor != COLOR_OVERRULED) {
				mCurrentColor = COLOR_OVERRULED;
			    setRGB(mOverruleForeground);
				}
			return;
			}

		mCurrentColor = COLOR_RGB;
		mRGBColor = rgb;
		setRGB(rgb);
		}


	private void setColor_(int theColor) {
		if (mIsValidatingView)
			return;

		if (theColor == COLOR_INITIALIZE) {
			mCurrentColor = -999;
			theColor = mStandardForegroundColor;
			}

		if (theColor != COLOR_HILITE_BOND_BG
		 && theColor != COLOR_EXCLUDE_GROUP_BG
		 && mOverruleForeground != 0)
			theColor = COLOR_OVERRULED;

		if (theColor == mCurrentColor)
			return;

		// if we have COLOR_EXCLUDE_GROUP_FG, then don't change until we get COLOR_RESTORE_PREVIOUS
		if (mCurrentColor == COLOR_EXCLUDE_GROUP_FG
		 && theColor != COLOR_RESTORE_PREVIOUS)
			return;

		if (theColor == COLOR_EXCLUDE_GROUP_FG)
			mPreviousColor = mCurrentColor;
		if (theColor == COLOR_RESTORE_PREVIOUS)
			theColor = mPreviousColor;

		mCurrentColor = theColor;

		switch (theColor) {
		case Molecule.cAtomColorNone:
			setRGB(mCustomForeground == 0 ? 0xFF000000 : mCustomForeground);
			break;
		case COLOR_CUSTOM_FOREGROUND:
			setRGB(mCustomForeground);
			break;
		case COLOR_OVERRULED:
			setRGB(mOverruleForeground);
			break;
		case COLOR_HILITE_BOND_BG:
			setRGB(mBondBGHiliteColor);
			break;
		case COLOR_HILITE_BOND_FG:
			setRGB(mBondFGHiliteColor);
			break;
		case COLOR_EXCLUDE_GROUP_BG:
			setRGB(mExcludeGroupBGColor);
			break;
		case COLOR_EXCLUDE_GROUP_FG:
			setRGB(mExcludeGroupFGColor);
			break;
		case COLOR_RGB:
			setRGB(mRGBColor);
			break;
		case Molecule.cAtomColorBlue:
			setRGB(COLOR_BLUE);
			break;
		case Molecule.cAtomColorRed:
		    setRGB(COLOR_RED);
			break;
		case Molecule.cAtomColorMagenta:
			setRGB(COLOR_MAGENTA);
			break;
		case Molecule.cAtomColorGreen:
			setRGB(COLOR_GREEN);
			break;
		case Molecule.cAtomColorOrange:
			setRGB(COLOR_ORANGE);
			break;
		case Molecule.cAtomColorDarkGreen:
			setRGB(COLOR_DARK_GREEN);
			break;
		case Molecule.cAtomColorDarkRed:
			setRGB(COLOR_DARK_RED);
			break;
		case cColorGray:
			setRGB(0xFF808080);
			break;
		default:
		    setRGB(0xFF000000);
			break;
			}
		}


	/**
	 * Determines the ESR type and group that should be displayed on atom, which is either
	 * a tetrahedral stereo center, the central atom of a chiral allene or that atom of a
	 * BINAP kind of axial chiral bond, which is carrying the indicating stereo bond.
	 * @param atom
	 * @return ESR-type + 256* ESR-group; -1 if nothing to display
	 */
	private int getESRTypeToDisplayAt(int atom) {
		int type = -1;
		int group = -1;

		if ((mDisplayMode & cDModeSuppressESR) != 0)
			return type;

		if (mMol.isAtomStereoCenter(atom)) {	// covers tetrahedral and axial allene type
			type = mMol.getAtomESRType(atom);
			group = mMol.getAtomESRGroup(atom);
			}

		int bond = mMol.findBINAPChiralityBond(atom);	// covers axial BINAP type
		if (bond != -1) {
			type = mMol.getBondESRType(bond);
			group = mMol.getBondESRGroup(bond);
			}

		if (type != -1 && type != Molecule.cESRTypeAbs)
			type |= (group << 8);

		return type;
		}


	/**
	 * Determines the ESR color, that is used to draw the ESR label or up/down bond.
	 * @param atom tetrahedral or allene stereocenter, allene end atom or atom of a BINAP bond
	 * @return
	 */
	private int getESRColor(int atom) {
		if ((mDisplayMode & (cDModeSuppressESR | cDModeNoColorOnESRAndCIP)) != 0)
			return mAtomColor[atom];

		int esrInfo = getESRTypeToDisplayAt(atom);
		if (esrInfo == -1) {
			int alleneCenter = mMol.findAlleneCenterAtom(atom);
			if (alleneCenter != -1) {
				atom = alleneCenter;
				esrInfo = getESRTypeToDisplayAt(atom);
				}
			}

		if (esrInfo == -1)
			return mAtomColor[atom];

        switch (esrInfo & 0x00FF) {
        case Molecule.cESRTypeAnd:
            return Molecule.cAtomColorDarkGreen;
        case Molecule.cESRTypeOr:
            return Molecule.cAtomColorBlue;
        default: // Molecule.cAtomESRTypeAbs
            return Molecule.cAtomColorDarkRed;
            }
        }

	protected abstract void drawBlackLine(DepictorLine theLine);
    protected abstract void drawDottedLine(DepictorLine theLine);
	protected abstract void drawPolygon(GenericPolygon p);
	protected abstract void drawString(String theString,double x,double y);
	protected abstract void fillCircle(double x, double y, double d);
	protected abstract double getLineWidth();
	protected abstract double getStringWidth(String theString);
    protected abstract int getTextSize();
	protected abstract void setTextSize(int theSize);
	protected abstract void setLineWidth(double lineWidth);
	protected abstract void setRGB(int rgb);

    public static class DepictorDot {
        public double x,y;
        public int color;

        DepictorDot(double x, double y, int color) {
            this.x = x;
            this.y = y;
            this.color = color;
        }
    }


    public static class DepictorLine {
        public double x1;
        public double y1;
        public double x2;
        public double y2;

        public DepictorLine(double x1, double y1, double x2, double y2)
        {
            this.x1 = x1;
            this.y1 = y1;
            this.x2 = x2;
            this.y2 = y2;
        }

        public DepictorLine()
        {

        }
    }

}


