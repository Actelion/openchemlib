/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.chem;

import com.actelion.research.util.ColorHelper;

import java.awt.*;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;

public abstract class AbstractDepictor {
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

	private static final Color BOND_BG_HILITE_COLOR = new Color(204,222,255);
    private static final Color BOND_FG_HILITE_COLOR = new Color(255, 85,  0);
    private static final float BOND_FG_HILITE_HUE = 20f/360f;	// orange-red

	// the minimum perceived brightness difference of color atom labels to the background
    private static final float MIN_COLOR_BG_CONTRAST = 0.4f;
   
    private static final int COLOR_UNDEFINED = -1;
    private static final int COLOR_HILITE_BOND_BG = -2;
    private static final int COLOR_HILITE_BOND_FG = -3;
    private static final int COLOR_OVERRULED = -4;
    private static final int RGB_COLOR = -5;

    public static final int cOptAvBondLen = 24;
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

	private static final int cDModeShowSymmetryAny = 0x0700;
	public static final int cDModeShowSymmetrySimple = 0x0100;
    public static final int cDModeShowSymmetryDiastereotopic = 0x0200;
    public static final int cDModeShowSymmetryEnantiotopic = 0x0400;
	public static final int	cDModeNoImplicitAtomLabelColors = 0x0800;
	public static final int	cDModeNoStereoProblem = 0x1000;

	private static final float cFactorTextSize = 0.6f;
	private static final float cFactorChiralTextSize = 0.5f;
	private static final float cFactorBondSpacing = 0.15f;
	private static final float cFactorBondHiliting = 0.75f;
	private static final float cFactorDotDiameter = 0.12f;
	private static final float cFactorQFDiameter = 0.40f;
	private static final float cFactorLineWidth = 0.06f;

	private boolean[]				mAtomIsConnected;
	private boolean[]				mAtomLabelDisplayed;
	private float					mpBondSpacing,mpBondHiliting,mpDotDiameter,mpLineWidth,mpQFDiameter,mFactorTextSize;
	private int						mpLabelSize,mDefaultColor,mDisplayMode,mCurrentColor;
	private ArrayList<Rectangle2D.Float> mpTabuZone;
    private ArrayList<DepictorDot>  mpDot;
	private StereoMolecule     		mMol;
	private Rectangle2D.Float		mBoundingRect = new Rectangle2D.Float();
	private DepictorTransformation	mTransformation;

	private Point2D.Float			mChiralTextLocation;
	private float					mChiralTextSize;
	private int[]					mAtomColor;
	private String[]				mAtomText;
	private Point2D.Float[]			mAlternativeCoords;
	private Color					mOverruleForeground,mOverruleBackground,mBondBGHiliteColor,mBondFGHiliteColor;

	protected Object				mG;

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


	public void setDefaultColor(int c) {
		mDefaultColor = c;
	    updateBondHiliteColor();
		}



	/**
	 * If the overrule color is set, the entire molecule is drawn in the foreground
	 * color neglecting any atom color information. The background color is used
	 * to construct a proper bond hilite color, if bond hiliting is used.
	 * @param foreground null or color to be used for molecule drawing
	 * @param background may be null
	 */
	public void setOverruleColor(Color foreground, Color background) {
	    mOverruleForeground = foreground;
	    mOverruleBackground = background;
	    updateBondHiliteColor();
		}

	
	public void setTransformation(DepictorTransformation t) {
		mTransformation = t;
		}


	/**
	 * Sets a multiplication factor to the text size of all labels. The default is 1.0.
	 * @param factor text size factor
	 */
	public void setFactorTextSize(float factor) {
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
	 * Returns full transformation that moves/scales original molecule into viewRect.
	 * This method considers atom labels when generating the bounding box. This includes
	 * atom labels of drawn implicit hydrogen atoms.
	 * @param g
	 * @param viewRect
	 * @param mode is typically (cModeInflateToMaxAVBL | maximum_desired_bond_length)
	 * @return
	 */
	public DepictorTransformation updateCoords(Graphics g, Rectangle2D.Float viewRect, int mode) {
		validateView(g, viewRect, mode);
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
	 * Returns full transformation that moves/scales original molecule into viewRect.
	 * This simple method creates a transformation that places all atom coordinates in
	 * the viewRect.
	 * @param viewRect
	 * @param mode
	 * @return
	 */
	public DepictorTransformation simpleUpdateCoords(Rectangle2D.Float viewRect, int mode) {
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


	/**
	 * A depictor maintains a DepictorTransformation object, which defines translation and scaling
	 * of molecule coordinates into the viewRect. This method updates the depictor's transformation
	 * such that the molecule is centered in viewRect and reduced in size if it doesn't fit.
	 * If mode contains cModeInflateToMaxAVBL, then it is scaled to reach the desired mean bond length,
	 * or to fill viewRect, whatever is reached first.
	 * @param g
	 * @param viewRect
	 * @param mode (<0> or cModeInflateToMaxAVBL) + <desired mean bond length>
	 * @return incremental transformation being applied to depictor's current transformation 
	 */
	public DepictorTransformation validateView(Object g, Rectangle2D.Float viewRect, int mode) {
		if (mMol.getAllAtoms() == 0)
			return null;

		DepictorTransformation t1 = simpleValidateView(viewRect, mode);

		mMol.ensureHelperArrays(requiredHelperArrays());
		markIsolatedAtoms();

		mpDot.clear();
		mpTabuZone.clear();

	    mG = g;
		calculateParameters();

		Font oldFont = (g instanceof Graphics) ? ((Graphics)g).getFont() : null;
		setTextSize(mpLabelSize);

		for (int i=0; i<mMol.getAllAtoms(); i++)
	    	mpDrawAtom(i, false);

		if (g instanceof Graphics)
		    ((Graphics)g).setFont(oldFont);

		float avbl = mTransformation.getScaling() * mMol.getAverageBondLength();
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


	protected DepictorTransformation simpleValidateView(Rectangle2D.Float viewRect, int mode) {
	// returns incremental transformation that moves/scales already transformed molecule into viewRect
		if (mMol.getAllAtoms() == 0)
			return null;

		simpleCalculateBounds();

		float avbl = mTransformation.getScaling() * mMol.getAverageBondLength();
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
    protected void onDrawBond(int atom1, int atom2,float x1, float y1,float x2, float y2)
    {
        // NOOP
    }

    // This might be overridden by subclasses (e.g. SVG Depictor)
    protected void onDrawAtom(int atom, String symbol,float x, float y)
    {
        // NOOP
    }


    private void simpleCalculateBounds() {
		float minx = getAtomX(0);	// determine size of molecule
		float maxx = getAtomX(0);
		float miny = getAtomY(0);
		float maxy = getAtomY(0);

		for (int i=1; i<mMol.getAllAtoms(); i++) {
			if (getAtomX(i) < minx) minx = getAtomX(i);
			if (getAtomX(i) > maxx) maxx = getAtomX(i);
			if (getAtomY(i) < miny) miny = getAtomY(i);
			if (getAtomY(i) > maxy) maxy = getAtomY(i);
			}

		mBoundingRect = new Rectangle2D.Float(minx, miny, maxx-minx, maxy-miny);
		}


	private void expandBoundsByTabuZones(float avbl) {
		for (int i=0; i<mpTabuZone.size(); i++)
			mBoundingRect = (Rectangle2D.Float)mBoundingRect.createUnion(mpTabuZone.get(i));

		float border = 0.1f * avbl;
		mBoundingRect.x -= border;
		mBoundingRect.y -= border;
		mBoundingRect.width += 2.0*border;
		mBoundingRect.height += 2.0*border;
		}


	private void setChiralTextLocation(Rectangle2D.Float viewRect, float avbl, int mode) {
		float spacing = avbl / 2.0f;
		switch (mode & cModeChiralTextLocation) {
		case cModeChiralTextOnFrameBottom:
			if (viewRect != null) {
				mChiralTextLocation.x = viewRect.x + viewRect.width/2.0f;
				mChiralTextLocation.y = viewRect.y + viewRect.height - spacing;
				break;
				}
		case cModeChiralTextBelowMolecule:
			mChiralTextLocation.x = mBoundingRect.x + mBoundingRect.width/2.0f;
			mChiralTextLocation.y = mBoundingRect.y + mBoundingRect.height + spacing;
			if (viewRect != null && mChiralTextLocation.y > viewRect.y + viewRect.height - spacing)
				mChiralTextLocation.y = viewRect.y + viewRect.height - spacing;
			break;
		case cModeChiralTextOnFrameTop:
			if (viewRect != null) {
				mChiralTextLocation.x = viewRect.x + viewRect.width/2.0f;
				mChiralTextLocation.y = viewRect.y + spacing;
				break;
				}
		case cModeChiralTextAboveMolecule:
			mChiralTextLocation.x = mBoundingRect.x + mBoundingRect.width/2.0f;
			mChiralTextLocation.y = mBoundingRect.y - spacing;
			if (viewRect != null && mChiralTextLocation.y < viewRect.y + spacing)
				mChiralTextLocation.y = viewRect.y + spacing;
			break;
			}
		}


	private float getAtomX(int atom) {
		return mTransformation.transformX(mMol.getAtomX(atom));
		}


	private float getAtomY(int atom) {
		return mTransformation.transformY(mMol.getAtomY(atom));
		}


	public final Rectangle2D.Float getBoundingRect() {
			// requires a prior call of updateCoords() or validateView()
			// returns the bounding rectangle in device coordinates (of the moved/scaled molecule)
		return mBoundingRect;
		}


	protected void init() {
		mFactorTextSize = 1.0f;
		mTransformation = new DepictorTransformation();
		mpTabuZone = new ArrayList<Rectangle2D.Float>();
		mpDot = new ArrayList<DepictorDot>();
		mAtomLabelDisplayed = new boolean[mMol.getAllAtoms()];
		mChiralTextLocation = new Point2D.Float();
		mDefaultColor = Molecule.cAtomColorBlack;
		mCurrentColor = COLOR_UNDEFINED;
		mBondBGHiliteColor = BOND_BG_HILITE_COLOR;
		mBondFGHiliteColor = BOND_FG_HILITE_COLOR;
		}


    private void updateBondHiliteColor() {
    	if (mOverruleBackground == null) {
    		mBondBGHiliteColor = BOND_BG_HILITE_COLOR;
    		mBondFGHiliteColor = BOND_FG_HILITE_COLOR;
    		}
    	else {
    		Color foreground = (mOverruleForeground == null) ? Color.BLACK : mOverruleForeground;
    		Color background = (mOverruleBackground == null) ? Color.WHITE : mOverruleBackground;

    		mBondBGHiliteColor = new Color((3*foreground.getRed()+7*background.getRed())/10,
    									   (3*foreground.getGreen()+7*background.getGreen())/10,
    									   (3*foreground.getBlue()+7*background.getBlue())/10);
    		float brightness = ColorHelper.perceivedBrightness(mOverruleBackground);
    		mBondFGHiliteColor = Color.getHSBColor(BOND_FG_HILITE_HUE, 1.0f, brightness > 0.5 ? 0.5f : 1.0f);
    		}
    	}


    private void calculateParameters() {
		float averageBondLength = mTransformation.getScaling() * mMol.getAverageBondLength();
		mpLineWidth = averageBondLength * cFactorLineWidth;
		mpBondSpacing = averageBondLength * cFactorBondSpacing;
		mpBondHiliting = averageBondLength * cFactorBondHiliting;
		mpLabelSize    = (int)(averageBondLength * mFactorTextSize * cFactorTextSize + 0.5);
		mpDotDiameter = averageBondLength * cFactorDotDiameter;
		mpQFDiameter = averageBondLength * cFactorQFDiameter;
		mChiralTextSize = averageBondLength * cFactorChiralTextSize + 0.5f;
		}

	public synchronized void paint(Object g) {
		if (mMol.getAllAtoms() == 0)
		    return;

		mMol.ensureHelperArrays(requiredHelperArrays());

		mG = g;
		calculateParameters();

		boolean explicitAtomColors = false;
		mAtomColor = new int[mMol.getAllAtoms()];
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			mAtomColor[atom] = mMol.getAtomColor(atom);
			if (mAtomColor[atom] != Molecule.cAtomColorBlack)
				explicitAtomColors = true;
			if (mMol.isSelectedAtom(atom))
				mAtomColor[atom] = Molecule.cAtomColorRed;
			if (mMol.getStereoProblem(atom) && (mDisplayMode & cDModeNoStereoProblem) == 0)
				mAtomColor[atom] = Molecule.cAtomColorMagenta;
			}

		Font oldFont = (g instanceof Graphics) ? ((Graphics)g).getFont() : null;

		hiliteBondBackgrounds();
		indicateQueryFeatures();
		addChiralInfo();

		setTextSize(mpLabelSize);
		setLineWidth(mpLineWidth);
		setColor(mDefaultColor);

		markIsolatedAtoms();

		mpDot.clear();
		mpTabuZone.clear();

		for (int i=0; i<mMol.getAllAtoms(); i++) {
			if (isHighlightedAtom(i)) {
				setColor(COLOR_HILITE_BOND_FG);
	    		mpDrawAtom(i, true);
				setColor(mDefaultColor);
				}
			else if (mAtomColor[i] != 0) {
				setColor(mAtomColor[i]);
	    		mpDrawAtom(i, true);
				setColor(mDefaultColor);
				}
			else if (!explicitAtomColors
				  && ((mDisplayMode & cDModeNoImplicitAtomLabelColors) == 0)
				  && mMol.getAtomicNo(i) < ATOM_LABEL_COLOR.length) {
				setRGBColor(getContrastColor(ATOM_LABEL_COLOR[mMol.getAtomicNo(i)]));
	    		mpDrawAtom(i, true);
				setColor(mDefaultColor);
				}
			else {
	    		mpDrawAtom(i, true);
				}
			}
		mpDrawAllDots();
        mpDrawBondQueryFeatures();
		mpDrawAllBonds();

		if (g instanceof Graphics)
		    ((Graphics)g).setFont(oldFont);
		}


	private Color getContrastColor(int rgb) {
		Color bg = (mOverruleBackground == null) ? Color.WHITE : mOverruleBackground;
		float bgb = ColorHelper.perceivedBrightness(bg);
		Color fg = new Color(rgb);
		float fgb = ColorHelper.perceivedBrightness(fg);

		// minConstrast is MIN_COLOR_BG_CONTRAST with white or black background and reduces
		// to 0.5*MIN_COLOR_BG_CONTRAST when background brightness goes to 0.5
		float minContrast = MIN_COLOR_BG_CONTRAST * (0.5f + Math.abs(0.5f - bgb));

		float b1 = bgb - minContrast;	// lower edge of brightness avoidance zone
		float b2 = bgb + minContrast;	// higher edge of brightness avoidance zone
		boolean darken = (b1 <= 0f) ? false : (b2 >= 1.0f) ? true : fgb < bgb;
		float factor = 1f / 255f * ((darken) ?
				  1f-minContrast * brightnessShiftFunction(Math.max((bgb - fgb) / bgb, 0f))
				: 1f+minContrast * brightnessShiftFunction(Math.max((fgb - bgb) / (1f - bgb), 0f)));

		return new Color(factor * fg.getRed(), factor * fg.getGreen(), factor * fg.getBlue());
		}


	/**
	 * Nonlinear function for adjusting color brightness depending on similarity to background brightness
	 * @param d distance between color brightness and background brightness (0 ... 1)
	 * @return
	 */
	private float brightnessShiftFunction(float d) {
		return (1.0f/(0.5f+1.5f*d)-0.5f)/1.5f;
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
	         : ((mDisplayMode & cDModeShowSymmetryDiastereotopic) != 0) ? Molecule.cHelperSymmetryDiastereotopic
             : ((mDisplayMode & cDModeShowSymmetryEnantiotopic) != 0) ? Molecule.cHelperSymmetryEnantiotopic
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
				float avbl = mTransformation.getScaling() * mMol.getAverageBondLength();
				simpleCalculateBounds();
				expandBoundsByTabuZones(avbl);
				setChiralTextLocation(null, avbl, 0);
				}

			setTextSize((int)mChiralTextSize);
			setColor(Molecule.cAtomColorRed);
			drawString(chiralText, mChiralTextLocation.x, mChiralTextLocation.y+0.3f*mChiralTextSize);
			}
		}


	private void hiliteBondBackgrounds() {
        setColor(COLOR_HILITE_BOND_BG);
        setLineWidth(mpBondHiliting);
        DepictorLine line = new DepictorLine();
        for (int bond=0; bond<mMol.getAllBonds(); bond++) {
        	if (mMol.isBondBackgroundHilited(bond)) {
	        	line.x1 = getAtomX(mMol.getBondAtom(0, bond));
	        	line.y1 = getAtomY(mMol.getBondAtom(0, bond));
	            line.x2 = getAtomX(mMol.getBondAtom(1, bond));
	            line.y2 = getAtomY(mMol.getBondAtom(1, bond));
	            drawBlackLine(line);
        		}
        	}
		}


	private void indicateQueryFeatures() {
		if (mMol.isFragment()) {
			setColor(Molecule.cAtomColorOrange);
			if (((mDisplayMode & cDModeHiliteAllQueryFeatures) != 0))
				for (int atom=0; atom<mMol.getAtoms(); atom++)
					if (mMol.getAtomQueryFeatures(atom) != 0)
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


	private void mpDrawAllBonds() {
		mAlternativeCoords = new Point2D.Float[mMol.getAllAtoms()];

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
					String infoStr;
					switch (mMol.getBondCIPParity(i)) {
					case Molecule.cBondCIPParityEorP:
						infoStr = (mMol.getBondOrder(i) == 2) ? "E" : mMol.isBondParityPseudo(i) ? "p" : "P";
						break;
					case Molecule.cBondCIPParityZorM:
						infoStr = (mMol.getBondOrder(i) == 2) ? "Z" : mMol.isBondParityPseudo(i) ? "m" : "M";
						break;
					default:
						infoStr = "?";
						break;
						}
					setTextSize((mpLabelSize*2+1)/3);
					setColor(mMol.isBondForegroundHilited(i) ? COLOR_HILITE_BOND_FG : Molecule.cAtomColorRed);
					int atom1 = mMol.getBondAtom(0,i);
					int atom2 = mMol.getBondAtom(1,i);
					float x = (getAtomX(atom1) + getAtomX(atom2)) / 2;
					float y = (getAtomY(atom1) + getAtomY(atom2)) / 2;
					float dx = (getAtomX(atom1) - getAtomX(atom2)) / 3;
					float dy = (getAtomY(atom1) - getAtomY(atom2)) / 3;
					mpDrawString(x+dy,y-dx,infoStr,true,true);
					setColor(mDefaultColor);
					setTextSize(mpLabelSize);
					}
				}
			}

		if ((mDisplayMode & cDModeBondNo) != 0) {
			setTextSize((mpLabelSize*2+1)/3);
			setColor(Molecule.cAtomColorDarkGreen);
			for (int i=0; i<mMol.getAllBonds(); i++) {
				int atom1 = mMol.getBondAtom(0,i);
				int atom2 = mMol.getBondAtom(1,i);
				float x = (getAtomX(atom1) + getAtomX(atom2)) / 2;
				float y = (getAtomY(atom1) + getAtomY(atom2)) / 2;
				mpDrawString(x,y,String.valueOf(i),true,true);
				}
			setColor(mDefaultColor);
			setTextSize(mpLabelSize);
			}
		}


	private void mpDrawBond(int bnd) {
		DepictorLine theLine = new DepictorLine();
		DepictorLine aLine = new DepictorLine();
		DepictorLine bLine = new DepictorLine();
		Point2D.Float piBondOffset = new Point2D.Float();
		Point2D.Float nextBondOffset = new Point2D.Float();

		int atom1 = mMol.getBondAtom(0,bnd);
		int atom2 = mMol.getBondAtom(1,bnd);

        onDrawBond(atom1,atom2,getAtomX(atom1),getAtomY(atom1),getAtomX(atom2),getAtomY(atom2));

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
            return;
            }

        int bondOrder = (mMol.getBondType(bnd) == Molecule.cBondTypeDelocalized) ? 0
                      : (mMol.getBondType(bnd) == Molecule.cBondTypeMetalLigand) ? 1
					  : mMol.getBondOrder(bnd);

		switch (bondOrder) {
		case 1:
            switch (mMol.getBondType(bnd)) {
			case Molecule.cBondTypeSingle:
				mpHandleLine(theLine, atom1, atom2);
				break;
			case Molecule.cBondTypeUp:
				mpHandleWedge(theLine, atom1, atom2);
				break;
			case Molecule.cBondTypeDown:
				float xdiff = theLine.x2 - theLine.x1;
				float ydiff = theLine.y2 - theLine.y1;

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
						setColor((i<9) ? color1 : color2);
						drawBlackLine(aLine);
						setColor(mDefaultColor);
						}
					}
				break;
            case Molecule.cBondTypeMetalLigand:
                mpHandleDashedLine(theLine, atom1, atom2);
                break;
			}
			break;
		case 0:	// bonds defined to be aromatic
		case 2:
			if ((mAtomLabelDisplayed[atom1] || mMol.getAtomPi(atom1) == 2)
		 	 && (mAtomLabelDisplayed[atom2] || mMol.getAtomPi(atom2) == 2)
		 	 && !mMol.isRingBond(bnd)
		 	 && bondOrder == 2) {
											// float bond connecting two atoms both being either
											// a hetero atom or the central atom of an allene
				if (!mpProperLine(theLine))
					break;

				mpCalcPiBondOffset(theLine.x2 - theLine.x1,
								   theLine.y2 - theLine.y1,piBondOffset);
				float xdiff = piBondOffset.x / 2;
				float ydiff = piBondOffset.y / 2;

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
											// C=X float bond with atm1 is the carbon
											// or R2C=C=CR2 with atm2 is central atom
				mpDBFromNonLabelToLabel(theLine, bnd, false);
				}
			else if ((mAtomLabelDisplayed[atom1] || mMol.getAtomPi(atom1) == 2)
		          && bondOrder == 2) {
											// C=X float bond with atm2 is the carbon
											// or R2C=C=CR2 with atm1 is central atom
				mpDBFromNonLabelToLabel(theLine, bnd, true);
				}
			else {
								// standard carbon-carbon float bond. Thus,one bond
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
								   theLine.y2 - theLine.y1,piBondOffset);
				aLine.x1 = theLine.x1 + piBondOffset.x;
				aLine.y1 = theLine.y1 + piBondOffset.y;
				aLine.x2 = theLine.x2 + piBondOffset.x;
				aLine.y2 = theLine.y2 + piBondOffset.y;
				drawLine(aLine, atom1, atom2);
				aLine.x1 = theLine.x1 - piBondOffset.x;
				aLine.y1 = theLine.y1 - piBondOffset.y;
				aLine.x2 = theLine.x2 - piBondOffset.x;
				aLine.y2 = theLine.y2 - piBondOffset.y;
				drawLine(aLine, atom1, atom2);
				}
			break;
            }
		}


	private void mpDBFromNonLabelToLabel(DepictorLine theLine, int bnd, boolean inverted) {
		DepictorLine aLine = new DepictorLine();
		DepictorLine bLine = new DepictorLine();
		Point2D.Float piBondOffset = new Point2D.Float();
		Point2D.Float nextBondOffset = new Point2D.Float();

		int atm1 = mMol.getBondAtom(0,bnd);
		int atm2 = mMol.getBondAtom(1,bnd);

		if (inverted) {
			float td = theLine.x1;
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

		if (mMol.isRingBond(bnd)) {	// don't draw a centered float bond when
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
			float xdiff = piBondOffset.x / 2;
			float ydiff = piBondOffset.y / 2;

			boolean aLineIsInnerLine = false;
//			boolean bLineIsInnerLine = false;

			aLine.x1 = theLine.x1 + xdiff;
			aLine.y1 = theLine.y1 + ydiff;
			aLine.x2 = theLine.x2 + xdiff;
			aLine.y2 = theLine.y2 + ydiff;

			if (mMol.getConnAtoms(atm1) > 1) {
				if (!mpCalcNextBondOffset(atm1,atm2,1,nextBondOffset)) {
					mAlternativeCoords[atm1] = new Point2D.Float(aLine.x1, aLine.y1);
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
					mAlternativeCoords[atm1] = new Point2D.Float(bLine.x1, bLine.y1);
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
		float temp;
		temp = aLine.x2;
		aLine.x2 = bLine.x2;
		bLine.x2 = temp;
		temp = aLine.y2;
		aLine.y2 = bLine.y2;
		bLine.y2 = temp;
		}


	private void mpCalcPiBondOffset(float dx, float dy, Point2D.Float piBondOffset) {
		if (dx == 0) {
			if (dy < 0)
				piBondOffset.x =   mpBondSpacing;
			else
				piBondOffset.x = - mpBondSpacing;
			piBondOffset.y = 0;
			return;
			}

		float alpha = (float)Math.atan(dy / dx);
		if (dx < 0)
			alpha += Math.PI;
		piBondOffset.x = - (mpBondSpacing * (float)Math.sin(alpha));
		piBondOffset.y =   (mpBondSpacing * (float)Math.cos(alpha));
		}


	private boolean mpProperLine(DepictorLine theLine) {
		// cuts line ends according to needs of involved atoms and returns
		// 'false' if line lies entirely in tabuZones,otherwise it returns 'true'
		boolean endsExchanged,retval;

		if (theLine.x1 == theLine.x2 && theLine.y1 == theLine.y2) {
			for (int i=0; i<mpTabuZone.size(); i++) {
				Rectangle2D.Float tabuZone = mpTabuZone.get(i);
				if (tabuZone.contains(theLine.x1, theLine.y1))
					return false;
				}
			return true;
			}

		Rectangle2D.Float theFrame = mpGetFrame(theLine);

		endsExchanged = false;
		if (theLine.x1 > theLine.x2) {  // first point is the one with smaller x
			mpExchangeLineEnds(theLine);
			endsExchanged = true;
			}

		for (int i=0; i<mpTabuZone.size(); i++) {
			Rectangle2D.Float tabuZone = mpTabuZone.get(i);
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

	private boolean mpCalcNextBondOffset(int atm1,int atm2,int side,Point2D.Float nextBondOffset) {
		final float RO_LIMIT = 2.617993878f;	// right outer angle limit = 150 degrees
		final float LO_LIMIT = 3.665191429f;	// left  outer angle limit = 210 degrees
		final float RI_LIMIT = 0.523598776f;	// right inner angle limit =  30 degrees
		final float LI_LIMIT = 5.759586531f;	// left  inner angle limit = 330 degrees

		boolean retval;
		int i,remoteAtm,bnd;
		float bondAngle,theBondAngle,testAngle;
		float angleDiff,currentAngleDiff,distance;

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
								// a bond is leading away from float bond's right side

				if (currentAngleDiff > RO_LIMIT)
					currentAngleDiff = RO_LIMIT;

				if (currentAngleDiff < RI_LIMIT)
					currentAngleDiff = RI_LIMIT;

				if (currentAngleDiff <= angleDiff) {
					angleDiff = currentAngleDiff;
					distance = (float)mpBondSpacing * (float)Math.tan(angleDiff - Math.PI/2) / 2;
					nextBondOffset.x = - (distance * (float)Math.sin(bondAngle));
					nextBondOffset.y = - (distance * (float)Math.cos(bondAngle));
					}
				}
			else {
				if (currentAngleDiff >= Math.PI)
					retval = true;
								// a bond is leading away from float bond's left side

				if (currentAngleDiff < LO_LIMIT)
					currentAngleDiff = LO_LIMIT;

				if (currentAngleDiff > LI_LIMIT)
					currentAngleDiff = LI_LIMIT;

				if (currentAngleDiff >= angleDiff) {
					angleDiff = currentAngleDiff;
					distance = (float)mpBondSpacing * (float)Math.tan(4.712388981 - angleDiff) / 2;
					nextBondOffset.x = - (distance * (float)Math.sin(bondAngle));
					nextBondOffset.y = - (distance * (float)Math.cos(bondAngle));
					}
				}
			}
		return retval;
		}


	private void mpExchangeLineEnds(DepictorLine theLine) {
		float temp;
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

		Rectangle2D.Float theFrame = mpGetFrame(theWedge);

		for (int i=0; i<mpTabuZone.size(); i++) {
			Rectangle2D.Float tabuZone = mpTabuZone.get(i);
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


	private Rectangle2D.Float mpGetFrame(DepictorLine theLine) {
		Rectangle2D.Float theFrame = new Rectangle2D.Float();
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


	private boolean mpInTabuZone(float x,float y,int tabuZoneNo) {
		if ((mDisplayMode & cDModeNoTabus) != 0)
			return false;

		Rectangle2D.Float tabuZone = mpTabuZone.get(tabuZoneNo);

			// cannot use tabuZone.contains() because points on edge would be considered to be within the react
		return (x > tabuZone.x
			 && x < tabuZone.x + tabuZone.width
			 && y > tabuZone.y
			 && y < tabuZone.y + tabuZone.height);
		}


	private void mpShortenLine(DepictorLine theLine,int pointNo,int tabuZoneNo) {
		float x1,y1,x2,y2,dx,dy,tabuX,tabuY,sx,sy;

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

		Rectangle2D.Float tabuZone = mpTabuZone.get(tabuZoneNo);
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
		float[] angle = new float[ExtendedMolecule.cMaxConnAtoms];
		float[] bondAngle = new float[2];

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
			bondAngle[1] = bondAngle[0] + (float)Math.PI;
			changed = false;
			}
		else {
			bondAngle[1] = bondAngle[0];
			bondAngle[0] = bondAngle[1] - (float)Math.PI;
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


	private void mpDrawAtom(int atom, boolean drawAtoms) {
		float chax,chay,xdiff,ydiff,x,y;

        if (drawAtoms)
            onDrawAtom(atom,mMol.getAtomLabel(atom), getAtomX(atom), getAtomY(atom));


		String propStr = null;
		if (mMol.getAtomCharge(atom) != 0) {
			String valStr = (Math.abs(mMol.getAtomCharge(atom)) == 1) ? ""
			                    : String.valueOf(Math.abs(mMol.getAtomCharge(atom)));
			propStr = (mMol.getAtomCharge(atom) < 0) ? valStr + "-" : valStr + "+";
			}
		if (mAtomText != null && (atom < mAtomText.length) && mAtomText[atom] != null)
			propStr = append(propStr, mAtomText[atom]);

		String isoStr = null;
		int queryFeatures = mMol.getAtomQueryFeatures(atom);
		if (queryFeatures != 0) {
			if ((queryFeatures & Molecule.cAtomQFAromatic) != 0)
				isoStr = append(isoStr, "a");
			if ((queryFeatures & Molecule.cAtomQFNotAromatic) != 0)
				isoStr = append(isoStr, "!a");
			if ((queryFeatures & Molecule.cAtomQFMoreNeighbours) != 0)
				isoStr = append(isoStr, "s");
            if ((queryFeatures & Molecule.cAtomQFHydrogen) != 0) {
                int hydrogens = (queryFeatures & Molecule.cAtomQFHydrogen);
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
                else if (hydrogens == Molecule.cAtomQFNot3Hydrogen)
                    isoStr = append(isoStr, "h<3");
                else if (hydrogens == Molecule.cAtomQFNot2Hydrogen+Molecule.cAtomQFNot3Hydrogen)
                    isoStr = append(isoStr, "h<2");
                }
            if ((queryFeatures & Molecule.cAtomQFCharge) != 0) {
                int charge = (queryFeatures & Molecule.cAtomQFCharge);
    			if (charge == Molecule.cAtomQFNotChargePos+Molecule.cAtomQFNotChargeNeg)
    				isoStr = append(isoStr, "c0");
    			else if (charge == Molecule.cAtomQFNotCharge0+Molecule.cAtomQFNotChargeNeg)
    				isoStr = append(isoStr, "c+");
    			else if (charge == Molecule.cAtomQFNotCharge0+Molecule.cAtomQFNotChargePos)
    				isoStr = append(isoStr, "c-");
                }
            if ((queryFeatures & Molecule.cAtomQFPiElectrons) != 0) {
                int piElectrons = (queryFeatures & Molecule.cAtomQFPiElectrons);
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
                int neighbours = (queryFeatures & Molecule.cAtomQFNeighbours);
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
            if ((queryFeatures & Molecule.cAtomQFRingState) != 0) {
                int ringBonds = (queryFeatures & Molecule.cAtomQFRingState);
                if (ringBonds == Molecule.cAtomQFNot2RingBonds+Molecule.cAtomQFNot3RingBonds+Molecule.cAtomQFNot4RingBonds)
                    isoStr = append(isoStr, "c");
                else if (ringBonds == Molecule.cAtomQFNotChain)
                    isoStr = append(isoStr, "r");
                else if (ringBonds == Molecule.cAtomQFNotChain+Molecule.cAtomQFNot3RingBonds+Molecule.cAtomQFNot4RingBonds)
                    isoStr = append(isoStr, "rb2");
                else if (ringBonds == Molecule.cAtomQFNotChain+Molecule.cAtomQFNot2RingBonds+Molecule.cAtomQFNot4RingBonds)
                    isoStr = append(isoStr, "rb3");
                else if (ringBonds == Molecule.cAtomQFNotChain+Molecule.cAtomQFNot2RingBonds+Molecule.cAtomQFNot3RingBonds)
                    isoStr = append(isoStr, "rb4");
                }
            if ((queryFeatures & Molecule.cAtomQFRingSize) != 0) {
                isoStr = append(isoStr, "rs"+((queryFeatures & Molecule.cAtomQFRingSize)>>Molecule.cAtomQFRingSizeShift));
                }
            if ((queryFeatures & Molecule.cAtomQFFlatNitrogen) != 0) {
                isoStr = append(isoStr, "sp2");
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

		String infoStr = null;
		if ((mDisplayMode & cDModeSuppressCIPParity) == 0) {
			if (mMol.isAtomConfigurationUnknown(atom))
				infoStr = "?";
			else if (mMol.getAtomCIPParity(atom) != 0) {
				if (mMol.getConnAtoms(atom) == 2) {
					switch (mMol.getAtomCIPParity(atom)) {
					case Molecule.cAtomCIPParitySorP:
						infoStr = mMol.isAtomParityPseudo(atom) ? "p" : "P";
						break;
					case Molecule.cAtomCIPParityRorM:
						infoStr = mMol.isAtomParityPseudo(atom) ? "m" : "M";
						break;
					default:
						infoStr = "*";
						break;
						}
					}
				else {
					switch (mMol.getAtomCIPParity(atom)) {
					case Molecule.cAtomCIPParityRorM:
						infoStr = mMol.isAtomParityPseudo(atom) ? "r" : "R";
						break;
					case Molecule.cAtomCIPParitySorP:
						infoStr = mMol.isAtomParityPseudo(atom) ? "s" : "S";
						break;
					default:
						infoStr = "*";
						break;
						}
					}
				}
			}
        if ((mDisplayMode & cDModeShowSymmetryAny) != 0)
            infoStr = append(infoStr, String.valueOf(mMol.getSymmetryRank(atom)));

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
		if (mMol.getAtomicNo(atom) != 6
		 || !mAtomIsConnected[atom]
		 || (mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFNoMoreNeighbours) != 0
		 || mMol.getAtomRadical(atom) != 0)
			hydrogensToAdd = mMol.getImplicitHydrogens(atom);

		String atomStr = mMol.getAtomCustomLabel(atom);
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

		float labelWidth = 0.0f;

		if (atomStr != null) {
			labelWidth = getStringWidth(atomStr);
			mpDrawString(getAtomX(atom),getAtomY(atom),atomStr,drawAtoms,true);
			mAtomLabelDisplayed[atom] = true;
			}
		else if (mpAlleneCenter(atom))
			mpDrawDot(getAtomX(atom),getAtomY(atom),atom,drawAtoms);

		if (propStr != null) {
			setTextSize((mpLabelSize*2+1)/3);
			x = getAtomX(atom) + ((labelWidth + getStringWidth(propStr)) / 2.0f + 1);
			y = getAtomY(atom) - ((getTextSize()*4-4)/8);
			mpDrawString(x,y,propStr,drawAtoms,true);
			setTextSize(mpLabelSize);
			}

        if ((mDisplayMode & cDModeAtomNo) != 0)
            isoStr = String.valueOf(atom);

		if (isoStr != null) {
			setTextSize((mpLabelSize*2+1)/3);
			x = getAtomX(atom) - ((labelWidth + getStringWidth(isoStr)) / 2.0f);
			y = getAtomY(atom) - ((getTextSize()*4-4)/8);
			mpDrawString(x,y,isoStr,drawAtoms,true);
			setTextSize(mpLabelSize);
			}

		if (infoStr != null) {
			setTextSize((mpLabelSize*2+1)/3);
			x = getAtomX(atom) - ((labelWidth + getStringWidth(infoStr)) / 2.0f);
			y = getAtomY(atom) + ((getTextSize()*4+4)/8);
			int theColor = mCurrentColor;
			setColor(Molecule.cAtomColorRed);
			mpDrawString(x,y,infoStr,drawAtoms,false);
			setColor(theColor);
			setTextSize(mpLabelSize);
			}

		if (mapStr != null) {
			setTextSize((mpLabelSize*2+1)/3);
			x = getAtomX(atom) + ((labelWidth + getStringWidth(mapStr)) / 2.0f + 1);
			y = getAtomY(atom) + ((getTextSize()*4+4)/8);
			int theColor = mCurrentColor;
			setColor(mMol.isAutoMappedAtom(atom) ? Molecule.cAtomColorDarkGreen : Molecule.cAtomColorDarkRed);
			mpDrawString(x,y,mapStr,drawAtoms,true);
			setColor(theColor);
			setTextSize(mpLabelSize);
			}

        if (esrStr != null) {
            float angle = mpGetFreeSpaceAngle(atom);
            setTextSize((mpLabelSize*2+1)/3);
            x = getAtomX(atom) + (float)getTextSize()*0.7f*(float)Math.sin(angle);
            y = getAtomY(atom) + (float)getTextSize()*0.7f*(float)Math.cos(angle);
            int theColor = mCurrentColor;
            setColor(getESRColor(atom));
            mpDrawString(x,y,esrStr,drawAtoms,false);
            setColor(theColor);
            setTextSize(mpLabelSize);
            }

        if (hydrogensToAdd == 0 && unpairedElectrons == 0)
			return;

		float hindrance[] = new float[4];
		for (int i=0; i<mMol.getAllConnAtoms(atom); i++) {
			int bnd = mMol.getConnBond(atom,i);
			for (int j=0; j<2; j++) {
				if (mMol.getBondAtom(j,bnd) == atom) {
					float theAngle = mMol.getBondAngle(mMol.getBondAtom(j,bnd),mMol.getBondAtom(1-j,bnd));
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
		if (isoStr != null || infoStr != null) hindrance[3] += 10;
								// make sure not to add H's to the left of the
								// hetero atom where they would mess up the atom mass

		String hNoStr = "";
		if (hydrogensToAdd != 0) {
			float hydrogenWidth = getStringWidth("H");
			float hNoWidth = 0.0f;
			if (hydrogensToAdd == -1) {
				hNoStr = "n";
				setTextSize((mpLabelSize*2+1)/3);
				hNoWidth = getStringWidth(hNoStr);
				}
			else if (hydrogensToAdd > 1) {
				hNoStr = String.valueOf(hydrogensToAdd);
				setTextSize((mpLabelSize*2+1)/3);
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
					chay = getAtomY(atom) - getTextSize();
					}
				else {
					hindrance[2] += 10;
					chay = getAtomY(atom) + getTextSize();
					}
				}

			if (hNoWidth > 0) {
				x = chax + ((hydrogenWidth + hNoWidth) / 2.0f);
				y = chay + ((getTextSize()*4+4)/8);
				mpDrawString(x,y,hNoStr,drawAtoms,true);
				setTextSize(mpLabelSize);
				}

			mpDrawString(chax,chay,"H",drawAtoms,true);
			}

		int bestSide = 0;
		if (unpairedElectrons != 0) {
			float minHindrance = 50.0f;
			float counterHindrance = 0.0f;
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
				mpDrawDot(chax,chay,atom,drawAtoms);
				return;
				}

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

			mpDrawDot(chax,chay,atom,drawAtoms);
			mpDrawDot(chax + xdiff,chay + ydiff,atom,drawAtoms);
			}
		}

    private float mpGetFreeSpaceAngle(int atom) {
            // returns the angle from the given atom that is furthest away
            // from any bond and stereo label
        float[] angle = new float[mMol.getAllConnAtoms(atom)];
        for (int i=0; i<mMol.getAllConnAtoms(atom); i++)
            angle[i] = mMol.getBondAngle(atom, mMol.getConnAtom(atom, i));
        Arrays.sort(angle);
        float maxMean = mpGetMeanAngle(angle, 0);
        float maxVal = mpGetAngleESRScore(angle, 0, maxMean);
        for (int i=1; i<angle.length; i++) {
            float mean = mpGetMeanAngle(angle, i);
            float val = mpGetAngleESRScore(angle, i, mean);
            if (maxVal < val) {
                maxVal = val;
                maxMean = mean;
                }
            }
        return maxMean;
        }

    private float mpGetAngleESRScore(float[] angleList, int index, float meanAngle) {
            // initial score is the angle difference between associated bonds
        float score = (index == 0) ? 
                2.0f * (float)Math.PI + angleList[0] - angleList[angleList.length-1]
              : angleList[index] - angleList[index-1];

            // subtract a penalty for colliding with the CIP label
        if (meanAngle > -Math.PI*2/3 && meanAngle < Math.PI/3.0)
            score -= 2*Math.cos(meanAngle + Math.PI/6.0);
        else
            // add a small score for angles on the top right
            score -= 0.5*Math.cos(meanAngle + Math.PI/6.0);

            return score;
        }
    
    private float mpGetMeanAngle(float[] angle, int index) {
        if (index > 0)
            return (angle[index] + angle[index-1]) / 2.0f;

        float mean = (float)Math.PI + (angle[0] + angle[angle.length-1]) / 2.0f;
        return (mean > Math.PI) ? mean - 2.0f * (float)Math.PI : mean;
        }

    private String append(String a, String b) {
        return (a == null) ? b : a+","+b;
        }

	private void mpDrawString(float x,float y,String str,
							  boolean drawAtom,boolean withTabu) {
		if (withTabu) {
			float strWidth,xdiff,ydiff;

			strWidth = getStringWidth(str);
			xdiff = strWidth / 2 + getTextSize() / 8;
			ydiff = getTextSize() / 2;
			if (str == "+" || str == "-")
				ydiff = ydiff * 2 / 3;

			mpTabuZone.add(new Rectangle2D.Float(x-xdiff, y-ydiff, 2*xdiff, 2*ydiff));
			}

		if (drawAtom)
			drawString(str,x,y);
		}


	private void mpDrawDot(float x,float y,int atm,boolean drawDot) {
		mpTabuZone.add(new Rectangle2D.Float(x-mpDotDiameter, y-mpDotDiameter,
											 2*mpDotDiameter, 2*mpDotDiameter));

		if (drawDot) {
			mpDot.add(new DepictorDot(x, y, isHighlightedAtom(atm) ? COLOR_HILITE_BOND_FG : mAtomColor[atm]));
			}
		}


	private void mpDrawAllDots() {
		for (DepictorDot dot:mpDot) {
			setColor(dot.color);
			drawDot(dot.x, dot.y);
			}
		setColor(mDefaultColor);
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
	                setTextSize((mpLabelSize*2+1)/3);
	                textSizeChanged = true;
	                }
	            float x = (getAtomX(atom1) + getAtomX(atom2)) / 2f;
	            float y = (getAtomY(atom1) + getAtomY(atom2)) / 2f;
	            float dx = getAtomX(atom2) - getAtomX(atom1);
	            float dy = getAtomY(atom2) - getAtomY(atom1);
	            float d = (float)Math.sqrt(dx*dx+dy*dy);
	            float hw = 0.60f * getStringWidth(label);	// slightly larger than 0.5f to increase label distance from bond
	            float hh = 0.55f * getTextSize();
	            if (d != 0f) {
	            	if (dx > 0)
	            		mpDrawString(x+hw*dy/d, y-hh*dx/d, label, true, true);
	            	else
	            		mpDrawString(x-hw*dy/d, y+hh*dx/d, label, true, true);
	            	}
            	}
            }

        if (textSizeChanged)
            setTextSize(mpLabelSize);
        }


    private void drawLine(DepictorLine theLine, int atom1, int atom2) {
    	if (mMol.isBondForegroundHilited(mMol.getBond(atom1, atom2))) {
			setColor(COLOR_HILITE_BOND_FG);
			drawBlackLine(theLine);
			setColor(mDefaultColor);
    		}
    	else if (mAtomColor[atom1] != mAtomColor[atom2]) {
    		drawColorLine(theLine, atom1, atom2);
    		}
    	else if (mAtomColor[atom1] != Molecule.cAtomColorBlack) {
			setColor(mAtomColor[atom1]);
			drawBlackLine(theLine);
			setColor(mDefaultColor);
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
			setColor(mAtomColor[atm1]);
			drawBlackLine(line1);
			}
		if (mpProperLine(line2)) {
			setColor(mAtomColor[atm2]);
			drawBlackLine(line2);
			}
		setColor(mDefaultColor);
		}


	private void drawDashedLine(DepictorLine theLine, int atom1, int atom2) {
		float xinc = (theLine.x2 - theLine.x1) / 10;
		float yinc = (theLine.y2 - theLine.y1) / 10;

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

		setColor(color1);
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

		setColor(color2);
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

		setColor(mDefaultColor);
		}


	private void drawWedge(DepictorLine theWedge,int atom1, int atom2) {
		float p1x[],p1y[],p2x[],p2y[];
		float xdiff,ydiff;

		xdiff = (theWedge.y1 - theWedge.y2) / 9;
		ydiff = (theWedge.x2 - theWedge.x1) / 9;
		p1x = new float[3];
		p1y = new float[3];
		p2x = new float[4];
		p2y = new float[4];
		p1x[0] = theWedge.x1;
		p1y[0] = theWedge.y1;
		p2x[2] = (theWedge.x2 + xdiff);
		p2y[2] = (theWedge.y2 + ydiff);
		p2x[3] = (theWedge.x2 - xdiff);
		p2y[3] = (theWedge.y2 - ydiff);
		p1x[1] = (p1x[0] + p2x[2]) / 2;
		p1y[1] = (p1y[0] + p2y[2]) / 2;
		p1x[2] = (p1x[0] + p2x[3]) / 2;
		p1y[2] = (p1y[0] + p2y[3]) / 2;
		p2x[0] = p1x[2];
		p2y[0] = p1y[2];
		p2x[1] = p1x[1];
		p2y[1] = p1y[1];

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

        setColor(color1);
		drawPolygon(p1x,p1y,3);
		setColor(color2);
		drawPolygon(p2x,p2y,4);
		setColor(mDefaultColor);
		}


	protected void drawDot(float x, float y) {
		fillCircle(x-mpDotDiameter/2, y-mpDotDiameter/2, mpDotDiameter);
		}


	private void setRGBColor(Color rgbColor) {
		if (mOverruleForeground != null) {
			if (mCurrentColor != COLOR_OVERRULED) {
				mCurrentColor = COLOR_OVERRULED;
			    setColor(mOverruleForeground);
				}
			return;
			}

		mCurrentColor = RGB_COLOR;
		setColor(rgbColor);
		}


	public void setColor(int theColor) {
		if (theColor != COLOR_HILITE_BOND_BG && mOverruleForeground != null)
			theColor = COLOR_OVERRULED;

		if (theColor == mCurrentColor)
			return;

		mCurrentColor = theColor;

		switch (theColor) {
		case COLOR_OVERRULED:
		    setColor(mOverruleForeground);
			break;
		case COLOR_HILITE_BOND_BG:
		    setColor(mBondBGHiliteColor);
			break;
		case COLOR_HILITE_BOND_FG:
		    setColor(mBondFGHiliteColor);
			break;
		case Molecule.cAtomColorBlue:
		    setColor(Color.blue);
			break;
		case Molecule.cAtomColorRed:
		    setColor(Color.red);
			break;
		case Molecule.cAtomColorMagenta:
		    setColor(Color.magenta);
			break;
		case Molecule.cAtomColorGreen:
		    setColor(Color.green);
			break;
		case Molecule.cAtomColorOrange:
		    setColor(Color.orange);
			break;
		case Molecule.cAtomColorDarkGreen:
		    setColor(new Color(0,160,0));
			break;
		case Molecule.cAtomColorDarkRed:
		    setColor(new Color(160,0,0));
			break;
		case cColorGray:
		    setColor(Color.gray);
			break;
		default:
		    setColor(Color.black);
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
		if ((mDisplayMode & cDModeSuppressESR) != 0)
			return mAtomColor[atom];

		int alleneCenter = mMol.findAlleneCenterAtom(atom);
		if (alleneCenter != -1)
			atom = alleneCenter;

		int esrInfo = getESRTypeToDisplayAt(atom);
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

/*	private void drawBlackArc(MyRect theArc,float arcAngle,boolean inverted) {
		float xdif = (float)(theArc.x2 - theArc.x1);
		float ydif = (float)(theArc.y2 - theArc.y1);
		float length = Math.sqrt(xdif*xdif + ydif*ydif);
		float theAngle = (inverted) ? (Math.PI - cArcAngle) / 2
									 : (cArcAngle - Math.PI) / 2;
		float radius = (length / 2) / Math.cos(theAngle);
		float centerX = (float)theArc.x1 + radius * Math.sin(arcAngle + theAngle);
		float centerY = (float)theArc.y1 + radius * Math.cos(arcAngle + theAngle);

		mG.drawArc((int)(centerX - radius + 0.5),(int)(centerY - radius + 0.5),
				   (int)(2 * radius + 0.5),(int)(2 * radius + 0.5),
				   (int)((arcAngle - theAngle) * 180 / Math.PI) - 90,
				   (inverted) ? (int)(-cArcAngle * 180 / Math.PI + 0.5)
				              : (int)(cArcAngle * 180 / Math.PI + 0.5));
		}	*/


	protected abstract void drawBlackLine(DepictorLine theLine);
    protected abstract void drawDottedLine(DepictorLine theLine);
	protected abstract void drawPolygon(float[] x, float[] y, int count);
	protected abstract void drawString(String theString,float x,float y);
	protected abstract void fillCircle(float x, float y, float r);
	protected abstract float getLineWidth();
	protected abstract float getStringWidth(String theString);
    protected abstract int getTextSize();
	protected abstract void setTextSize(int theSize);
	protected abstract void setLineWidth(float lineWidth);
	protected abstract void setColor(Color theColor);

    public static class DepictorDot {
        public float x,y;
        public int color;

        DepictorDot(float x, float y, int color) {
            this.x = x;
            this.y = y;
            this.color = color;
        }
    }


    public static class DepictorLine {
        public float x1;
        public float y1;
        public float x2;
        public float y2;

        public DepictorLine(float x1, float y1, float x2, float y2)
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


