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

import java.awt.*;
import java.awt.geom.*;

import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionArrow;




public class ExtendedDepictor {
    private StereoMolecule[]		mMolecule;
    private AbstractDepictor[]		mDepictor;
    private DrawingObjectList		mDrawingObjectList;
    private int						mDisplayMode,mReactantOrCoreCount;
    private boolean					mUseGraphics2D,mDoLayoutMolecules,mIsMarkushStructure;
    private DepictorTransformation	mTransformation;
    private Color                   mFragmentNoColor;

    public ExtendedDepictor(StereoMolecule mol, DrawingObjectList drawingObjectList, boolean useGraphics2D) {
        if (mol != null) {
            mMolecule = new StereoMolecule[1];
            mMolecule[0] = mol;
            }
        mIsMarkushStructure = false;
        mDrawingObjectList = drawingObjectList;
        mUseGraphics2D = useGraphics2D;
        mReactantOrCoreCount = -1;
        initialize();
        }

    public ExtendedDepictor(StereoMolecule[] mol, DrawingObjectList drawingObjectList, boolean useGraphics2D) {
        mMolecule = mol;
        mIsMarkushStructure = false;
        mDrawingObjectList = drawingObjectList;
        mUseGraphics2D = useGraphics2D;
        mReactantOrCoreCount = -1;
        initialize();
        }

    /**
     * Use this constructor for markush structures. The first fragments in the list
     * are the Markush core structures (typically only one), decorated with R1,R2,R3,...
     * The remaining fragments need to contain one atom with atomicNo=0 each, that
     * indicates the attachment point. They also may contain Rn atoms.
     * Any of the fragments may contain query features.
     * @param mol
     * @param markushCoreCount
     * @param drawingObjectList
     * @param useGraphics2D
     */
    public ExtendedDepictor(StereoMolecule[] mol, int markushCoreCount, DrawingObjectList drawingObjectList, boolean useGraphics2D) {
        mMolecule = mol;
        mIsMarkushStructure = true;
        mDrawingObjectList = drawingObjectList;
        mUseGraphics2D = useGraphics2D;
        mReactantOrCoreCount = markushCoreCount;
        initialize();
        }

    public ExtendedDepictor(Reaction reaction, DrawingObjectList drawingObjectList, boolean layoutReaction, boolean useGraphics2D) {
        if (reaction != null) {
            mMolecule = new StereoMolecule[reaction.getMolecules()];
            for (int i=0; i<reaction.getMolecules(); i++)
                mMolecule[i] = reaction.getMolecule(i);
            mReactantOrCoreCount = reaction.getReactants();
            mDoLayoutMolecules = layoutReaction;
            }
        mIsMarkushStructure = false;
        mDrawingObjectList = drawingObjectList;
        mUseGraphics2D = useGraphics2D;
        initialize();
        }

    private void initialize() {
            // for reactions and sets of molecules the availability of coordinates
            // is mandatory. However, every individual molecule may have its first
            // atom at coords 0.0/0.0, e.g. if they are encoded idcodes or if
        mTransformation = new DepictorTransformation();
        if (mMolecule != null) {
            mDepictor = new AbstractDepictor[mMolecule.length];
            for (int i=0; i<mMolecule.length; i++) {
                if (mUseGraphics2D)
                    mDepictor[i] = new Depictor2D(mMolecule[i]);
                else
                    mDepictor[i] = new Depictor(mMolecule[i]);
                }
            }
        }

    public void setDisplayMode(int displayMode) {
        mDisplayMode = displayMode;
        }

    public void setFragmentNoColor(Color c) {
        // use setFragmentNoColor(null) if you don't want fragment numbers to be shown
        mFragmentNoColor = c;
        }

    public int getMoleculeCount() {
        return mMolecule.length;
        }

    public ExtendedMolecule getMolecule(int i) {
        return mMolecule[i];
        }

    public AbstractDepictor getMoleculeDepictor(int i) {
        return mDepictor[i];
        }

	public void setOverruleColor(Color foreGround, Color background) {
		if (mDepictor != null)
			for (AbstractDepictor d:mDepictor)
				d.setOverruleColor(foreGround, background);
		}

    public void paint(Graphics g)
    {
        Color saveColor = g.getColor();
        Font saveFont = g.getFont();
        try {
            paintFragmentNumbers(g);
            paintStructures(g);
            paintDrawingObjects(g);
        } finally {
            g.setColor(saveColor);
            g.setFont(saveFont);
        }
    }

    public void paintFragmentNumbers(Graphics g) {
        if (mFragmentNoColor != null && mMolecule != null) {
        	float averageBondLength = calculateAverageBondLength();
            g.setColor(mFragmentNoColor);
            g.setFont(new Font("Helvetica",Font.BOLD, (int)(1.6*averageBondLength)));
            for (int i=0; i<mMolecule.length; i++) {
                if (mMolecule[i].getAllAtoms() != 0) {
                    Point cog = new Point();
                    for (int atom=0; atom<mMolecule[i].getAllAtoms(); atom++) {
                        cog.x += mMolecule[i].getAtomX(atom);
                        cog.y += mMolecule[i].getAtomY(atom);
                        }
                    cog.x /= mMolecule[i].getAllAtoms();
                    cog.y /= mMolecule[i].getAllAtoms();
                    cog.x = (int)mDepictor[i].getTransformation().transformX(cog.x);
                    cog.y = (int)mDepictor[i].getTransformation().transformY(cog.y);

                    String str = (mReactantOrCoreCount == -1) ? "F"+(i+1)
                                 : (i < mReactantOrCoreCount) ? ""+(char)('A'+i)
                                 : (mIsMarkushStructure) ? "R"+(i+1-mReactantOrCoreCount)
                                 :                        "P"+(i+1-mReactantOrCoreCount);
                    int width = g.getFontMetrics().stringWidth(str);
                    g.drawString(str, cog.x-width/2, cog.y+(int)(0.3*g.getFontMetrics().getHeight()));

                    /* this would require an mDepictor[i].validateView(...)
                    if (mIsMarkushStructure && i>= mReactantOrCoreCount) {
                        Rectangle2D.Float r = mDepictor[i].getBoundingRect();
                        if (r != null) {
                            g.drawRect((int)r.x-8, (int)r.y-8, (int)r.width+16, (int)r.height+16);
                            g.drawRect((int)r.x-7, (int)r.y-7, (int)r.width+14, (int)r.height+14);
                            }
                        } */
                    }
                }
            }
        }

    public void paintStructures(Object g) {
        if (mDepictor != null) {
            for (int i=0; i<mDepictor.length; i++) {
                mDepictor[i].setDisplayMode(mDisplayMode);
                mDepictor[i].paint(g);
/*
Rectangle2D.Float r = mDepictor[i].getBoundingRect();
if (r != null) {
g.setColor(Color.magenta);
g.drawRect((int)r.x, (int)r.y, (int)r.width, (int)r.height);
}*/
                }
            }
        }

    public void paintDrawingObjects(Graphics g) {
        if (mDrawingObjectList != null) {
            for (int i=0; i<mDrawingObjectList.size(); i++) {
                AbstractDrawingObject object = (AbstractDrawingObject)mDrawingObjectList.get(i);
                if (mUseGraphics2D)
                    object.draw2D((Graphics2D)g, mTransformation);
                else
                    object.draw(g, mTransformation);
/*
Rectangle2D.Float r = object.getBoundingRect();
mTransformation.applyTo(r);
g.setColor(Color.magenta);
g.drawRect((int)r.x, (int)r.y, (int)r.width, (int)r.height);*/
                }
            }
        }

    public DepictorTransformation updateCoords(Object g, Rectangle2D.Float viewRect, int mode) {
    // returns full transformation that moves/scales original molecules/objects into viewRect
        validateView(g, viewRect, mode);

        if (mTransformation.isVoidTransformation()) {
            return null;
            }
        else {
            if (mMolecule != null)
                for (int i=0; i<mMolecule.length; i++)
                    mTransformation.applyTo(mMolecule[i]);

            if (mDrawingObjectList != null)
                for (int i=0; i<mDrawingObjectList.size(); i++)
                    mTransformation.applyTo((AbstractDrawingObject)mDrawingObjectList.get(i));

            if (mDepictor != null)
                for (int i=0; i<mDepictor.length; i++)
                    mDepictor[i].getTransformation().clear();

            DepictorTransformation t = mTransformation;
            mTransformation = new DepictorTransformation();
            return t;
            }
        }

    public DepictorTransformation validateView(Object g, Rectangle2D.Float viewRect, int mode) {
    // returns incremental transformation that moves/scales already transformed molecules/objects into viewRect
        if (mDoLayoutMolecules)
            doLayoutMolecules(g);

        Rectangle2D.Float boundingRect = null;
        if (mDepictor != null) {
            for (int i=0; i<mDepictor.length; i++) {
                mDepictor[i].validateView(g, null, 0);
                boundingRect = (boundingRect == null) ? mDepictor[i].getBoundingRect()
                            : (Rectangle2D.Float)boundingRect.createUnion(mDepictor[i].getBoundingRect());
                }
            }
        if (mDrawingObjectList != null) {
            for (int i=0; i<mDrawingObjectList.size(); i++) {
                Rectangle2D.Float objectBounds = ((AbstractDrawingObject)mDrawingObjectList.get(i)).getBoundingRect();
                mTransformation.applyTo(objectBounds);
                boundingRect = (boundingRect == null) ? objectBounds
                        : (Rectangle2D.Float)boundingRect.createUnion(objectBounds);
                }
            }

        if (boundingRect == null)
            return null;

        float avbl = calculateAverageBondLength();

        DepictorTransformation t = new DepictorTransformation(boundingRect, viewRect, avbl, mode);

        if (!t.isVoidTransformation()) {
            t.applyTo(mTransformation);

            if (mDepictor != null)
                for (int i=0; i<mDepictor.length; i++)
                    mDepictor[i].applyTransformation(t);

            return t;
            }

        return null;
        }

    private float calculateAverageBondLength() {
    	float averageBondLength = 0.0f;
        int bondCount = 0;
        if (mMolecule != null) {
            for (int i=0; i<mMolecule.length; i++) {
                if (mMolecule[i].getAllAtoms() != 0) {
                    averageBondLength += mDepictor[i].getTransformation().getScaling()
                                       * mMolecule[i].getAllBonds() * mMolecule[i].getAverageBondLength();
                    bondCount += mMolecule[i].getAllBonds();
                    }
                }
            }
        return (bondCount == 0) ? AbstractDepictor.cOptAvBondLen
                                : mTransformation.getScaling() * averageBondLength / bondCount;
        }

    private void doLayoutMolecules(Object g) {
        Rectangle2D.Float[] boundingRect = new Rectangle2D.Float[mMolecule.length];
        float totalWidth = 0.0f;
        float totalHeight = 0.0f;
        for (int i=0; i<mMolecule.length; i++) {
            mDepictor[i].validateView(g, null, AbstractDepictor.cModeInflateToMaxAVBL);
            boundingRect[i] = mDepictor[i].getBoundingRect();
            totalWidth += boundingRect[i].width;
            totalHeight = Math.max(totalHeight, boundingRect[i].height);
            }

        float spacing = 1.5f * AbstractDepictor.cOptAvBondLen;
        float arrowWidth = 2f * AbstractDepictor.cOptAvBondLen;

        int arrow = -1;
        if (mDrawingObjectList == null) {
            mDrawingObjectList = new DrawingObjectList();
            mDrawingObjectList.add(new ReactionArrow());
            arrow = 0;
            }
        else {
            for (int i=0; i<mDrawingObjectList.size(); i++) {
                if (mDrawingObjectList.get(i) instanceof ReactionArrow) {
                    arrow = i;
                    break;
                    }
                }
            if (arrow == -1) {
                arrow = mDrawingObjectList.size();
                mDrawingObjectList.add(new ReactionArrow());
                }
            }

        float rawX = 0.5f * spacing;
        for (int i=0; i<mMolecule.length; i++) {
            if (i == mReactantOrCoreCount) {
                ((ReactionArrow)mDrawingObjectList.get(arrow)).setCoordinates(
                        rawX-spacing/2, totalHeight/2, rawX-spacing/2+arrowWidth, totalHeight/2);
                rawX += arrowWidth;
                }

            float dx = rawX - boundingRect[i].x;
            float dy = 0.5f * (totalHeight - boundingRect[i].height) - boundingRect[i].y;
            mDepictor[i].applyTransformation(new DepictorTransformation(1.0f, dx, dy));

            rawX += spacing + boundingRect[i].width;
            }

        mDoLayoutMolecules = false;
        }
    }
