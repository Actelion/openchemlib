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
*/

package com.actelion.research.chem;

import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionArrow;

import java.awt.*;
import java.awt.geom.Rectangle2D;



public class ExtendedDepictor {
    private StereoMolecule[]		mMolecule,mCatalyst;
    private AbstractDepictor[]		mDepictor,mCatalystDepictor;
    private DrawingObjectList		mDrawingObjectList;
    private int						mDisplayMode,mReactantOrCoreCount;
    private boolean					mUseGraphics2D,mReactionLayoutNeeded,mIsMarkushStructure;
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
			mCatalyst = new StereoMolecule[reaction.getCatalysts()];
	        for (int i=0; i<reaction.getCatalysts(); i++)
                mCatalyst[i] = reaction.getCatalyst(i);
            mReactionLayoutNeeded = layoutReaction;
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
		if (mCatalyst != null) {
			mCatalystDepictor = new AbstractDepictor[mCatalyst.length];
			for (int i=0; i<mCatalyst.length; i++) {
				if (mUseGraphics2D)
					mCatalystDepictor[i] = new Depictor2D(mCatalyst[i]);
				else
					mCatalystDepictor[i] = new Depictor(mCatalyst[i]);
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

    public StereoMolecule getMolecule(int i) {
        return mMolecule[i];
        }

    public AbstractDepictor getMoleculeDepictor(int i) {
        return mDepictor[i];
        }

    public void setForegroundColor(Color foreGround, Color background) {
        if (mDepictor != null)
            for (AbstractDepictor d:mDepictor)
                d.setForegroundColor(foreGround, background);

		if (mCatalystDepictor != null)
			for (AbstractDepictor d:mCatalystDepictor)
				d.setForegroundColor(foreGround, background);
        }

	public void setOverruleColor(Color foreGround, Color background) {
		if (mDepictor != null)
			for (AbstractDepictor d:mDepictor)
				d.setOverruleColor(foreGround, background);

		if (mCatalystDepictor != null)
			for (AbstractDepictor d:mCatalystDepictor)
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
            double averageBondLength = calculateAverageBondLength();
            g.setColor(mFragmentNoColor);
            g.setFont(g.getFont().deriveFont(Font.BOLD, (int)(1.6*averageBondLength)));
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
		if (mCatalystDepictor != null) {
			for (int i=0; i<mCatalystDepictor.length; i++) {
//				mCatalystDepictor[i].setDisplayMode(mDisplayMode);
				mCatalystDepictor[i].paint(g);
/*
Rectangle2D.Float r = mCatalystDepictor[i].getBoundingRect();
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

    public DepictorTransformation updateCoords(Object g, Rectangle2D.Double viewRect, int mode) {
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

			if (mCatalystDepictor != null)
				for (int i=0; i<mCatalystDepictor.length; i++)
					mCatalystDepictor[i].getTransformation().clear();

			DepictorTransformation t = mTransformation;
            mTransformation = new DepictorTransformation();
            return t;
            }
        }

    public DepictorTransformation validateView(Object g, Rectangle2D.Double viewRect, int mode) {
    // returns incremental transformation that moves/scales already transformed molecules/objects into viewRect
        if (mReactionLayoutNeeded)
            layoutReaction(g);

        Rectangle2D.Double boundingRect = null;
        if (mDepictor != null) {
            for (int i=0; i<mDepictor.length; i++) {
                mDepictor[i].validateView(g, null, 0);
                boundingRect = (boundingRect == null) ? mDepictor[i].getBoundingRect()
                            : (Rectangle2D.Double)boundingRect.createUnion(mDepictor[i].getBoundingRect());
                }
            }
		if (mCatalystDepictor != null) {
			for (int i=0; i<mCatalystDepictor.length; i++) {
				mCatalystDepictor[i].validateView(g, null, 0);
				boundingRect = (boundingRect == null) ? mCatalystDepictor[i].getBoundingRect()
						: (Rectangle2D.Double)boundingRect.createUnion(mCatalystDepictor[i].getBoundingRect());
				}
			}
        if (mDrawingObjectList != null) {
            for (int i=0; i<mDrawingObjectList.size(); i++) {
                Rectangle2D.Double objectBounds = mDrawingObjectList.get(i).getBoundingRect();
                mTransformation.applyTo(objectBounds);
                boundingRect = (boundingRect == null) ? objectBounds
                        : (Rectangle2D.Double)boundingRect.createUnion(objectBounds);
                }
            }

        if (boundingRect == null)
            return null;

        double avbl = calculateAverageBondLength();

        DepictorTransformation t = new DepictorTransformation(boundingRect, viewRect, avbl, mode);

        if (!t.isVoidTransformation()) {
            t.applyTo(mTransformation);

            if (mDepictor != null)
                for (int i=0; i<mDepictor.length; i++)
                    mDepictor[i].applyTransformation(t);

			if (mCatalystDepictor != null)
				for (int i=0; i<mCatalystDepictor.length; i++)
					mCatalystDepictor[i].applyTransformation(t);

			return t;
            }

        return null;
        }

    private double calculateAverageBondLength() {
    	float averageBondLength = 0.0f;
        int bondCount = 0;
        if (mMolecule != null) {
            for (int i=0; i<mMolecule.length; i++) {
                if (mMolecule[i].getAllAtoms() != 0) {
                    if (mMolecule[i].getAllBonds() != 0) {
                        averageBondLength += mDepictor[i].getTransformation().getScaling()
                                * mMolecule[i].getAllBonds() * mMolecule[i].getAverageBondLength();
                        bondCount += mMolecule[i].getAllBonds();
                        }
                    else {
						averageBondLength += mDepictor[i].getTransformation().getScaling()
								* mMolecule[i].getAverageBondLength();
						bondCount ++;
                        }
                    }
                }
            }
        return (bondCount == 0) ? AbstractDepictor.cOptAvBondLen
                                : mTransformation.getScaling() * averageBondLength / bondCount;
        }

    private void layoutReaction(Object g) {
        Rectangle2D.Double[] boundingRect = new Rectangle2D.Double[mMolecule.length];
        double totalWidth = 0.0;
        double totalHeight = 0.0;
        for (int i=0; i<mMolecule.length; i++) {
            mDepictor[i].validateView(g, null, AbstractDepictor.cModeInflateToMaxAVBL);
            boundingRect[i] = mDepictor[i].getBoundingRect();
            totalWidth += boundingRect[i].width;
            totalHeight = Math.max(totalHeight, boundingRect[i].height);
            }

        final double catalystScale = 0.7;
		double catalystSpacing = 0.5 * AbstractDepictor.cOptAvBondLen;
		Rectangle2D.Double[] catalystBoundingRect = new Rectangle2D.Double[mCatalyst.length];
		double totalCatalystWidth = 0.0;
		double totalCatalystHeight = 0.0;
		for (int i=0; i<mCatalyst.length; i++) {
			mCatalystDepictor[i].validateView(g, null, AbstractDepictor.cModeInflateToMaxAVBL+(int)(catalystScale*AbstractDepictor.cOptAvBondLen));
			catalystBoundingRect[i] = mCatalystDepictor[i].getBoundingRect();
			totalCatalystWidth = Math.max(totalCatalystWidth, catalystBoundingRect[i].width);
			totalCatalystHeight += catalystBoundingRect[i].height + catalystSpacing;
			}

		double spacing = 1.5 * AbstractDepictor.cOptAvBondLen;
        double arrowWidth = Math.max(2 * AbstractDepictor.cOptAvBondLen, totalCatalystWidth + AbstractDepictor.cOptAvBondLen);

		totalHeight = Math.max(totalHeight, AbstractDepictor.cOptAvBondLen + 2 * totalCatalystHeight);

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

        double rawX = 0.5 * spacing;
        for (int i=0; i<mMolecule.length; i++) {
            if (i == mReactantOrCoreCount) {
                ((ReactionArrow)mDrawingObjectList.get(arrow)).setCoordinates(
                        rawX-spacing/2, totalHeight/2, rawX-spacing/2+arrowWidth, totalHeight/2);

                double catX = rawX + 0.5 * (AbstractDepictor.cOptAvBondLen - spacing);
                double catY = 0.5 * (totalHeight - catalystSpacing) - totalCatalystHeight;
                for (int j=0; j<mCatalyst.length; j++) {
					double dx = catX + 0.5 * (totalCatalystWidth - catalystBoundingRect[j].width) - catalystBoundingRect[j].x;
					double dy = catY - catalystBoundingRect[j].y;
					mCatalystDepictor[j].applyTransformation(new DepictorTransformation(1.0, dx, dy));

					catY += catalystSpacing + catalystBoundingRect[j].height;
					}

                rawX += arrowWidth;
                }

            double dx = rawX - boundingRect[i].x;
            double dy = 0.5 * (totalHeight - boundingRect[i].height) - boundingRect[i].y;
            mDepictor[i].applyTransformation(new DepictorTransformation(1.0, dx, dy));

            rawX += spacing + boundingRect[i].width;
            }

        mReactionLayoutNeeded = false;
        }
    }
