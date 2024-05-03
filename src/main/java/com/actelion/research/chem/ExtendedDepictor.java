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

import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionArrow;
import com.actelion.research.gui.generic.GenericDepictor;
import com.actelion.research.gui.generic.GenericDrawContext;
import com.actelion.research.gui.generic.GenericRectangle;

import java.awt.*;
import java.util.PriorityQueue;


public class ExtendedDepictor {
	public static final int TYPE_MOLECULES = 0;
	public static final int TYPE_REACTION = 1;
	public static final int TYPE_MARKUSH = 2;

    private StereoMolecule[]		mMolecule,mCatalyst;
    private Reaction				mReaction;
    private GenericDepictor[]		mDepictor,mCatalystDepictor;
    private DrawingObjectList		mDrawingObjectList;
    private int						mDisplayMode,mReactantCount,mMarkushCoreCount,mChemistryType;
    private boolean					mReactionLayoutNeeded;
    private DepictorTransformation	mTransformation;
    private int                     mFragmentNoColor,mDefaultAVBL;

    public ExtendedDepictor(StereoMolecule mol, DrawingObjectList drawingObjectList) {
        if (mol != null) {
            mMolecule = new StereoMolecule[1];
            mMolecule[0] = mol;
            }
		mChemistryType = TYPE_MOLECULES;
        mDrawingObjectList = drawingObjectList;
        initialize();
        }

    public ExtendedDepictor(StereoMolecule[] mol, DrawingObjectList drawingObjectList) {
        mMolecule = mol;
		mChemistryType = TYPE_MOLECULES;
        mDrawingObjectList = drawingObjectList;
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
     */
    public ExtendedDepictor(StereoMolecule[] mol, int markushCoreCount, DrawingObjectList drawingObjectList) {
        mMolecule = mol;
		mChemistryType = TYPE_MARKUSH;
        mDrawingObjectList = drawingObjectList;
        mMarkushCoreCount = markushCoreCount;
        initialize();
        }

    public ExtendedDepictor(Reaction reaction, DrawingObjectList drawingObjectList, boolean layoutReaction) {
		mReaction = reaction;
        if (reaction != null) {
            mMolecule = new StereoMolecule[reaction.getMolecules()];
            for (int i=0; i<reaction.getMolecules(); i++)
                mMolecule[i] = reaction.getMolecule(i);
            mReactantCount = reaction.getReactants();
			mCatalyst = new StereoMolecule[reaction.getCatalysts()];
	        for (int i=0; i<reaction.getCatalysts(); i++)
                mCatalyst[i] = reaction.getCatalyst(i);
            mReactionLayoutNeeded = layoutReaction;
            }
		mChemistryType = TYPE_REACTION;
        mDrawingObjectList = drawingObjectList;
        initialize();
        }

    public boolean isFragment() {
    	if (mChemistryType == TYPE_REACTION)
    		return mReaction == null ? false : mReaction.isFragment();

    	if (mMolecule == null)
    		return false;

    	for (StereoMolecule mol:mMolecule)
    		if (mol.isFragment())
    			return true;

     	return false;
		}

    private void initialize() {
            // for reactions and sets of molecules the availability of coordinates
            // is mandatory. However, every individual molecule may have its first
            // atom at coords 0.0/0.0, e.g. if they are encoded idcodes or if
        mTransformation = new DepictorTransformation();
        if (mMolecule != null) {
            mDepictor = new GenericDepictor[mMolecule.length];
            for (int i=0; i<mMolecule.length; i++)
                mDepictor[i] = new GenericDepictor(mMolecule[i]);
            }
		if (mCatalyst != null) {
			mCatalystDepictor = new GenericDepictor[mCatalyst.length];
			for (int i=0; i<mCatalyst.length; i++)
				mCatalystDepictor[i] = new GenericDepictor(mCatalyst[i]);
			}
		mDefaultAVBL = AbstractDepictor.cOptAvBondLen;
        }

    public void setDisplayMode(int displayMode) {
        mDisplayMode = displayMode;
        }

    public void setDefaultAVBL(int avbl) {
    	mDefaultAVBL = avbl;
    }

    public void setFactorTextSize(double factor) {
    	if (mDepictor != null)
	    	for (GenericDepictor d:mDepictor)
	    		d.setFactorTextSize(factor);
        }

    public void setFragmentNoColor(int argb) {
        // use setFragmentNoColor(null) if you don't want fragment numbers to be shown
        mFragmentNoColor = argb;
        }

    public int getMoleculeCount() {
        return mMolecule == null ? 0 : mMolecule.length;
        }

    public StereoMolecule getMolecule(int i) {
        return mMolecule[i];
        }

    public Reaction getReaction() {
    	// as long as we don't change the input reaction, we can just buffer it and return it
		return mReaction;
		}

    public AbstractDepictor getMoleculeDepictor(int i) {
        return mDepictor[i];
        }

    @Deprecated
    // Use rgb version of this method instead
	public void setForegroundColor(Color foreGround, Color background) {
	    setForegroundColor(foreGround.getRGB(), background.getRGB());
		}

	public void setForegroundColor(int foreground, int background) {
        if (mDepictor != null)
            for (GenericDepictor d:mDepictor)
                d.setForegroundColor(foreground, background);

		if (mCatalystDepictor != null)
			for (GenericDepictor d:mCatalystDepictor)
				d.setForegroundColor(foreground, background);
        }

	@Deprecated
	// Use rgb version of this method instead
	public void setOverruleColor(Color foreground, Color background) {
		if (mDepictor != null)
			for (GenericDepictor d:mDepictor)
				d.setOverruleColor(foreground, background);

		if (mCatalystDepictor != null)
			for (GenericDepictor d:mCatalystDepictor)
				d.setOverruleColor(foreground, background);
		}

	public void setOverruleColor(int foreground, int background) {
		if (mDepictor != null)
			for (GenericDepictor d:mDepictor)
				d.setOverruleColor(foreground, background);

		if (mCatalystDepictor != null)
			for (GenericDepictor d:mCatalystDepictor)
				d.setOverruleColor(foreground, background);
	}

	public void paint(GenericDrawContext context) {
        int saveRGB = context.getRGB();
        int fontSize = context.getFontSize();
        try {
            paintFragmentNumbers(context);
            paintStructures(context);
            paintDrawingObjects(context);
        } finally {
            context.setRGB(saveRGB);
            context.setFont(fontSize, false, false);
        }
    }

    public void paintFragmentNumbers(GenericDrawContext context) {
        if (mFragmentNoColor != 0 && mMolecule != null) {
            double averageBondLength = calculateMedianBondLength() / mTransformation.getScaling();
            context.setRGB(mFragmentNoColor);
            context.setFont((int)(1.6*averageBondLength), true, false);
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


                    String str = (mChemistryType == TYPE_MOLECULES) ? ""+(i+1)
                               : (mChemistryType == TYPE_MARKUSH) ? ((i < mMarkushCoreCount) ? ""+(char)('A'+i) : "R"+(i+1-mMarkushCoreCount))
                               : (mChemistryType == TYPE_REACTION) ? ((i < mReactantCount) ? ""+(char)('A'+i) : "P"+(i+1-mReactantCount)) : "?"+(i+1);
                    context.drawCenteredString(cog.x, cog.y, str);

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

    public void paintStructures(GenericDrawContext context) {
        if (mDepictor != null) {
	        double avbl = calculateMedianBondLength() / mTransformation.getScaling();  // this still contains individual depictor scaling
            for (GenericDepictor d:mDepictor) {
                d.setDisplayMode(mDisplayMode);
				d.setAtomLabelAVBL(avbl);
                d.paint(context);
/*
Rectangle2D.Float r = mDepictor[i].getBoundingRect();
if (r != null) {
g.setColor(Color.magenta);
g.drawRect((int)r.x, (int)r.y, (int)r.width, (int)r.height);
}*/
                }
            }
		if (mCatalystDepictor != null) {
			for (GenericDepictor d:mCatalystDepictor) {
//				d.setDisplayMode(mDisplayMode);
//				d.setFactorTextSize(mFactorTextSize);
				d.paint(context);
/*
Rectangle2D.Float r = mCatalystDepictor[i].getBoundingRect();
if (r != null) {
g.setColor(Color.magenta);
g.drawRect((int)r.x, (int)r.y, (int)r.width, (int)r.height);
}*/
				}
			}
        }

    public void paintDrawingObjects(GenericDrawContext context) {
        if (mDrawingObjectList != null) {
            for (AbstractDrawingObject object:mDrawingObjectList) {
                object.draw(context, mTransformation);
/*
Rectangle2D.Float r = object.getBoundingRect();
mTransformation.applyTo(r);
g.setColor(Color.magenta);
g.drawRect((int)r.x, (int)r.y, (int)r.width, (int)r.height);*/
	            }
            }
        }

	public DepictorTransformation updateCoords(GenericDrawContext context, GenericRectangle viewRect, int mode) {
    // returns full transformation that moves/scales original molecules/objects into viewRect
        validateView(context, viewRect, mode);

		if (mTransformation.isVoidTransformation()) {
            return null;
            }
        else {
            if (mMolecule != null)
                for (StereoMolecule mol:mMolecule)
                    mTransformation.applyTo(mol);

            if (mDrawingObjectList != null)
                for (AbstractDrawingObject o:mDrawingObjectList)
                    mTransformation.applyTo(o);

            if (mDepictor != null)
                for (GenericDepictor d:mDepictor)
                    d.getTransformation().clear();

			if (mCatalystDepictor != null)
				for (GenericDepictor d:mCatalystDepictor)
					d.getTransformation().clear();

			DepictorTransformation t = mTransformation;
            mTransformation = new DepictorTransformation();
            return t;
            }
        }

    public DepictorTransformation validateView(GenericDrawContext context, GenericRectangle viewRect, int mode) {
    // returns incremental transformation that moves/scales already transformed molecules/objects into viewRect
        if (mReactionLayoutNeeded)
            layoutReaction(context);

        GenericRectangle boundingRect = null;
        if (mDepictor != null) {
            for (GenericDepictor d:mDepictor) {
                d.validateView(context, null, 0);
                boundingRect = (boundingRect == null) ? d.getBoundingRect() : boundingRect.union(d.getBoundingRect());
                }
            }
		if (mCatalystDepictor != null) {
			for (GenericDepictor d:mCatalystDepictor) {
				d.validateView(context, null, 0);
				boundingRect = (boundingRect == null) ? d.getBoundingRect() : boundingRect.union(d.getBoundingRect());
				}
			}
        if (mDrawingObjectList != null) {
            for (AbstractDrawingObject o:mDrawingObjectList) {
	            GenericRectangle objectBounds = o.getBoundingRect(context);
                mTransformation.applyTo(objectBounds);
                boundingRect = (boundingRect == null) ? objectBounds : boundingRect.union(objectBounds);
                }
            }

        if (boundingRect == null)
            return null;

        double avbl = calculateMedianBondLength();

        DepictorTransformation t = new DepictorTransformation(boundingRect, viewRect, avbl, mode);

        if (!t.isVoidTransformation()) {
            t.applyTo(mTransformation);

            if (mDepictor != null)
                for (GenericDepictor d:mDepictor)
	                d.applyTransformation(t);

			if (mCatalystDepictor != null)
				for (GenericDepictor d:mCatalystDepictor)
					d.applyTransformation(t);

			return t;
            }

        return null;
        }

	private double calculateMedianBondLength() {
		PriorityQueue<Double> maxHeap = new PriorityQueue<>((a, b) -> (a > b) ? -1 : (a < b) ? 1 : 0);
		PriorityQueue<Double> minHeap = new PriorityQueue<>();

		if (mMolecule != null) {
			for (int i=0; i<mMolecule.length; i++) {
				for (int bond=0; bond<mMolecule[i].getAllBonds(); bond++) {
					maxHeap.offer(mDepictor[i].getTransformation().getScaling() * mMolecule[i].getBondLength(bond));
					minHeap.offer(maxHeap.poll());
					if (maxHeap.size() < minHeap.size())
						maxHeap.offer(minHeap.poll());
					}
				}
			}

		int bondCount = maxHeap.size() + minHeap.size();
		return (bondCount == 0) ? calculatePseudoBondLengthFromBounds() : mTransformation.getScaling() *
			((bondCount % 2 == 0) ? (maxHeap.peek() + minHeap.peek()) / 2.0 : maxHeap.peek());
		}

	private double calculatePseudoBondLengthFromBounds() {
		double x1 = Double.MAX_VALUE;
		double x2 = -Double.MAX_VALUE;
		double y1 = Double.MAX_VALUE;
		double y2 = -Double.MAX_VALUE;
		int count = 0;
		if (mMolecule != null) {
			for (int i=0; i<mMolecule.length; i++) {
				for (int atom=0; atom<mMolecule[i].getAllAtoms(); atom++) {
					double x = mDepictor[i].getTransformation().transformX(mMolecule[i].getCoordinates(atom).x);
					double y = mDepictor[i].getTransformation().transformY(mMolecule[i].getCoordinates(atom).y);
					x1 = Math.min(x1, x);
					x2 = Math.max(x2, x);
					y1 = Math.min(y1, y);
					y2 = Math.max(y2, y);
					count ++;
					}
				}
			}

		if (count <= 1)
			return mDefaultAVBL;

		double dx = x2 - x1;
		double dy = y2 - y1;
		double meanEdgeLength = (dx + dy) / 2;
		double targetArea = 3.0 * count;  // we assume 2/3 of empty space
		double area = dx * dy;
		double b = meanEdgeLength / (1 - targetArea);
		return Math.sqrt(b * b - area / (1 - targetArea)) - b;
		}

	private double calculateAverageBondLength() {
    	float averageBondLength = 0.0f;
        int bondCount = 0;
        if (mMolecule != null) {
            for (int i=0; i<mMolecule.length; i++) {
                if (mMolecule[i].getAllAtoms() != 0) {
                    if (mMolecule[i].getAllBonds() != 0) {
                        averageBondLength += (float)(mDepictor[i].getTransformation().getScaling()
			                                * mMolecule[i].getAllBonds() * mMolecule[i].getAverageBondLength());
                        bondCount += mMolecule[i].getAllBonds();
                        }
                    else {
						averageBondLength += (float)(mDepictor[i].getTransformation().getScaling()
											* mMolecule[i].getAverageBondLength());
						bondCount ++;
                        }
                    }
                }
            }
        return (bondCount == 0) ? mDefaultAVBL
                                : mTransformation.getScaling() * averageBondLength / bondCount;
        }

    private void layoutReaction(GenericDrawContext context) {
        GenericRectangle[] boundingRect = new GenericRectangle[mMolecule.length];
        double totalWidth = 0.0;
        double totalHeight = 0.0;
        for (int i=0; i<mMolecule.length; i++) {
            mDepictor[i].validateView(context, null, AbstractDepictor.cModeInflateToMaxAVBL);
            boundingRect[i] = mDepictor[i].getBoundingRect();
            totalWidth += boundingRect[i].width;
            totalHeight = Math.max(totalHeight, boundingRect[i].height);
            }

        final double catalystScale = 0.7;
		double catalystSpacing = 0.5 * AbstractDepictor.cOptAvBondLen;
	    GenericRectangle[] catalystBoundingRect = new GenericRectangle[mCatalyst.length];
		double totalCatalystWidth = 0.0;
		double totalCatalystHeight = 0.0;
		for (int i=0; i<mCatalyst.length; i++) {
			mCatalystDepictor[i].validateView(context, null, AbstractDepictor.cModeInflateToMaxAVBL+(int)(catalystScale*AbstractDepictor.cOptAvBondLen));
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
            if (i == mReactantCount) {
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
