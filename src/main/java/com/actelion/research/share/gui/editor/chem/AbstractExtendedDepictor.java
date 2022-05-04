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

package com.actelion.research.share.gui.editor.chem;

import com.actelion.research.chem.AbstractDepictor;
import com.actelion.research.chem.DepictorTransformation;
import com.actelion.research.chem.ExtendedMolecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.share.gui.Arrow;
import com.actelion.research.share.gui.DrawConfig;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;


public abstract class AbstractExtendedDepictor<T, C>
{
    protected StereoMolecule[] mMolecule;
    protected AbstractDepictor[] mDepictor;
    protected List<IDrawingObject> mDrawingObjectList;
    protected int mDisplayMode, mReactantOrCoreCount;
    protected boolean mUseGraphics2D, mDoLayoutMolecules, mIsMarkushStructure;
    protected DepictorTransformation mTransformation;
    private DrawConfig drawConfig;
    protected C mFragmentNoColor;

    public AbstractExtendedDepictor(StereoMolecule mol, java.util.List<IDrawingObject> drawingObjectList, boolean useGraphics2D, DrawConfig cfg)
    {
        this(new StereoMolecule[] {mol},drawingObjectList,useGraphics2D,cfg);
    }

    public AbstractExtendedDepictor(StereoMolecule[] mol, java.util.List<IDrawingObject> drawingObjectList, boolean useGraphics2D, DrawConfig cfg)
    {
        initialize(mol,drawingObjectList,useGraphics2D,cfg);
        mReactantOrCoreCount = -1;
    }

    /**
     * Use this constructor for markush structures. The first fragments in the list
     * are the Markush core structures (typically only one), decorated with R1,R2,R3,...
     * The remaining fragments need to contain one atom with atomicNo=0 each, that
     * indicates the attachment point. They also may contain Rn atoms.
     * Any of the fragments may contain query features.
     *
     * @param mol
     * @param markushCoreCount
     * @param drawingObjectList
     * @param useGraphics2D
     */
    public AbstractExtendedDepictor(StereoMolecule[] mol, int markushCoreCount, java.util.List<IDrawingObject> drawingObjectList, boolean useGraphics2D, DrawConfig cfg)
    {
        this(mol,drawingObjectList,useGraphics2D,cfg);
        mReactantOrCoreCount = markushCoreCount;
    }

    public AbstractExtendedDepictor(Reaction reaction, java.util.List<IDrawingObject> drawingObjectList, boolean layoutReaction, boolean useGraphics2D, DrawConfig cfg)
    {
//        System.out.printf("AED ctor %s\n",reaction.getReactants());
        StereoMolecule[] mol = null;
        if (reaction != null) {
            mol = new StereoMolecule[reaction.getMolecules()];
            for (int i = 0; i < reaction.getMolecules(); i++)
               mol[i] = reaction.getMolecule(i);
            mReactantOrCoreCount = reaction.getReactants();
            mDoLayoutMolecules = layoutReaction;
        }
        initialize(mol,drawingObjectList,useGraphics2D,cfg);
    }

    public abstract AbstractDepictor createDepictor(StereoMolecule stereoMolecule, boolean mUseGraphics2D,DrawConfig cfg);

    protected abstract void paintFragmentNumbers(T g);

    protected abstract void paintStructures(T g);

    protected abstract void paintDrawingObjects(T g);

    private void initialize(
        StereoMolecule[] mol,
        java.util.List<IDrawingObject> drawingObjectList,
        boolean useGraphics2D,
        DrawConfig cfg
    )
    {

        mMolecule = mol;
        mIsMarkushStructure = false;
        mDrawingObjectList = drawingObjectList;
        mUseGraphics2D = useGraphics2D;
        drawConfig = cfg;


        // for reactions and sets of molecules the availability of coordinates
        // is mandatory. However, every individual molecule may have its first
        // atom at coords 0.0/0.0, e.g. if they are encoded idcodes or if
        mTransformation = new DepictorTransformation();
        if (mMolecule != null) {
            mDepictor = new AbstractDepictor[mMolecule.length];
            for (int i = 0; i < mMolecule.length; i++) {
                mDepictor[i] = createDepictor(mMolecule[i], mUseGraphics2D, drawConfig);
/*
                if (mUseGraphics2D)
                mDepictor[i] = new Depictor2D(mMolecule[i]);
            else
                mDepictor[i] = new Depictor(mMolecule[i]);
*/
            }
        }
    }


    public void setDisplayMode(int displayMode)
    {
        mDisplayMode = displayMode;
    }

    public void setFragmentNoColor(C c)
    {
        // use setFragmentNoColor(null) if you don't want fragment numbers to be shown
        mFragmentNoColor = c;
    }

    public int getMoleculeCount()
    {
        return mMolecule.length;
    }

    public ExtendedMolecule getMolecule(int i)
    {
        return mMolecule[i];
    }

    public AbstractDepictor getMoleculeDepictor(int i)
    {
        return mDepictor[i];
    }

    public void setOverruleColor(Color foreGround, Color background)
    {
        if (mDepictor != null)
            for (AbstractDepictor d : mDepictor)
                d.setOverruleColor(foreGround, background);
    }

    public void setForegroundColor(Color foreGround, Color background) {
        if (mDepictor != null)
            for (AbstractDepictor d:mDepictor)
                d.setForegroundColor(foreGround, background);
    }


    public void paint(T g)
    {
        paintFragmentNumbers(g);
        paintStructures(g);
        paintDrawingObjects(g);
    }


    public DepictorTransformation updateCoords(T g, GenericRectangle viewRect, int mode)
    {
        // returns full transformation that moves/scales original molecules/objects into viewRect
        validateView(g, viewRect, mode);

//        System.out.printf("UpdateCoordinates %s\n",mTransformation);
        if (mTransformation.isVoidTransformation()) {
            return null;
        } else {
            if (mMolecule != null)
                for (int i = 0; i < mMolecule.length; i++)
                    mTransformation.applyTo(mMolecule[i]);

            if (mDrawingObjectList != null)
                for (int i = 0; i < mDrawingObjectList.size(); i++) {
                    IDrawingObject d = mDrawingObjectList.get(i);
                    GenericRectangle b  = d.getBoundingRect();
                    GenericRectangle t = new GenericRectangle(b.getX(), b.getY(), b.getWidth(), b.getHeight());
                    mTransformation.applyTo(t);
                }
            if (mDepictor != null)
                for (int i = 0; i < mDepictor.length; i++)
                    mDepictor[i].getTransformation().clear();

            DepictorTransformation t = mTransformation;
            mTransformation = new DepictorTransformation();
            return t;
        }
    }

    public DepictorTransformation validateView(T g, GenericRectangle viewRect, int mode)
    {
        // returns incremental transformation that moves/scales already transformed molecules/objects into viewRect
        if (mDoLayoutMolecules)
            doLayoutMolecules(g);

        GenericRectangle boundingRect = null;
        if (mDepictor != null) {
            for (int i = 0; i < mDepictor.length; i++) {
                mDepictor[i].validateView(g, null, 0);
                boundingRect = (boundingRect == null) ? mDepictor[i].getBoundingRect()
                        : boundingRect.union(mDepictor[i].getBoundingRect());
            }
        }
        if (mDrawingObjectList != null) {
            for (int i = 0; i < mDrawingObjectList.size(); i++) {
                GenericRectangle b = mDrawingObjectList.get(i).getBoundingRect();
                GenericRectangle objectBounds = new GenericRectangle(b.getX(), b.getY(), b.getWidth(), b.getHeight());
                mTransformation.applyTo(objectBounds);
                boundingRect = (boundingRect == null) ? objectBounds
                        : boundingRect.union(objectBounds);
            }
        }

        if (boundingRect == null)
            return null;

        double avbl = calculateAverageBondLength();

        DepictorTransformation t = new DepictorTransformation(boundingRect, viewRect, avbl, mode);

        if (!t.isVoidTransformation()) {
            t.applyTo(mTransformation);

            if (mDepictor != null)
                for (int i = 0; i < mDepictor.length; i++)
                    mDepictor[i].applyTransformation(t);

            return t;
        }

        return null;
    }

    protected double calculateAverageBondLength()
    {
        float averageBondLength = 0.0f;
        int bondCount = 0;
        if (mMolecule != null) {
            for (int i = 0; i < mMolecule.length; i++) {
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

    private void doLayoutMolecules(T g)
    {
        GenericRectangle[] boundingRect = new GenericRectangle[mMolecule.length];
        double totalWidth = 0.0f;
        double totalHeight = 0.0f;
        for (int i = 0; i < mMolecule.length; i++) {
            mDepictor[i].validateView(g, null, AbstractDepictor.cModeInflateToMaxAVBL);
            boundingRect[i] = mDepictor[i].getBoundingRect();
            totalWidth += boundingRect[i].width;
            totalHeight = Math.max(totalHeight, boundingRect[i].height);
        }

        float spacing = 1.5f * AbstractDepictor.cOptAvBondLen;
        float arrowWidth = 2f * AbstractDepictor.cOptAvBondLen;

        int arrow = -1;
        if (mDrawingObjectList == null) {
            mDrawingObjectList = new ArrayList<IDrawingObject>();
            mDrawingObjectList.add(new Arrow(drawConfig, 0,0,0,0));
            arrow = 0;
        } else {
//            for (int i = 0; i < mDrawingObjectList.size(); i++) {
//                if (mDrawingObjectList.get(i) instanceof ReactionArrow) {
//                    arrow = i;
//                    break;
//                }
//            }
            if (mDrawingObjectList.size()== 0)
            //if (arrow == -1)
            {
                arrow = mDrawingObjectList.size();
                mDrawingObjectList.add(new Arrow(drawConfig,0,0,0,0));
            }
        }

        float rawX = 0.5f * spacing;
        for (int i = 0; i < mMolecule.length; i++) {
            if (i == mReactantOrCoreCount) {
//                ((ReactionArrow) mDrawingObjectList.get(arrow)).setCoordinates(
//                        rawX - spacing / 2, totalHeight / 2, rawX - spacing / 2 + arrowWidth, totalHeight / 2);
                rawX += arrowWidth;
            }

            double dx = rawX - boundingRect[i].x;
            double dy = 0.5 * (totalHeight - boundingRect[i].height) - boundingRect[i].y;
            mDepictor[i].applyTransformation(new DepictorTransformation(1.0, dx, dy));

            rawX += spacing + boundingRect[i].width;
        }

        mDoLayoutMolecules = false;
    }

    public DepictorTransformation getmTransformation()
    {
        return mTransformation;
    }
}
