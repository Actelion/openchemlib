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

package com.actelion.research.gui.editor;

import com.actelion.research.gui.generic.*;
import com.actelion.research.gui.hidpi.HiDPIHelper;

public class GenericEditorToolbar implements GenericEventListener<GenericMouseEvent> {
    protected static final int cButtonsPerColumn = 17;

	protected static final int cButtonClear = 0;
	protected static final int cButtonCleanStructure = 1;
	public static final int cToolLassoPointer = 2;
	protected static final int cToolUnknownParity = 3;
	protected static final int cToolDelete = 4;
	protected static final int cToolStdBond = 5;
	protected static final int cToolUpBond = 6;
	protected static final int cTool3Ring = 7;
	protected static final int cTool5Ring = 8;
	protected static final int cTool7Ring = 9;
	protected static final int cToolPosCharge = 10;
	protected static final int cToolAtomC = 11;
	protected static final int cToolAtomN = 12;
	protected static final int cToolAtomO = 13;
	protected static final int cToolAtomF = 14;
	protected static final int cToolAtomBr = 15;
	protected static final int cToolAtomH = 16;
	protected static final int cButtonUndo = 17;
	protected static final int cToolZoom = 18;
	protected static final int cToolMapper = 19;
    protected static final int cToolESR = 20;
	protected static final int cToolText = 21;
	protected static final int cToolChain = 22;
	protected static final int cToolDownBond = 23;
	protected static final int cTool4Ring = 24;
	protected static final int cTool6Ring = 25;
	protected static final int cToolAromRing = 26;
	protected static final int cToolNegCharge = 27;
	protected static final int cToolAtomSi = 28;
	protected static final int cToolAtomP = 29;
	protected static final int cToolAtomS = 30;
	protected static final int cToolAtomCl = 31;
	protected static final int cToolAtomI = 32;
	protected static final int cToolAtomOther = 33;

    protected static final int cESRMenuBorder = HiDPIHelper.scale(4);
	protected static final int cESRMenuOffsetX = HiDPIHelper.scale(-5);
	protected static final int cESRMenuOffsetY = HiDPIHelper.scale(-5);
	protected static final float cButtonBorder = HiDPIHelper.getUIScaleFactor()*3f;
	protected static final float cButtonSize = HiDPIHelper.getUIScaleFactor()*21f;
    protected static final int cToolESRAbs = 101;
    protected static final int cToolESROr = 102;
    protected static final int cToolESRAnd  = 103;

	private GenericCanvas mToolbarCanvas;
	private GenericDrawArea mArea;
	private GenericImage mImageUp,mImageDown,mESRImageUp,mESRImageDown;
	protected int		mWidth,mHeight,mCurrentTool,mPressedButton,mMode,mESRSelected,mESRHilited;
    protected boolean   mESRMenuVisible;

	public GenericEditorToolbar(GenericCanvas toolbarCanvas, GenericDrawArea theArea) {
		mToolbarCanvas = toolbarCanvas;
		mArea = theArea;
		init();
		}

	public GenericEditorToolbar(GenericCanvas toolbarCanvas, GenericDrawArea theArea, int mode) {
		mToolbarCanvas = toolbarCanvas;
		mArea = theArea;
		mMode = mode;
		if ((mMode & GenericDrawArea.MODE_REACTION) != 0)
			mMode |= GenericDrawArea.MODE_MULTIPLE_FRAGMENTS;
		init();
		}

	public int getWidth() {
		return mWidth;
		}

	public int getHeight() {
		return mHeight;
		}

	/* CXR added this 09/05/2011
    public void setReactionMode(boolean rxn) {
        if (rxn)  {
            mMode = GenericDrawArea.MODE_MULTIPLE_FRAGMENTS | GenericDrawArea.MODE_REACTION;
        } else
            mMode &= ~GenericDrawArea.MODE_REACTION;
        }*/

	private void init() {
		mImageDown = mArea.getUIHelper().createImage("drawButtonsDown.gif");
		mImageUp = mArea.getUIHelper().createImage("drawButtonsUp.gif");

		mWidth = mImageUp.getWidth();
		mHeight = mImageUp.getHeight();

		mESRImageDown = mArea.getUIHelper().createImage("ESRButtonsDown.gif");
		mESRImageUp = mArea.getUIHelper().createImage("ESRButtonsUp.gif");

		scaleImages();

		mCurrentTool = cToolStdBond;
		mPressedButton = -1;
		}

	private void scaleImages() {
		if (cButtonSize == 21)	// no scaling
			return;

		mWidth = HiDPIHelper.scale(mWidth);
		mHeight = HiDPIHelper.scale(mHeight);

		mImageDown.scale(mWidth, mHeight);
		mImageUp.scale(mWidth, mHeight);

		int width = HiDPIHelper.scale(mESRImageUp.getWidth());
		int height = HiDPIHelper.scale(mESRImageUp.getHeight());

		mESRImageDown.scale(width, height);
		mESRImageUp.scale(width, height);
		}

	public void setCurrentTool(int tool) {
		if (mCurrentTool != tool) {
			mCurrentTool = tool;
			mArea.toolChanged(tool);
			mToolbarCanvas.repaint();
			}
		}

	public void paintContent(GenericDrawContext context) {
		context.drawImage(mImageUp, 0,0);
        drawPressedButton(context, mCurrentTool);
        if (mPressedButton != -1)
            drawPressedButton(context, mPressedButton);

        if (mESRSelected != 0) {
	        double[] l = getButtonLocation(cToolESR);
            GenericImage esrButtons = (mCurrentTool == cToolESR) ? mESRImageDown : mESRImageUp;
	        context.drawImage(esrButtons, cESRMenuBorder, cESRMenuBorder+mESRSelected*cButtonSize, l[0], l[1], cButtonSize, cButtonSize);
            }

		if (mESRMenuVisible) {
	        double[] l = getButtonLocation(cToolESR);
	        double esrMenuX = l[0]-cESRMenuBorder+cESRMenuOffsetX;
	        double esrMenuY = l[1]-cESRMenuBorder+cESRMenuOffsetY-cButtonSize*mESRSelected;
	        context.drawImage(mESRImageUp, esrMenuX, esrMenuY);
            if (mESRHilited != -1) {
	            context.drawImage(mESRImageDown, cESRMenuBorder, cESRMenuBorder+mESRHilited*cButtonSize,
			            l[0]+cESRMenuOffsetX, l[1]+cESRMenuOffsetY+(mESRHilited-mESRSelected)*cButtonSize, cButtonSize, cButtonSize);
                }
            }
        }

	@Override
	public void eventHappened(GenericMouseEvent e) {
		if (e.getWhat() == GenericMouseEvent.MOUSE_PRESSED) {
			int b = getButtonNo(e);
			if (b == -1) return;

	        if (b == cToolESR) {
	            mESRMenuVisible = true;
	            validateESRHiliting(e);
		        mToolbarCanvas.repaint();
	            }

	        mPressedButton = b;

	        if (b != mCurrentTool)
		        mToolbarCanvas.repaint();
			}
		else if (e.getWhat() == GenericMouseEvent.MOUSE_RELEASED) {
	        int releasedButton = -1;
	        if (mESRMenuVisible) {
	            if (mESRHilited != -1) {
	                mESRSelected = mESRHilited;
	                releasedButton = cToolESR;
	                }
	            mESRMenuVisible = false;
	            }

	        if (mPressedButton == -1)
	            return;

	            // if user didn't change esr menu than require that the button beneath
	            // mouse pointer was the same when pressing and releasing the mouse
	        if (releasedButton == -1)
	            releasedButton = getButtonNo(e);

	        if (releasedButton != mPressedButton
			 || (mPressedButton == cToolMapper && (mMode & GenericDrawArea.MODE_REACTION) == 0)
			 || (mPressedButton == cToolText && (mMode & GenericDrawArea.MODE_DRAWING_OBJECTS) == 0)) {
				mPressedButton = -1;
		        mToolbarCanvas.repaint();
				return;
				}

	        mPressedButton = -1;
			if (releasedButton == cButtonClear
			 || releasedButton == cButtonCleanStructure
			 || releasedButton == cButtonUndo) {
	//		 || (releasedButton == cButtonChiral && (mMode & JDrawArea.MODE_MULTIPLE_FRAGMENTS) == 0)) {
				mToolbarCanvas.repaint();
				mArea.buttonPressed(releasedButton);
				return;
				}
			mCurrentTool = releasedButton;
			mToolbarCanvas.repaint();

	        if (mCurrentTool == cToolESR)
	            mArea.toolChanged((mESRSelected == 0) ? cToolESRAbs
	                            : (mESRSelected == 1) ? cToolESROr : cToolESRAnd);
	        else
	            mArea.toolChanged(releasedButton);
			}
 		else if (e.getWhat() == GenericMouseEvent.MOUSE_DRAGGED) {
	        int oldESRHilited = mESRHilited;
	        validateESRHiliting(e);
	        if (oldESRHilited != mESRHilited)
		        mToolbarCanvas.repaint();
	       }
		}

 	protected int getButtonNo(GenericMouseEvent e) {
        int x = e.getX();
        int y = e.getY();
		if (x<0 || x>=2*cButtonSize+cButtonBorder || y<0 || y>cButtonsPerColumn*cButtonSize) return -1;
		x -= cButtonBorder;
		y -= cButtonBorder;
		if ((x % cButtonSize) > cButtonSize-cButtonBorder) return -1;
		if ((y % cButtonSize) > cButtonSize-cButtonBorder) return -1;
		return cButtonsPerColumn * (int)(x/cButtonSize) + (int)(y/cButtonSize);
		}

    private void validateESRHiliting(GenericMouseEvent e) {
        int x = e.getX();
        int y = e.getY();
		mESRHilited = -1;
	    double[] l = getButtonLocation(cToolESR);
	    l[0] += cESRMenuOffsetX;
	    l[1] += cESRMenuOffsetX;
	    if (x>l[0] && x<l[0]+cButtonSize) {
            int b = (int)((y-l[1]+cButtonSize*mESRSelected) / cButtonSize);
            if (b >= 0 && b <= 2)
                mESRHilited = b;
            }
        }

	protected void drawPressedButton(GenericDrawContext context, int button) {
        double[] bl = getButtonLocation(button);
		context.drawImage(mImageDown, bl[0], bl[1], bl[0], bl[1], cButtonSize, cButtonSize);
		}

    private double[] getButtonLocation(int button) {
	    double[] p = new double[2];
		p[0] = cButtonSize * (button / cButtonsPerColumn) + cButtonBorder-2;
		p[1] = cButtonSize * (button % cButtonsPerColumn) + cButtonBorder-2;
        return p;
		}
	}
