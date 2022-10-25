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

import com.actelion.research.gui.LookAndFeelHelper;
import com.actelion.research.gui.generic.*;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.actelion.research.util.ColorHelper;

public class GenericEditorToolbar implements GenericEventListener<GenericMouseEvent> {
    protected static final int cButtonsPerColumn = 17;

	private static final int cImageOversize = 4;  // source image size in regard to target size with no UI-scaling

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
	protected static final int cToolCustomAtom = 33;

	protected static final int cBorder = Math.round(HiDPIHelper.scale(2f));
	protected static final float cSourceButtonSize = cImageOversize*21f;
	protected static final float cButtonSize = HiDPIHelper.getUIScaleFactor()*21f;
    protected static final int cToolESRAbs = 101;
    protected static final int cToolESROr = 102;
    protected static final int cToolESRAnd  = 103;

	private GenericCanvas mToolbarCanvas;
	private GenericEditorArea mArea;
	private GenericImage mImageNormal, mImageDisabled, mESRImageNormal;
	protected int mWidth,mHeight,mCurrentTool, mSelectedButton,mHighlightedButton,mESRSelected;
	private float mImageScaling;

	public GenericEditorToolbar(GenericCanvas toolbarCanvas, GenericEditorArea theArea) {
		mToolbarCanvas = toolbarCanvas;
		mArea = theArea;
		init();
		}

	public int getWidth() {
		return mWidth;
		}

	public int getHeight() {
		return mHeight;
		}

	private void init() {
		mImageNormal = mArea.getUIHelper().createImage("editorButtons.png");
		if (LookAndFeelHelper.isDarkLookAndFeel())
			HiDPIHelper.adaptForLookAndFeel(mImageNormal);
		mImageDisabled = mArea.getUIHelper().createImage("editorButtons.png");
		HiDPIHelper.disableImage(mImageDisabled);

		mWidth =  2*cBorder + HiDPIHelper.scale(mImageNormal.getWidth() / cImageOversize);
		mHeight = 2*cBorder + HiDPIHelper.scale(mImageNormal.getHeight() / cImageOversize);

		mESRImageNormal = mArea.getUIHelper().createImage("esrButtons.png");

		mImageScaling = cImageOversize / HiDPIHelper.getUIScaleFactor();

		mCurrentTool = cToolStdBond;
		mSelectedButton = -1;
		mHighlightedButton = -1;
	}

	public void setCurrentTool(int tool) {
		if (mCurrentTool != tool) {
			mCurrentTool = tool;
			mArea.toolChanged(tool);
			mToolbarCanvas.repaint();
			}
		}

	public void paintContent(GenericDrawContext context) {
		int background = mToolbarCanvas.getBackgroundRGB();
		boolean isDark = (ColorHelper.perceivedBrightness(background) < 0.5);
		int highlightBackground = isDark ? ColorHelper.brighter(background, 0.6f) : ColorHelper.darker(background, 0.6f);
		int selectedBackground = isDark ? ColorHelper.brighter(background, 0.8f) : ColorHelper.darker(background, 0.8f);

		int sw = mImageNormal.getWidth();
		int sh = mImageNormal.getHeight();
		context.drawImage(mImageNormal, 0, 0, sw, sh, cBorder, cBorder, sw/mImageScaling, sh/mImageScaling);

		// draw the selected ESR button
		double[] l = getButtonLocation(cToolESR);
		context.drawImage(mESRImageNormal, 0, mESRSelected*cSourceButtonSize,
				cSourceButtonSize, cSourceButtonSize, l[0], l[1], cButtonSize, cButtonSize);

		// draw disabled buttons just over originals
		if ((mArea.getMode() & GenericEditorArea.MODE_REACTION) == 0)
			drawButton(context, cToolMapper, -1, true);
		if ((mArea.getMode() & GenericEditorArea.MODE_DRAWING_OBJECTS) == 0)
			drawButton(context, cToolText, -1, true);

        drawButton(context, mCurrentTool, selectedBackground, false);

		if (mHighlightedButton != -1 && mHighlightedButton != mSelectedButton)
			drawButton(context, mHighlightedButton, highlightBackground, false);

        if (mSelectedButton != -1)
	        drawButton(context, mSelectedButton, 0x006D5FB4, false);
        }

	@Override
	public void eventHappened(GenericMouseEvent e) {
		if (e.getWhat() == GenericMouseEvent.MOUSE_PRESSED) {
			int b = getButtonNo(e);
			if (!isSelectableButton(b))
				return;

	        if (b == cToolESR && b == mCurrentTool)
		        mESRSelected = (++mESRSelected) % 3;

	        mSelectedButton = b;

	        mToolbarCanvas.repaint();
			}
		else if (e.getWhat() == GenericMouseEvent.MOUSE_RELEASED) {
	        if (mSelectedButton == -1)
	            return;

            int releasedButton = getButtonNo(e);

	        if (releasedButton != mSelectedButton
			 || (mSelectedButton == cToolMapper && (mArea.getMode() & GenericEditorArea.MODE_REACTION) == 0)
			 || (mSelectedButton == cToolText && (mArea.getMode() & GenericEditorArea.MODE_DRAWING_OBJECTS) == 0)) {
				mSelectedButton = -1;
		        mToolbarCanvas.repaint();
				return;
				}

	        mSelectedButton = -1;
			if (releasedButton == cButtonClear
			 || releasedButton == cButtonCleanStructure
			 || releasedButton == cButtonUndo) {
				mToolbarCanvas.repaint();
				mArea.buttonPressed(releasedButton);
				return;
				}
			mCurrentTool = releasedButton;
			mToolbarCanvas.repaint();

	        if (mCurrentTool == cToolESR) {
		        mArea.toolChanged((mESRSelected == 0) ? cToolESRAbs
						        : (mESRSelected == 1) ? cToolESRAnd : cToolESROr);
	            }
	        else if (mCurrentTool == cToolCustomAtom) {
	        	mArea.showCustomAtomDialog(-1);
		        mArea.toolChanged(releasedButton);
	            }
	        else {
		        mArea.toolChanged(releasedButton);
	            }
			}
		else if (e.getWhat() == GenericMouseEvent.MOUSE_MOVED
			  || e.getWhat() == GenericMouseEvent.MOUSE_EXITED) {
			int b = getButtonNo(e);
			if (b == mSelectedButton)   // cannot highlight selected button
				b = -1;
			if (b != mHighlightedButton) {
				mHighlightedButton = b;
				mToolbarCanvas.repaint();
				}
			}
		}

 	protected int getButtonNo(GenericMouseEvent e) {
        int x = e.getX() - cBorder;
        int y = e.getY() - cBorder;
		if (x<0 || x>=2*cButtonSize || y<0 || y>=cButtonsPerColumn*cButtonSize)
			return -1;
		int button = cButtonsPerColumn * (int)(x/cButtonSize) + (int)(y/cButtonSize);
		return isSelectableButton(button) ? button : -1;
		}

	private boolean isSelectableButton(int b) {
		return b >= 0 && b < 2*cButtonsPerColumn
			&& (b != mCurrentTool || b == cToolESR || b == cToolCustomAtom)
			&& (b != cToolMapper || (mArea.getMode() & GenericEditorArea.MODE_REACTION) != 0)
			&& (b != cToolText || (mArea.getMode() & GenericEditorArea.MODE_DRAWING_OBJECTS) != 0);
		}

	private void drawButton(GenericDrawContext context, int button, int background, boolean isDisabled) {
		background |= 0x80000000;
        double[] bl = getButtonLocation(button);
        if (background != -1) {
	        context.setRGB(background);
	        context.fillRectangle(bl[0], bl[1], cButtonSize, cButtonSize);
            }
        if (button == cToolESR)
	        context.drawImage(mESRImageNormal, 0, mESRSelected*cSourceButtonSize, cSourceButtonSize, cSourceButtonSize,
			        bl[0], bl[1], cButtonSize, cButtonSize);
        else
			context.drawImage(isDisabled ? mImageDisabled : mImageNormal, (bl[0] - cBorder)*mImageScaling, (bl[1] - cBorder)*mImageScaling,
					cSourceButtonSize, cSourceButtonSize, bl[0], bl[1], cButtonSize, cButtonSize);
		}

	private double[] getButtonLocation(int button) {
	    double[] p = new double[2];
		p[0] = cButtonSize * (button / cButtonsPerColumn) + cBorder;
		p[1] = cButtonSize * (button % cButtonsPerColumn) + cBorder;
        return p;
		}
	}
