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

import com.actelion.research.gui.generic.GenericCanvas;
import com.actelion.research.gui.generic.GenericDrawContext;
import com.actelion.research.gui.generic.GenericMouseEvent;
import com.actelion.research.gui.generic.GenericMouseListener;
import com.actelion.research.gui.hidpi.HiDPIHelper;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.net.URL;

public class GenericEditorToolbar implements GenericMouseListener {
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

    protected static final int cESRMenuBorder = 4;
	protected static final int cESRMenuX = HiDPIHelper.scale(20);
	protected static final int cESRMenuY = HiDPIHelper.scale(64);
	protected static final float cButtonBorder = HiDPIHelper.getUIScaleFactor()*3f;
	protected static final float cButtonSize = HiDPIHelper.getUIScaleFactor()*21f;
    protected static final int cToolESRAbs = 101;
    protected static final int cToolESROr = 102;
    protected static final int cToolESRAnd  = 103;

	private GenericCanvas mToolbarCanvas;
	private GenericDrawArea mArea;
	private Image		mImageUp,mImageDown,mESRImageUp,mESRImageDown;
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
		mImageDown = createImage("drawButtonsDown.gif");
		mImageUp = createImage("drawButtonsUp.gif");

		mWidth = ((BufferedImage)mImageUp).getWidth();
		mHeight = ((BufferedImage)mImageUp).getHeight();

		mESRImageDown = createImage("ESRButtonsDown.gif");
		mESRImageUp = createImage("ESRButtonsUp.gif");

		scaleImages();

		mCurrentTool = cToolStdBond;
		mPressedButton = -1;
		}

	private void scaleImages() {
		if (cButtonSize == 21)	// no scaling
			return;

		mWidth = HiDPIHelper.scale(mWidth);
		mHeight = HiDPIHelper.scale(mHeight);

		mImageDown = mImageDown.getScaledInstance(mWidth, mHeight, Image.SCALE_SMOOTH);
		mImageUp = mImageUp.getScaledInstance(mWidth, mHeight, Image.SCALE_SMOOTH);

		int width = HiDPIHelper.scale(((BufferedImage)mESRImageUp).getWidth());
		int height = HiDPIHelper.scale(((BufferedImage)mESRImageUp).getHeight());

		mESRImageDown = mESRImageDown.getScaledInstance(width, height, Image.SCALE_SMOOTH);
		mESRImageUp = mESRImageUp.getScaledInstance(width, height, Image.SCALE_SMOOTH);
		}

	public BufferedImage createImage(String fileName) {
		// Once double resolution is available use HiDPIHelper.createImage() !!!

		URL url = GenericEditorToolbar.class.getResource("/images/" + fileName);
		if (url == null)
			throw new RuntimeException("Could not find: " + fileName);

		try {
			return ImageIO.read(url);
			}
		catch (IOException ioe) { return null; }
		}

	public void setCurrentTool(int tool) {
		if (mCurrentTool != tool) {
			mCurrentTool = tool;
			mArea.toolChanged(tool);
			mToolbarCanvas.repaint();
			}
		}

	public void paintContent(GenericDrawContext context) {
		context.drawImage(0,0, mImageUp);
        drawPressedButton(context, mCurrentTool);
        if (mPressedButton != -1)
            drawPressedButton(context, mPressedButton);

        if (mESRSelected != 0) {
            setButtonClip(context, cToolESR);
            Point l = getButtonLocation(cToolESR);
            Image esrButtons = (mCurrentTool == cToolESR) ? mESRImageDown : mESRImageUp;
	        context.drawImage(l.x-cESRMenuBorder, l.y-cESRMenuBorder-cButtonSize*mESRSelected, esrButtons);
	        context.setClip(0, 0, Integer.MAX_VALUE, Integer.MAX_VALUE);
            }

        if (mESRMenuVisible) {
	        context.drawImage(cESRMenuX-cESRMenuBorder, cESRMenuY-cESRMenuBorder-cButtonSize*mESRSelected, mESRImageUp);
            if (mESRHilited != -1) {
	            context.setClip(cESRMenuX, Math.round(cESRMenuY-cButtonSize*mESRSelected+cButtonSize*mESRHilited), Math.round(cButtonSize-1), Math.round(cButtonSize-1));
	            context.drawImage(cESRMenuX-cESRMenuBorder, cESRMenuY-cESRMenuBorder-cButtonSize*mESRSelected, mESRImageDown);
	            context.setClip(0, 0, Integer.MAX_VALUE, Integer.MAX_VALUE);
               }
            }
        }

	@Override
	public void mouseActionHappened(GenericMouseEvent e) {
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
        if (x>=cESRMenuX && x<cESRMenuX+cButtonSize) {
            int b = (int)((y-cESRMenuY+cButtonSize*mESRSelected) / cButtonSize);
            if (b >= 0 && b <= 2)
                mESRHilited = b;
            }
        }

    protected void setButtonClip(GenericDrawContext context, int button) {
        Point l = getButtonLocation(button);
	    context.setClip(l.x, l.y,Math.round(cButtonSize-1),Math.round(cButtonSize-1));
        }

	protected void drawPressedButton(GenericDrawContext context, int button) {
        setButtonClip(context, button);
		context.drawImage(0,0,mImageDown);
		context.setClip(0,0,Integer.MAX_VALUE,Integer.MAX_VALUE);
		}

    protected Point getButtonLocation(int button) {
        return new Point(Math.round(cButtonSize * (button / cButtonsPerColumn) + cButtonBorder-2),
						 Math.round(cButtonSize * (button % cButtonsPerColumn) + cButtonBorder-2));
        }
	}
