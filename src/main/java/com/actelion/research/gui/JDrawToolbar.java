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

package com.actelion.research.gui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

public class JDrawToolbar extends JComponent
        implements MouseListener,MouseMotionListener {
    static final long serialVersionUID = 0x20090402;

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
//	protected static final int cButtonChiral = 20;	// in multi fragment mode this button is used as tool
//	protected static final int cToolChiral = 20;
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

    private static final int cESRMenuBorder = 4;
    private static final int cESRMenuX = 10;
    private static final int cESRMenuY = 54;
    protected static final int cToolESRAbs = 101;
    protected static final int cToolESROr = 102;
    protected static final int cToolESRAnd  = 103;

	private JDrawArea	mArea;
	private Image		mImageUp,mImageDown,mESRImageUp,mESRImageDown;
	private int			mCurrentTool,mPressedButton,mMode,mESRSelected,mESRHilited;
    private boolean     mESRMenuVisible;

	public JDrawToolbar(JDrawArea theArea) {
		mArea = theArea;
		init();
		}


	public JDrawToolbar(JDrawArea theArea, int mode) {
		mArea = theArea;
		mMode = mode;
		if ((mMode & JDrawArea.MODE_REACTION) != 0)
			mMode |= JDrawArea.MODE_MULTIPLE_FRAGMENTS;
		init();
		}

    // CXR added this 09/05/2011
    public void setReactionMode(boolean rxn) {
        if (rxn)  {
            mMode = JDrawArea.MODE_MULTIPLE_FRAGMENTS | JDrawArea.MODE_REACTION;
        } else
            mMode &= ~JDrawArea.MODE_REACTION;
        }


	private void init() {
		MediaTracker t = new MediaTracker(this);

		mImageDown = getImageResource("drawButtonsDown.gif");
		mImageUp = getImageResource("drawButtonsUp.gif");
        mESRImageDown = getImageResource("ESRButtonsDown.gif");
        mESRImageUp = getImageResource("ESRButtonsUp.gif");

		t.addImage(mImageUp,0);
		t.addImage(mImageDown,0);
        t.addImage(mESRImageDown,0);
        t.addImage(mESRImageDown,0);
		try {
			t.waitForAll();
			}
		catch (InterruptedException e) {}

		int w = mImageUp.getWidth(this);
		int h = mImageUp.getHeight(this);
		setMinimumSize(new Dimension(w,h));
		setPreferredSize(new Dimension(w,h));
		setSize(w,h);

		mCurrentTool = cToolStdBond;
		mPressedButton = -1;

		addMouseListener(this);
        addMouseMotionListener(this);
		}

	public void setCurrentTool(int tool) {
		if (mCurrentTool != tool) {
			mCurrentTool = tool;
			mArea.toolChanged(tool);
			repaint();
			}
		}

	private Image getImageResource(String imageFileName) {
		Image im=null;
		try {
            ImageIcon ic = new ImageIcon(getClass().getResource("/images/"+imageFileName));
            im = ic.getImage();
			}
		catch (Exception e) {
			System.out.println("Error loading images: "+e);
			}
		return im;
		}

	public void paintComponent(Graphics g) {
        super.paintComponent(g);

        g.drawImage(mImageUp,0,0,this);
        drawPressedButton(g, mCurrentTool);
        if (mPressedButton != -1)
            drawPressedButton(g, mPressedButton);

        if (mESRSelected != 0) {
            setButtonClip(g, cToolESR);
            Point l = getButtonLocation(cToolESR);
            Image esrButtons = (mCurrentTool == cToolESR) ? mESRImageDown : mESRImageUp;
            g.drawImage(esrButtons, l.x-cESRMenuBorder, l.y-cESRMenuBorder-21*mESRSelected, this);
            g.setClip(0, 0, Integer.MAX_VALUE, Integer.MAX_VALUE);
            }

        if (mESRMenuVisible) {
            g.drawImage(mESRImageUp, cESRMenuX-cESRMenuBorder, cESRMenuY-cESRMenuBorder, this);
            if (mESRHilited != -1) {
                g.setClip(cESRMenuX, cESRMenuY+21*mESRHilited, 20, 20);
                g.drawImage(mESRImageDown, cESRMenuX-cESRMenuBorder, cESRMenuY-cESRMenuBorder, this);
                g.setClip(0, 0, Integer.MAX_VALUE, Integer.MAX_VALUE);
               }
            }
        }

	public void update(Graphics g) {
		paint(g);
		}

	public void mousePressed(MouseEvent e) {
		int b = getButtonNo(e);
		if (b == -1) return;

        if (b == cToolESR) {
            mESRMenuVisible = true;
            validateESRHiliting(e);
            repaint();
            }

        mPressedButton = b;

        if (b != mCurrentTool)
            repaint();
		}

	public void mouseReleased(MouseEvent e) {
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
		 || (mPressedButton == cToolMapper && (mMode & JDrawArea.MODE_REACTION) == 0)
		 || (mPressedButton == cToolText && (mMode & JDrawArea.MODE_DRAWING_OBJECTS) == 0)) {
			mPressedButton = -1;
            repaint();
			return;
			}

        mPressedButton = -1;
		if (releasedButton == cButtonClear
		 || releasedButton == cButtonCleanStructure
		 || releasedButton == cButtonUndo) {
//		 || (releasedButton == cButtonChiral && (mMode & JDrawArea.MODE_MULTIPLE_FRAGMENTS) == 0)) {
            repaint();
			mArea.buttonPressed(releasedButton);
			return;
			}
		mCurrentTool = releasedButton;
        repaint();

        if (mCurrentTool == cToolESR)
            mArea.toolChanged((mESRSelected == 0) ? cToolESRAbs
                            : (mESRSelected == 1) ? cToolESROr : cToolESRAnd);
        else
            mArea.toolChanged(releasedButton);
		}

	public void mouseEntered(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}
	public void mouseClicked(MouseEvent e) {}

    public void mouseMoved(MouseEvent e) {}
    public void mouseDragged(MouseEvent e) {
        int oldESRHilited = mESRHilited;
        validateESRHiliting(e);
        if (oldESRHilited != mESRHilited)
            repaint();
        }

 	private int getButtonNo(MouseEvent e) {
        int x = e.getX();
        int y = e.getY();
		if (x<0 || x>44 || y<0 || y>cButtonsPerColumn*21) return -1;
		x -= 3;
		y -= 3;
		if ((x % 21) > 18) return -1;
		if ((y % 21) > 18) return -1;
		return cButtonsPerColumn * (x/21) + (y/21);
		}

    private void validateESRHiliting(MouseEvent e) {
        int x = e.getX();
        int y = e.getY();
        if (x>cESRMenuX && x<cESRMenuX+21) {
            mESRHilited = (y-cESRMenuY) / 21;
            if (mESRHilited < 0 || mESRHilited > 2)
                mESRHilited = -1;
            }
        }

    private void setButtonClip(Graphics g, int button) {
        Point l = getButtonLocation(button);
        g.setClip(l.x, l.y,20,20);
        }

	private void drawPressedButton(Graphics g, int button) {
        setButtonClip(g, button);
		g.drawImage(mImageDown,0,0,this);
        g.setClip(0,0,Integer.MAX_VALUE,Integer.MAX_VALUE);
		}

    private Point getButtonLocation(int button) {
        return new Point(21 * (button / cButtonsPerColumn) + 2,
                         21 * (button % cButtonsPerColumn) + 2);
        }
	}
