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
package com.actelion.research.gui.wmf;

import java.awt.*;


public abstract class MetaFile
{

	public abstract void arc(int left, int top, int right, int bottom, int xstart, int ystart, int xend, int yend);
	public abstract int createBrushIndirect(int style, Color color, int hatch);
	public abstract int createFont(Font font, int i, boolean flag, boolean flag1);
	public abstract int createFont(
		int height, int with, int esc, int orient, int weight,
		boolean italic, boolean underline, boolean strikeOut, byte charSet,
		byte outPrecision, byte clipPrecision, byte quality, byte pitchAndFamily, String s);
	public abstract int createPatternBrush(int[] ai, int i, int j);
	public abstract int createPenIndirect(int style, int width, Color color);
	public abstract void deleteObject(int i);
	public abstract void ellipse(int left, int top, int right, int bottom);
	public abstract void escape(int function, byte[] data);
	public abstract void intersectClipRect(int i, int j, int k, int l);
	public abstract void lineTo(int x, int y);
	public abstract void moveTo(int x, int y);
	public abstract void pie(int left, int top, int right, int bottom, int xR1, int yR1, int xR2, int yR2);
	public abstract void polygon(int[] ptx, int[] pty, int count);
	public abstract void polyline(int[] ptx, int[] pty, int count);
	public abstract void polypolygon(Polygon[] apolygon);
	public abstract void rectangle(int left, int top, int right, int bottom);
	public abstract void roundRect(int left, int top, int right, int bottom, int width, int height);
	public abstract void selectObject(int handle);
	public abstract void setBKColor(Color color);
	public abstract void setBKMode(int mode);
	public abstract void setClipRgn();
	public abstract void setMapMode(int mode);
	public abstract void setPixel(int x, int y, Color color);
	public abstract void setPolyFillMode(int mode);
	public abstract void setROP2(int mode);
	public abstract void setStretchBltMode(int mode);
	public abstract void setTextAlign(int i);
	public abstract void setTextCharacterExtra(int i);
	public abstract void setTextColor(Color color);
	public abstract void setViewportExt(int i, int j);
	public abstract void setWindowExt(int cx, int cy);
	public abstract void setWindowOrg(int x, int y);
	public abstract void stretchBlt(int xOrigDest, int yOrigDest, int widthDest, int heightDest,
									int xOrigSrc, int yOrigSrc, int widthSrc, int heightSrc, int rasterOp,
									int[] pixelData, int imageWidth, int imageHeight);
	public abstract void textOut(int x, int y, String s);
	public abstract String translateFontName(String s);
}
