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
import java.awt.image.ImageObserver;
import java.awt.image.PixelGrabber;
import java.text.AttributedCharacterIterator;


public class WMFGraphics extends Graphics
{
    public static final int MM_ANISOTROPIC = 8;
    private static final int BS_NULL = WMFConstants.BS_NULL;
    public static final int SRC_AND = WMFConstants.SRCAND;// 0x8800c6;
    public static final int SRC_PAINT = WMFConstants.SRCPAINT;// 0xee0086;
    public static final int SRC_COPY =  WMFConstants.SRCCOPY;// 0xcc0020;

    private MetaFile wmf;
    private Color foreground;
    private Color background;
    private Font font;
    private int penstyle;
    private int penwidth;
    private int brushfillstyle;
    private int brushhatch;
    private int fontescapement;
    private Image brushpattern;
    private int penhandle;
    private int brushhandle;
    private int fonthandle;
    private Rectangle clip;
    private Point origin;
    private GraphicsState state;


    public WMFGraphics(MetaFile w, int width, int height, Color fg, Color bg)
    {
        font = new Font("Helvetica", 0, 12);
        penstyle = 0;
        penwidth = 1;
        brushfillstyle = 0;
        brushhatch = 0;
        fontescapement = 0;
        brushpattern = null;
        origin = new Point(0, 0);
        foreground = fg;
        background = bg;
        state = new GraphicsState();
        state.increaseCount();
        setWMF(w, width, height);
        reset();
    }

    public Color getBackground()
    {
        return background;
    }

    public void setBackground(Color background)
    {
        this.background = background;
    }

    WMFGraphics(
        MetaFile w, GraphicsState graphicsstate, Color fg, Color bg, Font font, int penstyle, int penwidth,
        int brushFillStyle, int brushHatch, int fontEsc, Image image, Point origin, Rectangle rectangle)
    {
        this.font = new Font("Helvetica", 0, 12);
        this.wmf = w;
        this.state = graphicsstate;
        graphicsstate.increaseCount();
        this.foreground = fg;
        this.background = bg;
        this.font = font;
        this.penstyle = penstyle;
        this.penwidth = penwidth;
        this.brushfillstyle = brushFillStyle;
        this.brushhatch = brushHatch;
        this.fontescapement = fontEsc;
        this.brushpattern = image;
        this.origin = new Point(origin.x, origin.y);

        if (rectangle != null) {
            clip = new Rectangle(rectangle);
        }
        createHandles();
        setGDIPen();
        setGDIHollowBrush();
        setGDIFont();
    }

    public void GDIPolyPolygon(Polygon[] apolygon)
    {
        restore();
        setGDIFillBrush();
        wmf.polypolygon(apolygon);
        setGDIHollowBrush();
    }

    @Override
    public void clearRect(int i, int j, int k, int l)
    {
        restore();

        Color color = foreground;
        setColor(background);
        fillRect(i, j, k, l);
        setColor(color);
    }

    @Override
    public void clipRect(int i, int j, int k, int l)
    {
        restore();
        wmf.intersectClipRect(i, j, i + k + 1, j + l + 1);

        Rectangle rectangle = new Rectangle(i, j, k, l);

        if (clip != null) {
            clip = clip.intersection(rectangle);
        } else {
            clip = rectangle;
        }

        state.setClip(clip);
    }

    @Override
    public void copyArea(int i, int j, int k, int l, int i1, int j1)
    {
        System.err.println("copyArea not supported");
    }

    @Override
    public Graphics create()
    {
        WMFGraphics wmfgraphics = new WMFGraphics(
            wmf, state, foreground, background, font, penstyle, penwidth, brushfillstyle,
            brushhatch, fontescapement, brushpattern, origin, clip);

        return wmfgraphics;
    }

    public void createHandles()
    {
        penhandle = wmf.createPenIndirect(penstyle, penwidth, foreground);
        wmf.selectObject(penhandle);
        state.setPen(penhandle);
        brushhandle = wmf.createBrushIndirect(BS_NULL, foreground, brushhatch);
        wmf.selectObject(brushhandle);
        state.setPen(brushhandle);
        fonthandle = wmf.createFont(font, fontescapement, false, false);
        wmf.selectObject(fonthandle);
        state.setPen(fonthandle);
    }

    public void deleteHandles()
    {
        wmf.deleteObject(penhandle);
        wmf.deleteObject(brushhandle);
        wmf.deleteObject(fonthandle);
    }

    @Override
    public void dispose()
    {
        state.decreaseCount();
    }

    @Override
    public void drawArc(int x, int y, int width, int height, int startAngle, int arcAngle)
    {
        restore();

        int k1 = x + (width / 2);
        int l1 = y + (height / 2);
        wmf.arc(
            x, y, x + width + 1, y + height + 1,
            k1
                + (int) Math.round(
                (double) width * Math.sin((2 * Math.PI * (double) (startAngle + 90)) / 360D)),
            l1
                + (int) Math.round(
                (double) height * Math.cos((2 * Math.PI * (double) (startAngle + 90)) / 360D)),
            k1
                + (int) Math.round(
                (double) width * Math.sin((2 * Math.PI * (double) (startAngle + arcAngle + 90)) / 360D)),
            l1
                + (int) Math.round(
                (double) height * Math.cos((2 * Math.PI * (double) (startAngle + arcAngle + 90)) / 360D)));
    }

    @Override
    public boolean drawImage(
        Image image, int leftD, int topD, int rightD, int bottomD, int leftS, int topS, int rightS, int bottomS, Color color,
        ImageObserver imageobserver)
    {
        restore();

        int imagewidth = image.getWidth(imageobserver);
        int imageheight = image.getHeight(imageobserver);
        int[] pixbuffer = new int[imagewidth * imageheight];
        PixelGrabber pixelgrabber = new PixelGrabber(image, 0, 0, imagewidth, imageheight, pixbuffer, 0, imagewidth);

        try {
            pixelgrabber.grabPixels();
        } catch (InterruptedException _ex) {
            return false;
        }

        if ((pixelgrabber.status() & 0x80) != 0) {
            return false;
        }

        int diffx = rightD - leftD;
        int diffy = bottomD - topD;
        int sourcediffx = rightS - leftS;
        int sourcediffy = bottomS - topS;
        int k3 = bottomS;
        bottomS = imageheight - topS;
        topS = imageheight - k3;

        if ((diffx < 0) != (sourcediffx < 0)) {
            flipHorizontal(pixbuffer, imagewidth, imageheight);

            if (sourcediffx < 0) {
                leftS = imagewidth - leftS;
            } else {
                leftS = imagewidth - rightS;
            }
        }

        if (diffx < 0) {
            leftD = rightD;
            if (sourcediffx < 0) {
                leftS = rightS;
            }
            diffx = -diffx;
        }

        if (sourcediffx < 0) {
            sourcediffx = -sourcediffx;
        }

        if ((diffy < 0) != (sourcediffy < 0)) {
            flipVertical(pixbuffer, imagewidth, imageheight);

            if (sourcediffy < 0) {
                topS = imageheight - topS;
            } else {
                topS = imageheight - bottomS;
            }
        }

        if (diffy < 0) {
            topD = bottomD;

            if (sourcediffy < 0) {
                topS = bottomS;
            }

            diffy = -diffy;
        }

        if (sourcediffy < 0) {
            sourcediffy = -sourcediffy;
        }

        int l3 = color.getRGB();

        for (int i4 = 0; i4 < pixbuffer.length; i4++)
            if ((pixbuffer[i4] & 0xff000000) == 0) {
                pixbuffer[i4] = l3;
            }

        wmf.stretchBlt(leftD, topD, diffx, diffy, leftS, topS,
            sourcediffx, sourcediffy, SRC_COPY, pixbuffer, imagewidth, imageheight);

        return true;
    }

    @Override
    public boolean drawImage(Image image,
                             int leftD, int topD, int rightD, int bottomD, int leftS, int topS, int rightS, int bottomS,
                             ImageObserver imageobserver)
    {
        restore();

        int imagewidth = image.getWidth(imageobserver);
        int imageheight = image.getHeight(imageobserver);
        int[] pixarray = new int[imagewidth * imageheight];
        PixelGrabber pixelgrabber = new PixelGrabber(image, 0, 0, imagewidth, imageheight, pixarray, 0, imagewidth);

        try {
            pixelgrabber.grabPixels();
        } catch (InterruptedException _ex) {
            return false;
        }

        if ((pixelgrabber.status() & 0x80) != 0) {
            return false;
        }

        int ddiffx = rightD - leftD;
        int ddiffy = bottomD - topD;
        int sdiffx = rightS - leftS;
        int sdiffy = bottomS - topS;
        bottomS = imageheight - topS;
        topS = imageheight - bottomS;

        if ((ddiffx < 0) != (sdiffx < 0)) {
            flipHorizontal(pixarray, imagewidth, imageheight);

            if (sdiffx < 0) {
                leftS = imagewidth - leftS;
            } else {
                leftS = imagewidth - rightS;
            }
        }

        if (ddiffx < 0) {
            leftD = rightD;

            if (sdiffx < 0) {
                leftS = rightS;
            }

            ddiffx = -ddiffx;
        }

        if (sdiffx < 0) {
            sdiffx = -sdiffx;
        }

        if ((ddiffy < 0) != (sdiffy < 0)) {
            flipVertical(pixarray, imagewidth, imageheight);

            if (sdiffy < 0) {
                topS = imageheight - topS;
            } else {
                topS = imageheight - bottomS;
            }
        }

        if (ddiffy < 0) {
            topD = bottomD;

            if (sdiffy < 0) {
                topS = bottomS;
            }

            ddiffy = -ddiffy;
        }

        if (sdiffy < 0) {
            sdiffy = -sdiffy;
        }

        int[] ai1 = new int[pixarray.length];
        boolean flag = false;

        for (int l3 = 0; l3 < pixarray.length; l3++)
            if ((pixarray[l3] & 0xff000000) == 0) {
                ai1[l3] = -1;
                pixarray[l3] = 0;
                flag = true;
            } else {
                ai1[l3] = 0;
            }

        if (flag) {
            wmf.stretchBlt(leftD, topD, ddiffx, ddiffy, leftS, topS, sdiffx, sdiffy, SRC_AND, ai1, imagewidth, imageheight);
            wmf.stretchBlt(leftD, topD, ddiffx, ddiffy, leftS, topS, sdiffx, sdiffy, SRC_PAINT, pixarray, imagewidth, imageheight);
        } else {
            wmf.stretchBlt(leftD, topD, ddiffx, ddiffy, leftS, topS, sdiffx, sdiffy, SRC_COPY, pixarray, imagewidth, imageheight);
        }

        return true;
    }

    @Override
    public boolean drawImage(
        Image image, int left, int top, int width, int height, Color color, ImageObserver imageobserver)
    {
        restore();

        return drawImage(
            image, left, top, left + width, top + height, 0, 0, image.getWidth(imageobserver),
            image.getHeight(imageobserver), color, imageobserver);
    }

    @Override
    public boolean drawImage(Image image, int x, int y, int w, int h, ImageObserver imageobserver)
    {
        restore();

        return drawImage(
            image, x, y, x + w, y + h, 0, 0, image.getWidth(imageobserver),
            image.getHeight(imageobserver), imageobserver);
    }

    @Override
    public boolean drawImage(Image image, int x, int y, Color color, ImageObserver imageobserver)
    {
        restore();

        return drawImage(
            image, x, y, image.getWidth(imageobserver), image.getHeight(imageobserver), color,
            imageobserver);
    }

    @Override
    public boolean drawImage(Image image, int x, int y, ImageObserver imageobserver)
    {
        restore();

        return drawImage(
            image, x, y, image.getWidth(imageobserver), image.getHeight(imageobserver),
            imageobserver);
    }

    @Override
    public void drawLine(int x1, int y1, int x2, int y2)
    {
        restore();
        wmf.moveTo(x1, y1);
        wmf.lineTo(x2, y2);
        wmf.setPixel(x2, y2, getColor());
    }

    @Override
    public void drawOval(int x, int y, int width, int height)
    {
        restore();
        wmf.ellipse(x, y, x + width + 1, y + height + 1);
    }

    @Override
    public void drawPolygon(int[] ai, int[] ai1, int i)
    {
        restore();
        wmf.polygon(ai, ai1, i);
    }

    @Override
    public void drawPolyline(int[] xPoints, int[] yPoints, int i)
    {
        restore();
        wmf.polyline(xPoints, yPoints, i);
        wmf.setPixel(xPoints[i - 1], yPoints[i - 1], getColor());
    }

    @Override
    public void drawRect(int x, int y, int width, int height)
    {
        restore();
        wmf.rectangle(x, y, x + width + 1, y + height + 1);
    }

    @Override
    public void drawRoundRect(int x, int y, int width, int height, int i1, int j1)
    {
        restore();
        wmf.roundRect(x, y, x + width + 1, y + height + 1, i1, j1);
    }

    @Override
    public void drawString(String s, int x, int y)
    {
        restore();
        wmf.textOut(x, y, s);
    }

    @Override
    public void drawString(AttributedCharacterIterator attributedcharacteriterator, int i, int j)
    {
        System.err.println("drawString(java.text.AttributedCharacterIterator,..) not supported");
    }

    @Override
    public void fillArc(int x, int y, int width, int height, int startAngle, int endAngle)
    {
        restore();
        setGDIFillBrush();

        int cx = x + (width / 2);
        int cy = y + (height / 2);

        wmf.pie(
            x, y, x + width + 1, y + height + 1,
            cx
                + (int) Math.round(
                (double) width * Math.sin((Math.PI * 2 * (double) (startAngle + 90)) / 360D)),
            cy
                + (int) Math.round(
                (double) height * Math.cos((Math.PI * 2 * (double) (startAngle + 90)) / 360D)),
            cx
                + (int) Math.round(
                (double) width * Math.sin((Math.PI * 2 * (double) (startAngle + endAngle + 90)) / 360D)),
            cy
                + (int) Math.round(
                (double) height * Math.cos((Math.PI * 2 * (double) (startAngle + endAngle + 90)) / 360D)));
        setGDIHollowBrush();
    }

    @Override
    public void fillOval(int x, int y, int rx, int ry)
    {
        restore();
        setGDIFillBrush();
        drawOval(x, y, rx - 1, ry - 1);
        setGDIHollowBrush();
    }

    @Override
    public void fillPolygon(int[] ai, int[] ai1, int i)
    {
        restore();
        setGDIFillBrush();
        drawPolygon(ai, ai1, i);
        setGDIHollowBrush();
    }

    @Override
    public void fillRect(int i, int j, int k, int l)
    {
        restore();
        setGDIFillBrush();
        drawRect(i, j, k - 1, l - 1);
        setGDIHollowBrush();
    }

    @Override
    public void fillRoundRect(int i, int j, int k, int l, int i1, int j1)
    {
        restore();
        setGDIFillBrush();
        drawRoundRect(i, j, k - 1, l - 1, i1, j1);
        setGDIHollowBrush();
    }

    private void flipHorizontal(int[] ai, int i, int j)
    {
        for (int k = 0; k < j; k++) {
            int l = k * j;

            for (int i1 = 0; i1 < (i / 2); i1++) {
                int j1 = ai[l + i1];
                ai[l + i1] = ai[(l + i) - 1 - i1];
                ai[(l + i) - 1 - i1] = j1;
            }
        }
    }

    private void flipVertical(int[] ai, int i, int j)
    {
        int[] ai1 = new int[i];

        for (int k = 0; k < (j / 2); k++) {
            System.arraycopy(ai, k * i, ai1, 0, i);
            System.arraycopy(ai, (j - k - 1) * i, ai, k * i, i);
            System.arraycopy(ai1, 0, ai, (j - k - 1) * i, i);
        }
    }

    public int getBrushFillStyle()
    {
        return brushfillstyle;
    }

    public Shape getClip()
    {
        return getClipBounds();
    }

    public Rectangle getClipBounds()
    {
        if (clip != null) {
            return new Rectangle(clip);
        } else {
            return null;
        }
    }

    @Override
    public Color getColor()
    {
        return foreground;
    }

    @Override
    public Font getFont()
    {
        return font;
    }


    @Override
    @SuppressWarnings("deprecation")
    public FontMetrics getFontMetrics(Font font1)
    {
        return Toolkit.getDefaultToolkit().getFontMetrics(font1);
    }

    public int getPenStyle()
    {
        return penstyle;
    }


    private void reset()
    {
        setPenStyle(0);
        setPenWidth(1);
        setBrushFillStyle(0);
        setBrushHatch(0);
        setFontEscapement(0);
    }

    private void restore()
    {
        if (state.getCount() > 1) {
            if (penhandle != state.getPen()) {
                wmf.selectObject(penhandle);
                state.setPen(penhandle);
            }

            if (brushhandle != state.getBrush()) {
                wmf.selectObject(brushhandle);
                state.setBrush(brushhandle);
            }

            if (fonthandle != state.getFont()) {
                wmf.selectObject(fonthandle);
                state.setFont(fonthandle);
            }

            if (clip != state.getClip()) {
                setClip(clip);
                state.setClip(clip);
            }

            if (!origin.equals(state.getOrigin())) {
                translate(origin.x, origin.y);
                state.setOrigin(origin);
            }
        }
    }

    public void setBrushFillStyle(int i)
    {
        brushfillstyle = i;
    }

    public void setBrushHatch(int i)
    {
        brushhatch = i;
    }

    public void setBrushPattern(Image image)
    {
        brushpattern = image;
    }

    @Override
    public void setClip(int x, int y, int w, int h)
    {
        setClip(((new Rectangle(x, y, w, h))));
    }

    @Override
    public void setClip(Shape shape)
    {
        wmf.setClipRgn();
        clip = null;

        if (shape != null) {
            Rectangle rectangle = shape.getBounds();
            clipRect(rectangle.x, rectangle.y, rectangle.width, rectangle.height);
        }
    }

    @Override
    public void setColor(Color color)
    {
        restore();
        foreground = color;
        setGDIPen();
        wmf.setTextColor(foreground);
    }

    @Override
    public void setFont(Font font1)
    {
        restore();
        font = font1;
        setGDIFont();
    }

    public void setFontEscapement(int i)
    {
        fontescapement = i;
        setGDIFont();
    }

    public int setGDIFillBrush()
    {
        int i = brushhandle;

        if (brushfillstyle == 3) {
            if (brushpattern != null) {
                int j = brushpattern.getWidth(null);
                int k = brushpattern.getHeight(null);
                int[] ai = new int[j * k];
                PixelGrabber pixelgrabber = new PixelGrabber(brushpattern, 0, 0, j, k, ai, 0, j);

                try {
                    pixelgrabber.grabPixels();

                    if ((pixelgrabber.status() & 0x80) != 0) {
                        brushhandle = wmf.createBrushIndirect(0, foreground, brushhatch);
                    } else {
                        brushhandle = wmf.createPatternBrush(ai, j, k);
                    }
                } catch (InterruptedException _ex) {
                    brushhandle = wmf.createBrushIndirect(0, foreground, brushhatch);
                }
            } else {
                brushhandle = wmf.createBrushIndirect(0, foreground, brushhatch);
            }
        } else {
            brushhandle = wmf.createBrushIndirect(brushfillstyle, foreground, brushhatch);
        }

        wmf.selectObject(brushhandle);
        wmf.deleteObject(i);
        state.setBrush(brushhandle);

        return brushhandle;
    }

    public int setGDIFont()
    {
        int i = fonthandle;
        fonthandle = wmf.createFont(font, fontescapement, false, false);
        wmf.selectObject(fonthandle);
        wmf.deleteObject(i);
        state.setFont(fonthandle);

        return fonthandle;
    }

    public int setGDIHollowBrush()
    {
        int i = brushhandle;
        brushhandle = wmf.createBrushIndirect(1, foreground, brushhatch);
        wmf.selectObject(brushhandle);
        wmf.deleteObject(i);
        state.setBrush(brushhandle);

        return brushhandle;
    }

    public int setGDIPen()
    {
        int i = penhandle;
        penhandle = wmf.createPenIndirect(penstyle, penwidth, foreground);
        wmf.selectObject(penhandle);
        wmf.deleteObject(i);
        state.setPen(penhandle);

        return penhandle;
    }

    @Override
    public void setPaintMode()
    {
        System.err.println("setPaintMode not supported");
    }

    public void setPenStyle(int i)
    {
        penstyle = i;
        setGDIPen();
    }

    public void setPenWidth(int i)
    {
        penwidth = i;
        setGDIPen();
    }

    public void setWMF(MetaFile wmf1, int width, int height)
    {
        wmf = wmf1;
        createHandles();
        setup(width, height);
    }

    @Override
    public void setXORMode(Color color)
    {
        System.err.println("setXORMode not supported");
    }

    public void setSize(int width, int height)
    {
        wmf.setWindowExt(width, height);
    }

    private void setup(int width, int height)
    {
        wmf.setMapMode(MM_ANISOTROPIC);
        wmf.setWindowOrg(0, 0);
        wmf.setWindowExt(width, height);
        wmf.setViewportExt(width, height);
        wmf.setTextAlign(24);
        wmf.setBKMode(1);
        wmf.setBKColor(background);
        wmf.setTextColor(foreground);
        wmf.setPolyFillMode(1);
        wmf.setStretchBltMode(3);
        wmf.setROP2(13);
        wmf.setTextCharacterExtra(0);
    }

    @Override
    public void translate(int i, int j)
    {
        origin.translate(i, j);
        wmf.setWindowOrg(-origin.x, -origin.y);
        state.setOrigin(origin);
    }

    public void setPolyFillMode(int alternate)
    {
        wmf.setPolyFillMode(alternate);
    }
}
