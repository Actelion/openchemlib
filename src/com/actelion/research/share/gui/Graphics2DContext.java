/*
 * Project: Mercury.5
 * @(#)Graphic2DContext.java
 *
 * Copyright (c) 1997- 2015
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All Rights Reserved.
 *
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 *
 * Author: Christian Rufener
 */

package com.actelion.research.share.gui;

import com.actelion.research.share.gui.editor.geom.IDrawContext;
import com.actelion.research.share.gui.editor.geom.IPolygon;

import java.awt.*;

/**
 * Created by rufenec on 01/04/15.
 */
public class Graphics2DContext implements IDrawContext<Graphics2D>
{
    BasicStroke stroke = new BasicStroke((float) 1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
    Point pt = null;
    Graphics2D delegate;
    public Graphics2DContext(Graphics2D graphics)
    {
        delegate= graphics;
    }


    public static Color createColor(long color)
    {
        int r = (int) ((color & 0xFF000000l) >> 24);
        int g = (int) ((color & 0x00FF0000l) >> 16);
        int b = (int) ((color & 0x0000FF00l) >> 8);
        int a = (int) (color & 0x000000FFl);
        return new Color(r, g, b, a);
    }

    public void save()
    {

    }

    public void setStroke(long color)
    {
        delegate.setColor(createColor(color));

    }

    @Override
    public void fillElipse(double v, double v1, double highlightAtomDiameter, double highlightAtomDiameter1)
    {
        delegate.fillArc((int) v, (int) v1, (int) highlightAtomDiameter, (int) highlightAtomDiameter1, 0, 360);

    }

/*
    public void setLineWidth(int lineWidth)
    {
        stroke = new BasicStroke((float) lineWidth, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
        delegate.setStroke(stroke);
    }
*/

//    public void setLineCap(StrokeLineCap lineCap)
//    {
//    }
//
//    public void setLineJoin(StrokeLineJoin round)
//    {
//
//    }

    private void beginPath()
    {

    }

    private void moveTo(double x1, double y1)
    {
        pt = new Point((int) x1, (int) y1);
    }

    private void stroke()
    {

    }

    private void lineTo(double x2, double y2)
    {
        delegate.drawLine(pt.x, pt.y, (int) x2, (int) y2);
        pt = new Point((int) x2, (int) y2);
    }

    private void closePath()
    {
    }

    @Override
    public void restore()
    {
    }

    @Override
    public void drawRect(double x, double y, double width, double height)
    {
        delegate.drawRect((int) x, (int) y, (int) width, (int) height);

    }

    @Override
    public void drawText(String s, double x, double y, boolean centerHorz, boolean centerVert)
    {
        delegate.drawString(s, (float)x, (float)y);

    }

    @Override
    public void clearRect(double x, double y, double w, double h)
    {
        delegate.clearRect((int)x,(int)y,(int)w,(int)h);

    }

    @Override
    public Graphics2D getNative()
    {
        return delegate;
    }

    @Override
    public void drawLine(double x, double y, double x1, double y1)
    {
        delegate.draw(new java.awt.geom.Line2D.Float((int) x, (int) y, (int) x1, (int) y1));
    }

    @Override
    public void drawDashedLine(double srcx, double srcy, double targetx, double targety, int[] dashPattern)
    {
        Stroke dashed = new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[]{2}, 0);
        delegate.setStroke(dashed);
        delegate.drawLine((int)srcx, (int)srcy, (int)targetx, (int)targety);

    }

    @Override
    public void drawPolygon(IPolygon polygon)
    {
        beginPath();
        java.awt.geom.Point2D pt = polygon.get(0);
        moveTo(pt.getX(), pt.getY());
        for (int i = 1; i < polygon.size(); i++) {
            pt = polygon.get(i);
            lineTo(pt.getX(), pt.getY());
        }
        closePath();
        stroke();
    }

    @Override
    public java.awt.Dimension getBounds(String s)
    {
        double w = getStringWidth(s);
        double h = delegate.getFontMetrics().getHeight();
        return new Dimension((int)w,(int)h);
    }

    protected double getStringWidth(String theString)
    {
        java.awt.font.GlyphVector mCurrentGlyphVector =
                (delegate).getFont().createGlyphVector((delegate).getFontRenderContext(), theString);
        return mCurrentGlyphVector.getLogicalBounds().getWidth();
    }

    @Override
    public void setFont(String fontName, double size,boolean bold)
    {
        delegate.setFont(new Font(fontName,bold ? Font.BOLD : 0,(int)size));

    }

    @Override
    public void setFill(long color)
    {
        delegate.setColor(createColor(color));
    }

    public void strokeLine(double x, double y, double x1, double y1)
    {
        drawLine(x, y, x1,y1);
    }

    public void fillPolygon(double[] px, double[] py, int count)
    {
        int x[] = new int[px.length];
        int y[] = new int[py.length];
        for (int i = 0; i < count; i++) {
            x[i] = (int) px[i];
            y[i] = (int) py[i];
        }
        delegate.fillPolygon(x, y, count);
    }

    @Override
    public void setLineWidth(double i)
    {
        stroke = new BasicStroke((float) i, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
        delegate.setStroke(stroke);
    }

    @Override
    public void fillRect(double x, double y, double width, double height)
    {
        delegate.fillRect((int) x, (int) y, (int) width, (int) height);
    }



    public void drawImage(Image img, double dx1, double dy1, double dx2, double dy2,
                          double sx1, double sy1, double sx2, double sy2)
    {
        delegate.drawImage(img,
                (int) dx1, (int) dy1, (int) dx2, (int) dy2,
                (int) sx1, (int) sy1, (int) sx2, (int) sy2, null);
    }

    public void drawImage(Image image, double x, double y)
    {
        delegate.drawImage(image, (int) x, (int) y, null);
    }


    public void fillText(String s, double x, double y)
    {
        delegate.drawString(s, (float) x, (float) y);
    }

    public void setTextSize(int textSize)
    {
        Font f = delegate.getFont();
        Font ft = new Font(f.getName(), 0, textSize);
        delegate.setFont(ft);
    }

    public void setColor(Color fillColor)
    {
        delegate.setColor(fillColor);
    }
}
