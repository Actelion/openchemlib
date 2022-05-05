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

package com.actelion.research.jfx.gui;

import com.actelion.research.gui.generic.GenericPoint;
import com.actelion.research.share.gui.editor.geom.IDrawContext;
import com.actelion.research.share.gui.editor.geom.IPolygon;
import javafx.geometry.Bounds;
import javafx.geometry.VPos;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.paint.Color;
import javafx.scene.shape.StrokeLineCap;
import javafx.scene.shape.StrokeLineJoin;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
import javafx.scene.text.TextAlignment;

/**
 * Project:
 * User: rufenec
 * Date: 11/24/2014
 * Time: 6:24 PM
 */
public class GraphicsContextImpl implements IDrawContext<GraphicsContext>
{
    GraphicsContext ctx;

    public GraphicsContextImpl(GraphicsContext graphicsContext2D)
    {
        this.ctx = graphicsContext2D;
    }


    public GraphicsContext getContext()
    {
        return ctx;
    }
    public static Color createColor(long color)
    {
        double r = (double) ((color & 0xFF000000l) >> 24) / 255.0;
        double g = (double) ((color & 0x00FF0000l) >> 16) / 255.0;
        double b = (double) ((color & 0x0000FF00l) >> 8) / 255.0;
        double a = (double) (color & 0x000000FFl ) / 255.0;
        return new Color(r, g, b, a);
    }

    @Override
    public GraphicsContext getNative()
    {
        return ctx;
    }

    @Override
    public void drawLine(double x, double y, double x1, double y1)
    {
        ctx.setLineCap(StrokeLineCap.ROUND);
        ctx.setLineJoin(StrokeLineJoin.MITER);
        ctx.beginPath();
        ctx.moveTo(x,y);
        ctx.lineTo(x1,y1);
        ctx.stroke();
    }

    @Override
    public void drawDashedLine(double srcx, double srcy, double targetx, double targety, int[] dashPattern)
    {
        double dx = targetx - srcx;
        double dy = targety - srcy;
        double angle = Math.atan2(dy, dx);
        double x = srcx;
        double y = srcy;
        int idx = 0;
        boolean draw = true;

        ctx.beginPath();
        ctx.moveTo(srcx, srcy);

        while (!((dx < 0 ? x <= targetx : x >= targetx) && (dy < 0 ? y <= targety : y >= targety))) {
            double dashLength = dashPattern[idx++ % dashPattern.length];
            double nx = x + (Math.cos(angle) * dashLength);
            x = dx < 0 ? Math.max(targetx, nx) : Math.min(targetx, nx);
            double ny = y + (Math.sin(angle) * dashLength);
            y = dy < 0 ? Math.max(targety, ny) : Math.min(targety, ny);
            if (draw) {
                ctx.lineTo(x, y);
            } else {
                ctx.moveTo(x, y);
            }
            draw = !draw;
        }

        ctx.closePath();
        ctx.stroke();

    }

    @Override
    public void drawPolygon(IPolygon polygon)
    {
        ctx.beginPath();
        GenericPoint pt = polygon.get(0);
        ctx.moveTo(pt.getX(), pt.getY());
        for (int i = 1; i < polygon.size(); i++) {
            pt = polygon.get(i);
            ctx.lineTo(pt.getX(), pt.getY());
        }
        ctx.closePath();
        ctx.stroke();
    }


    @Override
    public java.awt.Dimension getBounds(String s)
    {
        Font currentFont = ctx.getFont();
        Text t = new Text(s);
        t.setFont(currentFont);
        Bounds b = t.getLayoutBounds();
        java.awt.Dimension bounds = new java.awt.Dimension((int)b.getWidth(), (int)b.getHeight());
        return bounds;
    }

    @Override
    public void setFont(String name, double size, boolean bold)
    {
        ctx.setFont(Font.font(name,bold ? FontWeight.BOLD : null, size));
    }

    @Override
    public String getFont() {
        return ctx.getFont().getName();
    }

    @Override
    public void setFill(long color)
    {
        ctx.setFill(createColor(color));
    }

    @Override
    public void fillText(String str, double x, double y)
    {
        ctx.fillText(str,x,y);
    }

    @Override
    public void save()
    {
        ctx.save();
    }

    @Override
    public void restore()
    {
        ctx.restore();
    }

    @Override
    public void drawRect(double x, double y, double width, double height)
    {
        ctx.strokeRect(x,y,width,height);
    }

    @Override
    public void drawText(String s, double x, double y, boolean centerHorz, boolean centerVert)
    {
        ctx.setTextAlign(centerHorz ? TextAlignment.CENTER : TextAlignment.LEFT);
        ctx.setTextBaseline(centerVert ? VPos.CENTER : VPos.TOP);
        ctx.strokeText(s, x, y);

    }

    @Override
    public void clearRect(double x, double y, double w, double h)
    {
        ctx.clearRect(x, y, w, h);
    }

    @Override
    public void setStroke(long color)
    {
        ctx.setStroke(createColor(color));
    }

    @Override
    public void fillElipse(double x, double y, double rx, double ry)
    {
        ctx.fillArc(x, y, rx, ry, 0, 360, javafx.scene.shape.ArcType.ROUND);
    }

    @Override
    public void fillRect(double x, double y, double w, double h)
    {
        ctx.fillRect(x, y, w, h);
    }

    @Override
    public void strokeLine(double x, double y, double x1, double y1)
    {
        ctx.beginPath();
        ctx.moveTo(x,y);
        ctx.lineTo(x1,y1);
        ctx.stroke();
    }

    @Override
    public void fillPolygon(double[] px, double[] py, int i)
    {
        ctx.fillPolygon(px, py, i);
    }

    @Override
    public void setLineWidth(double i)
    {
        ctx.setLineWidth(i);

    }

}
