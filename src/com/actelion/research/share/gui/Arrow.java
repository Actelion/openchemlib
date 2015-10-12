/*
 * Project: DD_jfx
 * @(#)Arrow.java
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

import java.awt.geom.Rectangle2D;

public class Arrow
{

     /*
                                      *
                                        **
                                          ***
                                            ****
                                              *****
                                                ******
     ***************************************************
                                                ******
                                              *****
                                            ****
                                          ***
                                        **
                                      *
      */


    protected Rectangle2D rect = null;

    public Arrow(double x, double y, double w, double h)
    {
        rect = new Rectangle2D.Double(x, y, w, h);
    }

//    public Arrow(Rectangle2D r)
//    {
//        rect = r;
//    }

    public void paint(IDrawContext ctx)
    {
        double dx = (rect.getMinX());
        double dy = (rect.getMinY());
        double dwidth = rect.getWidth();
        double dheight = rect.getHeight();
        double arrowEndX = dx + dwidth;
        double arrowEndY = dy + (dheight / 2);
        double xOffset = dwidth/15;
        double[] px = {
            (arrowEndX - (xOffset)),
            arrowEndX,
            (arrowEndX - (dwidth / 5))
        };
        double[] py = {
            arrowEndY,
            arrowEndY,
            (arrowEndY - (dwidth / 10))
        };
        if ((arrowEndX - dx) >= 5) {
            ctx.drawLine((int) dx, (int) arrowEndY, (int) arrowEndX-xOffset, (int) arrowEndY);
        }
        ctx.fillPolygon(px, py, 3);
        //ctx.drawLine(arrowEndX,arrowEndY,arrowEndX,arrowEndX-10);
        py[2] = (int) (arrowEndY + (dwidth / 10));
        ctx.fillPolygon(px, py, 3);
    }

}
