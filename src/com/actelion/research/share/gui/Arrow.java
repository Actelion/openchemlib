/*
* Copyright (c) 1997 - 2015
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
