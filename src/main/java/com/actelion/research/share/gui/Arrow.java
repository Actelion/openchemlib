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

package com.actelion.research.share.gui;

import com.actelion.research.chem.DepictorTransformation;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.share.gui.editor.chem.IDrawingObject;
import com.actelion.research.share.gui.editor.geom.IDrawContext;

public class Arrow implements IDrawingObject
{

    protected final DrawConfig gfxConfig;/* = GeomFactory.getGeomFactory();*/

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


    protected GenericRectangle rect = null;
    private boolean selected = false;

    public Arrow(DrawConfig cfg, double x, double y, double w, double h)
    {
        gfxConfig = cfg;
        rect = new GenericRectangle(x, y, w, h);
    }


    @Override
    public void setSelected(boolean b) {
        selected = b;
    }

    @Override
    public boolean isSelected()
    {
        return selected;
    }

    @Override
    public void move(float dx, float dy) {
        rect.set(rect.getX()+dx,rect.getY()+dy,rect.getWidth(),rect.getHeight());
    }

    @Override
    public GenericRectangle getBoundingRect() {
        return rect;
    }

    @Override
    public void setRect(float x, float y, float w, float h) {
        rect = new GenericRectangle(x, y, w, h);
    }

    @Override
    public void scale(float scaling) {
        rect.set(rect.getX()*scaling,rect.getY()*scaling,rect.getWidth()*scaling,rect.getHeight()*scaling);
    }

    @Override
    public void draw(IDrawContext ctx,DepictorTransformation t)
    {

        double dx =        t == null ? (rect.getX()) :   t.transformX((float)rect.getX()) ;
        double dy =        t == null ? (rect.getY()) :   t.transformY((float)rect.getY()) ;
        double dwidth =    t == null ? rect.getWidth() :    t.transformX((float)rect.getWidth());
        double dheight =   t == null ? rect.getHeight() :   t.transformY((float)rect.getHeight()) ;
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

        if (selected) {
            ctx.setStroke(gfxConfig.getSelectionColor());
            ctx.setFill(gfxConfig.getSelectionColor());
        }
        if ((arrowEndX - dx) >= 5) {
            ctx.drawLine((int) dx, (int) arrowEndY, (int) arrowEndX-xOffset, (int) arrowEndY);
        }
        ctx.fillPolygon(px, py, 3);
        //ctx.drawLine(arrowEndX,arrowEndY,arrowEndX,arrowEndX-10);
        py[2] = (int) (arrowEndY + (dwidth / 10));
        ctx.fillPolygon(px, py, 3);

//        System.out.printf("Drawing Arrow %s\n",selected);
/*
        if (selected) {
            int width = 5;
            ctx.save();
            ctx.setLineWidth(width);
            ctx.setStroke(builder.getHighLightColor());
            ctx.drawLine((int) dx, (int) arrowEndY, (int) (arrowEndX - (dwidth / 5)), (int) arrowEndY);
            ctx.restore();
        }
*/
    }

    @Override
    public boolean isMovable()
    {
        return false;
    }

    public boolean isOnProductSide(float x, float y) {
//        System.out.printf("ISOnProductSide %f %f %s\n",x,y,rect);
        return x > rect.getX() + rect.getWidth()/2;
/*
        double dx = rect.getWidth();
        double dy = rect.getHeight();
        double mx = (rect.getX() + rect.getWidth()) / 2.0f;
        double my = (rect.getY() + rect.getHeight()) / 2.0f;

        if (dx == 0.0)
            return (dy < 0.0) ^ (y > my);

        if (dy == 0.0)
            return (dx < 0.0) ^ (x > mx);

        double m = -dx/dy;	// m of orthogonal line through S

        double sx = (rect.getX() + m*m*x - m*y + m*rect.getY()) / (1 + m*m);
        return (dx < 0.0) ^ (sx > mx);
*/
    }

}
