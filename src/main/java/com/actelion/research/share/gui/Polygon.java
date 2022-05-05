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

import com.actelion.research.gui.generic.GenericPoint;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.share.gui.editor.geom.IPolygon;

import java.util.ArrayList;
import java.util.List;


public class Polygon implements IPolygon
{
    List<GenericPoint> list = new ArrayList<>();

    public void add(GenericPoint o)
    {
        list.add(o);
    }

    public void remove(GenericPoint o)
    {
        int idx = list.lastIndexOf(o);
        if (idx >= 0)
            list.remove(idx);
    }

    public int size()
    {
        return list.size();
    }

    public GenericPoint get(int i)
    {
        return list.get(i);
    }

    GenericRectangle calculateBounds()
    {
        double boundsMinX = Double.MAX_VALUE;
        double boundsMinY = Double.MIN_VALUE;
        double boundsMaxX = 0.0;
        double boundsMaxY = 0.0;

        for (int i = 0; i < list.size(); i++) {
            double x = list.get(i).getX();
            boundsMinX = Math.min(boundsMinX, x);
            boundsMaxX = Math.max(boundsMaxX, x);
            double y = list.get(i).getY();
            boundsMinY = Math.min(boundsMinY, y);
            boundsMaxY = Math.max(boundsMaxY, y);
        }
        GenericRectangle r = new GenericRectangle(boundsMinX, boundsMinY,
            boundsMaxX - boundsMinX,
            boundsMaxY - boundsMinY);
        return r;
    }

    public GenericRectangle getBoundingBox()
    {
        if (list.size() == 0) {
            return new GenericRectangle(0.0, 0.0, 0.0, 0.0);
        }
        return calculateBounds();
    }

    public boolean contains(GenericPoint pt)
    {
        return contains(pt.getX(), pt.getX());
    }

    public boolean contains(double x, double y)
    {
        GenericRectangle r = getBoundingBox();
        boolean contains = r.contains(x, y);
        if (list.size() <= 2 || !contains) {
            return false;
        }
        int hits = 0;
        double lastx = list.get(list.size() - 1).getX();
        double lasty = list.get(list.size() - 1).getY();
        double curx;
        double cury;

        // Walk the edges of the polygon
        for (int i = 0; i < list.size(); lastx = curx, lasty = cury, i++) {
            curx = list.get(i).getX();
            cury = list.get(i).getY();

            if (cury == lasty) {
                continue;
            }

            double leftx;
            if (curx < lastx) {
                if (x >= lastx) {
                    continue;
                }
                leftx = curx;
            } else {
                if (x >= curx) {
                    continue;
                }
                leftx = lastx;
            }

            double test1 = 0.0;
            double test2 = 0.0;
            if (cury < lasty) {
                if (y < cury || y >= lasty) {
                    continue;
                }
                if (x < leftx) {
                    hits++;
                    continue;
                }
                test1 = x - curx;
                test2 = y - cury;
            } else {
                if (y < lasty || y >= cury) {
                    continue;
                }
                if (x < leftx) {
                    hits++;
                    continue;
                }
                test1 = x - lastx;
                test2 = y - lasty;
            }

            if (test1 < (test2 / (lasty - cury) * (lastx - curx))) {
                hits++;
            }
        }
        return ((hits & 1) != 0);
    }

}
