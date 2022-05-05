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

package com.actelion.research.share.gui.editor.actions;

import com.actelion.research.gui.generic.GenericPoint;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.share.gui.editor.Model;
import com.actelion.research.share.gui.editor.chem.IArrow;
import com.actelion.research.share.gui.editor.geom.GeomFactory;
import com.actelion.research.share.gui.editor.geom.IDrawContext;
import com.actelion.research.share.gui.editor.io.IMouseEvent;

/**
 * Project:
 * User: rufenec
 * Date: 5/16/13
 * Time: 3:46 PM
 */
public class ArrowAction extends DrawAction
{
    GenericPoint origin,last;
    IArrow arrow = null;

    public ArrowAction(Model model)
    {
        super(model);
    }

    @Override
    public boolean onMouseDown(IMouseEvent ev)
    {
        origin = new GenericPoint(ev.getX(),ev.getY());
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean onMouseUp(IMouseEvent ev)
    {
        if (arrow != null) {
            model.addDrawingObject(arrow);
        }
        return true;
    }

    @Override
    public boolean onMouseMove(IMouseEvent ev, boolean drag)
    {
        if (drag) {
            GeomFactory factory = model.getGeomFactory();
            last = new GenericPoint(ev.getX(),ev.getY());
            GenericRectangle r = new GenericRectangle(
                (int) Math.min(last.getX(), origin.getX()),
                 (int) last.getY(),
                 (int) Math.abs(last.getX() - origin.getX()),0);
            arrow = factory.createArrow(r);
            return true;
        }
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean paint(IDrawContext ctx)
    {
//        if (arrow != null) {
//            arrow.draw(ctx);
//            return true;
//        }
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
