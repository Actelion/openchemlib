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

package com.actelion.research.share.gui.editor.chem;


import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.share.gui.editor.geom.IDrawContext;

/**
 * Project:
 * User: rufenec
 * Date: 5/16/13
 * Time: 12:31 PM
 */

public abstract class DrawingObject
        implements Cloneable,IDrawingObject
{

//    Object o;
    private boolean selected = false;

    public final boolean isSelected()
    {
        return selected;
    }
    public final void setSelected(boolean set)
    {
        selected = set;
    }

    private boolean layouted =false;
    public boolean isLayouted()
    {
        return layouted;
    }

    public void setLayouted(boolean layouted)
    {
        this.layouted = layouted;
    }

    public abstract void draw(IDrawContext ctx);
//    public abstract boolean pointWouldSelect(Point2D pt);
	public abstract GenericRectangle getBoundingRect();
    public abstract void move(float dx, float dy);
   	public abstract void scale(float f);

    public Object clone()
    {
        try {
            return super.clone();
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        return null;
    }
}

