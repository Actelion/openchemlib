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

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.share.gui.editor.Model;
import com.actelion.research.share.gui.editor.geom.ICursor;
import com.actelion.research.share.gui.editor.geom.IDrawContext;
import com.actelion.research.share.gui.editor.io.IMouseEvent;

import java.awt.geom.Point2D;

/**
 * Project:
 * User: rufenec
 * Date: 4/28/2014
 * Time: 12:39 PM
 */
public class ZoomRotateAction extends DrawAction {
    private java.awt.geom.Point2D origin = null;

    public ZoomRotateAction(Model m) {
        super(m);
    }

    @Override
    public void onActionEnter() {
        model.pushUndo();
    }

    @Override
    public boolean onMouseDown(IMouseEvent ev) {
        origin = new Point2D.Double(ev.getX(), ev.getY());
        StereoMolecule mol = model.getMolecule();
        mol.zoomAndRotateInit((float) origin.getX(), (float) origin.getY());
        return false;
    }

    @Override
    public boolean onMouseUp(IMouseEvent ev) {
        return true;
    }

    @Override
    public boolean onMouseMove(IMouseEvent ev, boolean drag) {
        if (drag) {
            boolean selectedOnly = false;
            StereoMolecule mol = model.getMolecule();
            for (int i = 0; i < mol.getAllAtoms(); i++) {
                if (mol.isSelectedAtom(i))
                    selectedOnly = true;
            }
            java.awt.geom.Point2D pt = new Point2D.Double(ev.getX(), ev.getY());
            float magnification = (Math.abs(pt.getY() - origin.getY()) < 20 ? 1.0f : (float) Math.exp((pt.getY() - origin.getY()) / 100f));
            float angleChange = (Math.abs(pt.getX() - origin.getX()) < 20 ? 0.0f : (float) (pt.getX() - origin.getX()) / 50.0f);
            mol.zoomAndRotate(magnification, angleChange, selectedOnly);
            return true;
        }
        return false;
    }

    @Override
    public boolean paint(IDrawContext ctx) {
        return false;
    }

    @Override
    public int getCursor() {
        return ICursor.TOOL_ZOOMCURSOR;
    }


}
