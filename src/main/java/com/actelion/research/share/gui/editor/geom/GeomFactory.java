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

package com.actelion.research.share.gui.editor.geom;


import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.share.gui.DrawConfig;
import com.actelion.research.share.gui.Polygon;
import com.actelion.research.share.gui.editor.chem.IArrow;
import com.actelion.research.share.gui.editor.dialogs.IAtomPropertiesDialog;
import com.actelion.research.share.gui.editor.dialogs.IAtomQueryFeaturesDialog;
import com.actelion.research.share.gui.editor.dialogs.IBondQueryFeaturesDialog;
import com.actelion.research.share.gui.editor.io.IKeyCode;

/**
 * Project:
 * User: rufenec
 * Date: 11/24/2014
 * Time: 3:15 PM
 */
public abstract class GeomFactory
{

    private DrawConfig drawConfig;
    public GeomFactory(DrawConfig cfg)
    {
        drawConfig = cfg;
    }


    public final IPolygon createPolygon()
    {
        return new Polygon();
    }
    public abstract IArrow createArrow(GenericRectangle r);
    public abstract IAtomQueryFeaturesDialog createAtomQueryFeatureDialog(StereoMolecule mol, int atom, boolean includeReactionHints);
    public abstract IBondQueryFeaturesDialog createBondFeaturesDialog(StereoMolecule mol, int bond);
    public abstract IAtomPropertiesDialog createAtomPropertiesDialog(StereoMolecule m, int atom);
    public abstract GenericRectangle getBoundingRect(StereoMolecule m);
    public abstract IKeyCode getDeleteKey();
    public abstract IKeyCode getEscapeKey();
    public abstract IKeyCode getBackSpaceKey();
    public abstract IKeyCode getEnterKey();

    public DrawConfig getDrawConfig()
    {
        return drawConfig;
    }

    public final long getHighLightColor() { return drawConfig.getHighLightColor();}
    public final long getMapToolColor() { return drawConfig.getMapToolColor();}
    public final long getSelectionColor() { return drawConfig.getSelectionColor();}
    public final long getForegroundColor() { return drawConfig.getForegroundColor();}
    public final long getBackgroundColor() { return drawConfig.getBackgroundColor();}

}
