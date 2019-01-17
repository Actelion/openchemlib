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
*/package com.actelion.research.jfx.gui.chem;import com.actelion.research.chem.AbstractDepictor;import com.actelion.research.chem.StereoMolecule;import com.actelion.research.jfx.dataformat.MoleculeDataFormats;import com.actelion.research.jfx.gui.misc.ClipboardHelper;import javafx.beans.property.BooleanProperty;import javafx.beans.property.SimpleBooleanProperty;import javafx.beans.property.StringProperty;import javafx.css.PseudoClass;import javafx.scene.control.Control;import javafx.scene.input.Clipboard;import javafx.scene.input.ClipboardContent;import javafx.scene.input.MouseEvent;import javafx.scene.paint.Color;/** * Project: * User: rufenec * Date: 10/12/11 * Time: 4:03 PM */public class MoleculeView extends Control implements IMoleculeView{    // TODO fix this hack!    // used to force class loader to load class MoleculeDataFormats    // (otherwise exception when class get loaded upon drop!!!)    static {        MoleculeDataFormats.DATA_FORMATS.toString();    }    private boolean transparent;    private boolean fillSize = false;    public MoleculeView()    {        this.setPrefSize(100, 100);        this.setSkin(new MoleculeViewSkin(this));    }    public MoleculeView(StereoMolecule mol)    {        this();        setValue(mol);    }    protected void layoutChildren()    {        super.layoutChildren();    }    @Override    public void setMolecule(StereoMolecule r)    {        ((MoleculeViewSkin)this.getSkin()).setMolecule(r);    }    public AbstractDepictor createDepictor(StereoMolecule value)    {        return new JFXCanvasDepictor(value);    }    public javafx.beans.property.ObjectProperty<StereoMolecule> valueProperty()    {        return ((MoleculeViewSkin)this.getSkin()).valueProperty();    }    public StringProperty idcodeProperty()    {        return ((MoleculeViewSkin)this.getSkin()).idcodeProperty();    }    public final void setValue(StereoMolecule t) {        ((MoleculeViewSkin)this.getSkin()).setValue(t);    }    public final StereoMolecule getValue() {        return ((MoleculeViewSkin)this.getSkin()).getValue();    }    public boolean isTransparent()    {        return transparent;    }    public void setTransparent(boolean transparent)    {        this.transparent = transparent;    }    protected ClipboardContent getContent(MouseEvent mouseEvent)    {        return ClipboardHelper.writeContent(getValue());    }    protected boolean putContent(Clipboard clipboard) {        StereoMolecule mol = ClipboardHelper.readContent(clipboard);        if (mol!=null) {            setValue(mol);            return true;        } else {            return false;        }    }    public Color getBackgroundColor()    {        return ((MoleculeViewSkin)this.getSkin()).getBackgroundColor();    }    public void setBackgroundColor(Color c)    {        ((MoleculeViewSkin)this.getSkin()).setBackgroundColor(c);    }    public boolean sizeContent()    {        return fillSize;    }    public void sizeContent(boolean sizeIt)    {        this.fillSize = sizeIt;    }    private static final PseudoClass PSEUDO_CLASS_READONLY            = PseudoClass.getPseudoClass("readonly");    private BooleanProperty editable = new SimpleBooleanProperty(this, "editable", true) {        @Override protected void invalidated() {            pseudoClassStateChanged(PSEUDO_CLASS_READONLY, ! get());        }    };    public final boolean isEditable() { return editable.getValue(); }    public final void setEditable(boolean value) { editable.setValue(value); }    public final BooleanProperty editableProperty() { return editable; }}