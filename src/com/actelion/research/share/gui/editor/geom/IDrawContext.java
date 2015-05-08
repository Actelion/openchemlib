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

package com.actelion.research.share.gui.editor.geom;


/**
 * Project:
 * User: rufenec
 * Date: 11/24/2014
 * Time: 3:22 PM
 */
public interface IDrawContext<T>
{

    T getNative();

    void drawLine(double x, double y, double x1, double y1);

    void drawDashedLine(double srcx, double srcy, double targetx, double targety, int[] dashPattern);

    void drawPolygon(IPolygon polygon);

//    void drawArrow(double[] px, double[] py, boolean selected);

    java.awt.Dimension getBounds(String s);

    void setFont(String helvetica, double size,boolean bold);

    void setFill(long color);

    void fillText(String str, double x, double y);

    void save();

    void restore();

    void drawRect(double x, double y, double width, double height);

    void drawText(String s, double x, double y, boolean centerHorz, boolean centerVert);

    void clearRect(double x, double y, double w, double h);

    void setStroke(long color);

    void fillElipse(double v, double v1, double highlightAtomDiameter, double highlightAtomDiameter1);

    void fillRect(double x, double y, double w, double h);

    void strokeLine(double x, double y, double x1, double y1);

    void fillPolygon(double[] px, double[] py, int i);

    void setLineWidth(double i);


//    void save();
//
//    void setStroke(Color highlightStyle);
//
//    void setLineWidth(int i);
//
//    void setLineCap(StrokeLineCap defaultStrokeLineCap);
//
//    void setLineJoin(StrokeLineJoin defaultStrokeLineJoin);
//
//    void beginPath();
//
//    void moveTo(double x1, double y1);
//
//    void lineTo(double x2, double y2);
//
//    void stroke();
//
//    void closePath();
//
//    void restore();
//
//    void setFill(Color highlightStyle);
//
//    void setTextSize(int i);
//
//    void fillText(String s, double x, double y);
//
//    void fillArc(double x, double y, double dx, double dy, int startAngle, int endAngle, ArcType defaultArcType);
//
//    void strokeLine(double x, double y, double x1, double y1);
}
