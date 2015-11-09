/*
 * Project: DD_jfx
 * @(#)IDrawContext.java
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


}
