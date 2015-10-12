/*
 * Project: DD_jfx
 * @(#)IPolygon.java
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
 * Time: 4:45 PM
 */
public interface IPolygon
{
    void add(java.awt.geom.Point2D pt);

    int size();

    java.awt.geom.Point2D get(int i);


    void remove(java.awt.geom.Point2D origin);

    boolean contains(double atomX, double atomY);
}
