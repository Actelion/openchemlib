/*
 * Project: DD_jfx
 * @(#)Delegator.java
 *
 * Copyright (c) 1997- 2014
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All Rights Reserved.
 *
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 *
 * Author:
 */

package com.actelion.research.share.gui;

public class Delegator<T>
{
    T delegate;

    public Delegator(T t)
    {
        delegate = t;
    }
    public T getNative()
    {
        return delegate;
    }
}
