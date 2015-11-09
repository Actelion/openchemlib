/*
 * Project: DD_jfx
 * @(#)NativeClipboardAccessor.java
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

package com.actelion.research.gui.clipboard;

public class NativeClipboardAccessor
{

    public static native boolean copyMoleculeToClipboard(String filname,byte[] sketch, byte[] serializedObject);
	public static native boolean copySizedMoleculeToClipboard(String filname, byte[] sketch, byte[] serializedObject, int cx,int cy);
	/* Formats are "MDLSK","MDLCT","MDL_MOL","CF_METAFILEPICT","CF_DIB" "ACT_MOLECULE" */
    public static native byte[] getClipboardData(String format);
    public static native boolean setClipBoardData(String format, byte[] buffer);

	public static boolean copyReactionToClipboard(String filname,byte[] sketch, byte[] serializedObject)
	{
		return copyMoleculeToClipboard(filname,sketch,serializedObject);
	}

    static {
        try {
            System.loadLibrary("actelionclip");
            System.out.println("actelionclip loaded");
        } catch (UnsatisfiedLinkError e) {
            e.printStackTrace();
        } catch (SecurityException e) {
        	e.printStackTrace();
        }
    }
}
