/*
* @(#)NativeClipboardAccessor.java	1.2 17.10.2003
*
* Copyright 2000-2003 Actelion Ltd. All Rights Reserved.
*
* This software is the proprietary information of Actelion Ltd.
* Use is subject to license terms.
*
 */

package com.actelion.research.gui.clipboard;

/**
 *
 * <p>Title: Actelion Library</p>
 * <p>Description: Actelion Java Library</p>
 * <p>Copyright: Copyright (c) 2002-2003</p>
 * <p>Company: Actelion Ltd</p>
 * @author Thomas Sander, Christian Rufener
 * @version 1.2
 */
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
