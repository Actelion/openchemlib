/*
 * Project: DD_jfx
 * @(#)TextClipboardHandler.java
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

import java.awt.Image;
import java.io.IOException;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.Toolkit;
import java.awt.datatransfer.DataFlavor;

/**
 * <p> </p>
 *
 * <p>Description: Actelion Electronic Lab Notebook</p>
 *
 * <p> </p>
 *
 * <p> </p>
 *
 * @author Christian Rufener
 * @version 4.0
 */
public class TextClipboardHandler
{

	public static String pasteText()
	{
		Clipboard systemClipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
		Transferable clipboardContents = systemClipboard.getContents(null);

		if (clipboardContents != null) {
			try {
				if (clipboardContents.isDataFlavorSupported(DataFlavor.stringFlavor)) {
					return (String)clipboardContents.getTransferData(DataFlavor.stringFlavor);
				}
			} catch (UnsupportedFlavorException ufe) {
				ufe.printStackTrace();
			} catch (IOException ioe) {
				ioe.printStackTrace();
			}
		}
		return null;
	}

	public static boolean copyText(String text)
	{
		boolean ok = false;
		try {
			TextSelection imageSelection = new TextSelection(text);
			Toolkit.getDefaultToolkit().getSystemClipboard().setContents(imageSelection,null);
			ok = true;
		} catch (Exception e) {

		}
		return ok;
	}

	private static class TextSelection implements Transferable
	{
		private String data;

		public TextSelection(String s) {
			data = s;
		}

		public DataFlavor[] getTransferDataFlavors() {
			return new DataFlavor[] { DataFlavor.stringFlavor };
		}

		public boolean isDataFlavorSupported(DataFlavor flavor) {
			return DataFlavor.stringFlavor.equals(flavor);
		}

		public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException,IOException {
			if (!DataFlavor.stringFlavor.equals(flavor)) {
				throw new UnsupportedFlavorException(flavor);
			}
			return data;
		}
	}

	private TextClipboardHandler()
	{
	}

}
