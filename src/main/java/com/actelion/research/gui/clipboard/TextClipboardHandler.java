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

package com.actelion.research.gui.clipboard;

import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.io.IOException;

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

		public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException {
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
