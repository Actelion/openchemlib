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

package com.actelion.research.gui.swing;

import com.actelion.research.gui.LookAndFeelHelper;
import com.actelion.research.gui.generic.GenericCursorHelper;
import com.actelion.research.util.Platform;

import java.awt.*;

public class SwingCursorHelper extends GenericCursorHelper {

	private static Cursor[]	sCursor;
	private static boolean sIsDarkLaF;

	public static Cursor getCursor(int cursor) {
		boolean isDarkLaF = LookAndFeelHelper.isDarkLookAndFeel();
		if (sIsDarkLaF != isDarkLaF) {
			sIsDarkLaF = isDarkLaF;
			sCursor = null;
			}

		if (sCursor == null)
			sCursor = new Cursor[cCursorCount];

		if (sCursor[cursor] == null)
			sCursor[cursor] = createCursor(cursor);

		return sCursor[cursor];
		}

	public static Cursor createCursor(int cursor) {
		Toolkit tk = Toolkit.getDefaultToolkit();
		Dimension size = tk.getBestCursorSize(32, 32);

		if (!Platform.isWindows()) {
			/* if (size.width == 32 && size.height == 32 && IMAGE_NAME_32[cursor] != null) {
				SwingImage image = new SwingImage("cursor/" + IMAGE_NAME_32[cursor]);
				adaptForLaF(image);
				return tk.createCustomCursor(image.get(), new Point(HOTSPOT_32[2*cursor], HOTSPOT_32[2*cursor+1]), "");
				} Swing seems to use a very small color palette. It accepts 32x32 cursors, but they seem a little large */
			if (size.width >= 24 && size.height >= 24 && IMAGE_NAME_32[cursor] != null) {
				SwingImage image = new SwingImage("cursor/" + IMAGE_NAME_32[cursor]);
				adaptForLaF(image);
				image.scale(24, 24);
				return tk.createCustomCursor(image.get(), new Point(HOTSPOT_32[2*cursor]*3/4, HOTSPOT_32[2*cursor+1]*3/4), "");
				}
			}
		if (size.width>15 && size.height>15 && cursor<IMAGE_DATA_16.length) {
			SwingImage image = new SwingImage(size.width, size.height);
			build16x16CursorImage(image, cursor);
			adaptForLaF(image);
			return tk.createCustomCursor(image.get(), new Point(HOTSPOT_16[2*cursor], HOTSPOT_16[2*cursor+1]), "");
			}

		if (cursor == cPointerCursor)
			return Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR);
		if (cursor == cTextCursor)
			return Cursor.getPredefinedCursor(Cursor.TEXT_CURSOR);
		if (cursor == cPointedHandCursor)
			return Cursor.getPredefinedCursor(Cursor.HAND_CURSOR);

		return Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR);
		}
	} 
