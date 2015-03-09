/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.util;

import java.awt.*;
import java.awt.image.*;
// Should this go into gui?

public class CursorHelper {
	public static final int cChainCursor = 0;
	public static final int cDeleteCursor = 1;
	public static final int cHandCursor = 2;
	public static final int cHandPlusCursor = 3;
	public static final int cFistCursor = 4;
	public static final int cLassoCursor = 5;
	public static final int cLassoPlusCursor = 6;
	public static final int cSelectRectCursor = 7;
	public static final int cSelectRectPlusCursor = 8;
	public static final int cZoomCursor = 9;
    public static final int cInvisibleCursor = 10;
	public static final int cPointerCursor = 11;
	public static final int cTextCursor = 12;

	private static final int cCursorCount = 13;

	private static final int[][] cCursorData = { {
		0xaaa2aaaa, 0xaa84aaaa, 0xaa152aa2, 0xa8514a84, 0xa1485214,
		0x852a1452, 0x54aa854a, 0x56aaa12a, 0x56aaa8aa },{
		0xaaaa5555, 0xaaa90005, 0xaaa40011, 0xaa900041, 0xaa400101,
		0xa9000406, 0xa400101a, 0x9000406a, 0x555501aa, 0x400106aa,
		0x40011aaa, 0x40016aaa, 0x5555aaaa },{
		0xaaa96aaa, 0xa96419aa, 0xa414146a, 0xa4141066, 0xa9041051,
		0xa9041041, 0x96400041, 0x41400001, 0x40400006, 0x90000006,
		0xa4000006, 0xa400001a, 0xa900001a, 0xaa40006a, 0xaa90006a,
		0xaa90006a },{
		0xaaa96aaa, 0xa96419aa, 0xa414146a, 0xa4141066, 0xa9041051,
		0xa9041041, 0x96400041, 0x41400001, 0x40404006, 0x90004006,
		0xa4055406, 0xa400401a, 0xa900401a, 0xaa40006a, 0xaa90006a,
		0xaa90006a },{
		0xaaaaaaaa, 0xaaaaaaaa, 0xaaa95aaa, 0xaa54156a, 0xa904105a,
		0xa9041046, 0x96400041, 0x91400001, 0x90400006, 0x90000006,
		0xa4000006, 0xa400001a, 0xa900001a, 0xaa40006a, 0xaa90006a,
		0xaa90006a, },{
		0xaaaaaaaa, 0xaaa5556a, 0xa9500016, 0xa4000001, 0x90000001,
		0x40000001, 0x40000016, 0x4000056a, 0x95015aaa, 0x4056aaaa,
		0x446aaaaa, 0x95aaaaaa, 0xa9aaaaaa, 0xa9aaaaaa, 0xa6aaaaaa },{
		0xaaaaaaaa, 0xaaa5556a, 0xa9500016, 0xa4010001, 0x90010001,
		0x40155001, 0x40010016, 0x4001056a, 0x95015aaa, 0x4056aaaa,
		0x446aaaaa, 0x95aaaaaa, 0xa9aaaaaa, 0xa9aaaaaa, 0xa6aaaaaa },{
		0x5555556a, 0x4000002a, 0x4000006a, 0x4000002a, 0x4000006a,
		0x4000002a, 0x4000006a, 0x4000002a, 0x4000006a, 0x4000002a,
		0x4000006a, 0x4000002a, 0x4444446a },{
		0x5555556a, 0x4000002a, 0x4000006a, 0x4000002a, 0x4000006a,
		0x4000002a, 0x4000406a, 0x4000402a, 0x4005546a, 0x4000402a,
		0x4000406a, 0x4000002a, 0x4444446a },{
		0xaaa55aaa, 0xaa5555aa, 0x695aa56a, 0x55aaaa5a, 0x55aaaa56,
		0x556aaa96, 0x556a6aa5, 0xaaaa6aa5, 0xaaa556a5, 0x556a6aa5,
		0x556a6a96, 0x55aaaa56, 0x55aaaa5a, 0x695aa56a, 0xaa5555aa,
		0xaaa55aaa, },
	    null };
	private static final Point[] cCursorHotSpot = {
		new Point( 1, 7),
		new Point( 7, 5),
		new Point( 8, 7),
		new Point( 8, 7),
		new Point( 8, 7),
		new Point( 2, 14),
		new Point( 2, 14),
		new Point( 0, 0),
        new Point( 0, 0),
        new Point( 8, 8),
		new Point( 0, 0) };

	private static Cursor[]	sCursor;

	public static Cursor getCursor(int cursor) {
		if (sCursor == null)
			sCursor = new Cursor[cCursorCount];

		if (sCursor[cursor] == null)
			sCursor[cursor] = createCursor(cursor);

		return sCursor[cursor];
		}

	private static Cursor createCursor(int cursor) {
		if (cursor<cCursorData.length) {
			Toolkit tk = Toolkit.getDefaultToolkit();
			Dimension size = tk.getBestCursorSize(16, 16);
			if (size != null && size.width>15 && size.height>15) {
				BufferedImage image = new BufferedImage(size.width, size.height, BufferedImage.TYPE_INT_ARGB);
                if (cCursorData[cursor] != null) {
    				Graphics2D g2 = image.createGraphics();
    				g2.setColor(new Color(255, 255, 255, 255));	// white, opaque
    				for (int j=0; j<cCursorData[cursor].length; j++) {
    					int row = cCursorData[cursor][j];
    					for (int k=15; k>=0; k--) {
    						if ((row&3) == 0)
    							g2.fillRect(k, j, 1, 1);
    						row >>= 2;
    						}
    					}
    				g2.setColor(new Color(0, 0, 0, 255));		// black, opaque
    				for (int j=0; j<cCursorData[cursor].length; j++) {
    					int row = cCursorData[cursor][j];
    					for (int k=15; k>=0; k--) {
    						if ((row&3) == 1)
    							g2.fillRect(k, j, 1, 1);
    						row >>= 2;
    						}
    					}
                    }
				return tk.createCustomCursor(image, cCursorHotSpot[cursor], "");
				}
			}
		else {
			switch (cursor) {
			case cPointerCursor:
				return Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR);
			case cTextCursor:
				return Cursor.getPredefinedCursor(Cursor.TEXT_CURSOR);
				}
			}
		return null;
		}
	} 
