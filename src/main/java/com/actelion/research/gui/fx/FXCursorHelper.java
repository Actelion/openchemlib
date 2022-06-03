package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.GenericCursorHelper;
import javafx.geometry.Dimension2D;
import javafx.scene.Cursor;
import javafx.scene.ImageCursor;

public class FXCursorHelper extends GenericCursorHelper {
	private static Cursor[]	sCursor;

	public static Cursor getCursor(int cursor) {
		if (sCursor == null)
			sCursor = new Cursor[cCursorCount];

		if (sCursor[cursor] == null)
			sCursor[cursor] = createCursor(cursor);

		return sCursor[cursor];
	}

	public static Cursor createCursor(int cursor) {
		if (cursor<IMAGE_DATA_16.length) {
			Dimension2D size = ImageCursor.getBestSize(16, 16);
			if (size.getWidth() == 0 || size.getHeight() < 16)
				return Cursor.DEFAULT;

			FXImage image = new FXImage((int)size.getWidth(), (int)size.getHeight());
			build16x16CursorImage(image, cursor);
			return new ImageCursor(image.get());
		}
		else {
			switch (cursor) {
			case cPointerCursor:
				return Cursor.DEFAULT;
			case cTextCursor:
				return Cursor.TEXT;
			case cPointedHandCursor:
				return Cursor.HAND;
			}
		}
		return null;

//		Dimension2D size = ImageCursor.getBestSize(32, 32);
//		if (size.getWidth() == 0 || size.getHeight() == 0)
//			return Cursor.DEFAULT;
//
//		String fileName = cImageName[cursor];
//		if (fileName == null)
//			return Cursor.DEFAULT;
//
//		return new ImageCursor(new FXImage("cursor/"+fileName).get());
	}
}
