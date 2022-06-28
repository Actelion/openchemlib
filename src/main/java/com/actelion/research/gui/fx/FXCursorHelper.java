package com.actelion.research.gui.fx;

import com.actelion.research.gui.LookAndFeelHelper;
import com.actelion.research.gui.generic.GenericCursorHelper;
import com.actelion.research.util.Platform;
import javafx.geometry.Dimension2D;
import javafx.scene.Cursor;
import javafx.scene.ImageCursor;

public class FXCursorHelper extends GenericCursorHelper {
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
		Dimension2D size = ImageCursor.getBestSize(32, 32);
		if (!Platform.isWindows()) {
			/* if (size.getWidth() >= 32 && size.getHeight() >= 32 && IMAGE_NAME_32[cursor] != null) {
				FXImage image = new FXImage("cursor/" + IMAGE_NAME_32[cursor]);
				adaptForLaF(image);
				return new ImageCursor(image.get(), HOTSPOT_32[2*cursor], HOTSPOT_32[2*cursor+1]);
			} FX seems to prefer 24x24 cursors anyway */
			if (size.getWidth() >= 24 && size.getHeight() >= 24 && IMAGE_NAME_32[cursor] != null) {
				FXImage image = new FXImage("cursor/" + IMAGE_NAME_32[cursor]);
				adaptForLaF(image);
				image.scale(24, 24);
				return new ImageCursor(image.get(), HOTSPOT_32[2*cursor]*3/4, HOTSPOT_32[2*cursor+1]*3/4);
			}
		}

		if (size.getWidth() >= 16 && size.getHeight() >= 16 && cursor<IMAGE_DATA_16.length) {
			FXImage image = new FXImage((int)size.getWidth(), (int)size.getHeight());
			build16x16CursorImage(image, cursor);
			adaptForLaF(image);
			return new ImageCursor(image.get(), HOTSPOT_16[2*cursor], HOTSPOT_16[2*cursor+1]);
		}

		if (cursor == cPointerCursor)
			return Cursor.DEFAULT;
		if (cursor == cTextCursor)
			return Cursor.TEXT;
		if (cursor == cPointedHandCursor)
			return Cursor.HAND;

		return Cursor.DEFAULT;
	}
}
