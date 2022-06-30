package com.actelion.research.gui.generic;

import com.actelion.research.gui.LookAndFeelHelper;

public abstract class GenericCursorHelper {
	public static final int cCursorCount = 14;

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
	public static final int cPointedHandCursor = 13;

	public static final int[][] IMAGE_DATA_16 = { {
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
	public static final int[] HOTSPOT_16 = {
			1, 7,
			7, 5,
			8, 7,
			8, 7,
			8, 7,
			2, 14,
			2, 14,
			0, 0,
			0, 0,
			8, 8,
			0, 0 };
	public static final int[] HOTSPOT_32 = {
			4, 20,
			10, 29,
			18, 18,
			18, 18,
			18, 18,
			16, 11,
			16, 11,
			4, 5,
			4, 5,
			16, 16,
			16, 16,
			0, 0,
			0, 0,
			12, 2 };

	public static final String[] IMAGE_NAME_32 = {
			"chain.png", "eraser.png", "hand.png", "handPlus.png", "fist.png", "lasso.png", "lassoPlus.png", "rect.png", "rectPlus.png", "zoom.png", "invisible.png", null, null, "pointingHand.png" };

	public static void adaptForLaF(GenericImage image) {
		if (LookAndFeelHelper.isDarkLookAndFeel())
			for (int x=0; x<image.getWidth(); x++)
				for (int y=0; y<image.getHeight(); y++)
					image.setRGB(x, y, image.getRGB(x, y) ^ 0x00FFFFFF);
		}

	public static void build16x16CursorImage(GenericImage image, int cursor) {
		if (IMAGE_DATA_16[cursor] != null) {
			for (int j = 0; j<IMAGE_DATA_16[cursor].length; j++) {
				int row = IMAGE_DATA_16[cursor][j];
				for (int k=15; k>=0; k--) {
					if ((row&3) == 0)
						image.setRGB(k, j, 0xFFFFFFFF);
					row >>= 2;
				}
			}
			for (int j = 0; j<IMAGE_DATA_16[cursor].length; j++) {
				int row = IMAGE_DATA_16[cursor][j];
				for (int k=15; k>=0; k--) {
					if ((row&3) == 1)
						image.setRGB(k, j, 0xFF000000);
					row >>= 2;
				}
			}
		}
	}
}
