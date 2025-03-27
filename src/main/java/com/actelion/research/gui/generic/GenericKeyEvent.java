package com.actelion.research.gui.generic;

public class GenericKeyEvent extends GenericEvent {
	public static final int KEY_CTRL = -1;
	public static final int KEY_ALT = -2;
	public static final int KEY_SHIFT = -3;
	public static final int KEY_DELETE = -4;
	public static final int KEY_BACK_SPACE = -5;
	public static final int KEY_HELP = -6;
	public static final int KEY_ESCAPE = -7;
	public static final int KEY_ENTER = -8;

	public static final int KEY_PRESSED = 1;
	public static final int KEY_RELEASED = 2;

	private int mKey;
	private boolean mIsAltDown,mIsCtrlDown,mIsShiftDown,mIsMenuShortcut;

	public GenericKeyEvent(int what, int key, boolean isAltDown, boolean isCtrlDown, boolean isShiftDown, boolean isMenuShortcut, Object source) {
		super(what, source);
		mKey = key;
		mIsAltDown = isAltDown;
		mIsCtrlDown = isCtrlDown;
		mIsShiftDown = isShiftDown;
		mIsMenuShortcut = isMenuShortcut;
	}

	public int getKey() {
		return mKey;
	}

	public boolean isCtrlDown() {
		return mIsCtrlDown;
	}

	public boolean isAltDown() {
		return mIsAltDown;
	}

	public boolean isShiftDown() {
		return mIsShiftDown;
	}

	public boolean isMenuShortcut() {
		return mIsMenuShortcut;
	}
}
