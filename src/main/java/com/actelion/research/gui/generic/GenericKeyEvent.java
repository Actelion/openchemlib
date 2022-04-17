package com.actelion.research.gui.generic;

public class GenericKeyEvent {
	public static final int KEY_CODE_CTRL = 1;
	public static final int KEY_CODE_ALT = 2;
	public static final int KEY_CODE_SHIFT = 3;

	public static final int KEY_PRESSED = 1;
	public static final int KEY_RELEASED = 2;
	public static final int KEY_TYPED = 3;

	private int mWhat;
	private int mKeyCode;
	private char mKeyChar;
	private boolean mIsMenuShortcut;

	public GenericKeyEvent(int what, int keyCode, char keyChar, boolean isMenuShortcut) {
		mWhat = what;
		mKeyCode = keyCode;
		mKeyChar = keyChar;
		mIsMenuShortcut = isMenuShortcut;
	}

	public int getWhat() {
		return mWhat;
	}

	public int getKeyCode() {
		return mKeyCode;
	}

	public char getKeyChar() {
		return mKeyChar;
	}

	public boolean isIsMenuShortcut() {
		return mIsMenuShortcut;
	}
}
