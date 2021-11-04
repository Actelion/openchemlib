package com.actelion.research.gui.generic;

public class GenericMouseEvent {
	public static final int MOUSE_PRESSED = 1;
	public static final int MOUSE_RELEASED = 2;
	public static final int MOUSE_CLICKED = 3;
	public static final int MOUSE_ENTERED = 4;
	public static final int MOUSE_EXITED = 5;
	public static final int MOUSE_MOVED = 6;
	public static final int MOUSE_DRAGGED = 7;

	private int mWhat,mButton,mClickCount,mX,mY;
	private boolean mShiftDown,mControlDown,mAltDown,mIsPopupTrigger;

	public GenericMouseEvent(int what, int button, int clickCount, int x, int y, boolean shiftDown, boolean controlDown, boolean altDown, boolean isPopupTrigger) {
		mWhat = what;
		mButton = button;
		mClickCount = clickCount;
		mX = x;
		mY = y;
		mShiftDown = shiftDown;
		mControlDown = controlDown;
		mAltDown = altDown;
		mIsPopupTrigger = isPopupTrigger;
	}

	public int getWhat() {
		return mWhat;
	}

	public int getButton() {
		return mButton;
	}

	public int getClickCount() {
		return mClickCount;
	}

	public int getX() {
		return mX;
	}

	public int getY() {
		return mY;
	}

	public boolean isShiftDown() {
		return mShiftDown;
	}

	public boolean isControlDown() {
		return mControlDown;
	}

	public boolean isAltDown() {
		return mAltDown;
	}

	public boolean isPopupTrigger() {
		return mIsPopupTrigger;
	}
}
