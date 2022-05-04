package com.actelion.research.gui.generic;

public class GenericMouseEvent extends GenericEvent {
	public static final int MOUSE_PRESSED = 1;
	public static final int MOUSE_RELEASED = 2;
	public static final int MOUSE_CLICKED = 3;
	public static final int MOUSE_ENTERED = 4;
	public static final int MOUSE_EXITED = 5;
	public static final int MOUSE_MOVED = 6;
	public static final int MOUSE_DRAGGED = 7;
	private static final String[] WHAT_MESSAGE = { "none", "pressed", "released", "clicked", "entered", "exited", "moved", "dragged" };

	private int mButton,mClickCount,mX,mY;
	private boolean mShiftDown,mControlDown,mAltDown,mIsPopupTrigger;

	public GenericMouseEvent(int what, int button, int clickCount, int x, int y, boolean shiftDown, boolean controlDown, boolean altDown, boolean isPopupTrigger, Object source) {
		super(what, source);
		mButton = button;
		mClickCount = clickCount;
		mX = x;
		mY = y;
		mShiftDown = shiftDown;
		mControlDown = controlDown;
		mAltDown = altDown;
		mIsPopupTrigger = isPopupTrigger;
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

	@Override
	public String toString() {
		return WHAT_MESSAGE[getWhat()]+" x:"+mX+" y:"+mY+" button:"+mButton+" clicks:"+mClickCount
				+(mShiftDown?" shift":"")+(mControlDown?" ctrl":"")+(mAltDown?" alt":"")+" isPopup:"+(mIsPopupTrigger?"y":"n");
	}
}
