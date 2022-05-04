package com.actelion.research.gui.generic;

public class GenericActionEvent extends GenericEvent {
	public static final int WHAT_OK = 0;
	public static final int WHAT_CANCEL = 1;
	public static final int WHAT_ITEM_SELECTED = 2;
	public static final int WHAT_STATE_TOGGLED = 3;

	private int mValue;
	private String mMessage;

	public GenericActionEvent(Object source, int what, int value) {
		super(what, source);
		mValue = value;
		}

	public GenericActionEvent(Object source, int what, String message) {
		super(what, source);
		mMessage = message;
		}

	public int getValue() {
		return mValue;
	}

	public String getMessage() {
		return mMessage;
	}
	}
