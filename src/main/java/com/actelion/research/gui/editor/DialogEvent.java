package com.actelion.research.gui.editor;

import com.actelion.research.gui.generic.GenericComponent;

public class DialogEvent {
	public static final int WHAT_OK = 0;
	public static final int WHAT_CANCEL = 1;
	public static final int WHAT_ITEM_SELECTED = 2;
	public static final int WHAT_STATE_TOGGLED = 3;

	private GenericComponent mSource;
	private int mWhat;
	private int mValue;
	private String mMessage;

	public DialogEvent(GenericComponent source, int what, int value) {
		mSource = source;
		mWhat = what;
		mValue = value;
		}

	public DialogEvent(GenericComponent source, int what, String message) {
		mSource = source;
		mWhat = what;
		mMessage = message;
		}

	public GenericComponent getSource() {
		return mSource;
	}

	public int getWhat() {
		return mWhat;
		}

	public int getValue() {
		return mValue;
	}

	public String getMessage() {
		return mMessage;
	}
	}
