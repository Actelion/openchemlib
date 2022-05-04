package com.actelion.research.gui.generic;

import java.util.EventObject;

public class GenericEvent extends EventObject {
	private int mWhat;

	public GenericEvent(int what, Object source) {
		super(source);
		mWhat = what;
	}

	public int getWhat() {
		return mWhat;
	}
}
