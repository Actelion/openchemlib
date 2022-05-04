package com.actelion.research.gui.generic;

import java.util.ArrayList;

public abstract class GenericEventHandler<T extends GenericEvent> {
	private ArrayList<GenericEventListener<T>> mListeners;
	private Object mSource;

	public GenericEventHandler(Object source) {
		mSource = source;
		mListeners = new ArrayList<>();
	}

	public Object getSource() {
		return mSource;
	}

	public void addListener(GenericEventListener<T> e) {
		mListeners.add(e);
	}

	public void removeListener(GenericEventListener<T> e) {
		mListeners.remove(e);
	}

	public void fireEvent(T e) {
		for (GenericEventListener<T> gel:mListeners)
			gel.eventHappened(e);
	}
}
