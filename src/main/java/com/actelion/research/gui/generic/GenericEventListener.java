package com.actelion.research.gui.generic;

public interface GenericEventListener<T extends GenericEvent> {
	void eventHappened(T e);
}
