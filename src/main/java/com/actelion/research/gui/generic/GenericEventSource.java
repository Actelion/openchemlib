package com.actelion.research.gui.generic;

public interface GenericEventSource<T extends GenericEvent> {
	void addEventConsumer(GenericEventListener<T> consumer);
	void removeEventConsumer(GenericEventListener<T> consumer);
	void fireEvent(T event);
}
