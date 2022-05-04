package com.actelion.research.gui.swing;

import com.actelion.research.gui.generic.GenericActionEvent;
import com.actelion.research.gui.generic.GenericComponent;
import com.actelion.research.gui.generic.GenericEventListener;

import javax.swing.*;
import java.util.ArrayList;

public class SwingComponent implements GenericComponent {
	private JComponent mComponent;
	private ArrayList<GenericEventListener<GenericActionEvent>> mConsumerList;

	public SwingComponent(JComponent c) {
		mComponent = c;
		mConsumerList = new ArrayList<>();
	}

	@Override
	public void addEventConsumer(GenericEventListener<GenericActionEvent> consumer) {
		mConsumerList.add(consumer);
	}

	@Override
	public void removeEventConsumer(GenericEventListener<GenericActionEvent> consumer) {
		mConsumerList.remove(consumer);
	}

	@Override
	public void setEnabled(boolean b) {
		mComponent.setEnabled(b);
	}

	public JComponent getComponent() {
		return mComponent;
	}

	@Override
	public void fireEvent(GenericActionEvent event) {
		for (GenericEventListener<GenericActionEvent> consumer:mConsumerList)
			consumer.eventHappened(event);
	}
}
