package com.actelion.research.gui.swing;

import com.actelion.research.gui.editor.DialogEventConsumer;
import com.actelion.research.gui.generic.GenericComponent;

import javax.swing.*;

public class SwingComponent implements GenericComponent {
	private JComponent mComponent;
	private DialogEventConsumer mConsumer;

	public SwingComponent(JComponent c) {
		mComponent = c;
	}

	@Override
	public void setEventConsumer(DialogEventConsumer consumer) {
		mConsumer = consumer;
	}

	@Override
	public void setEnabled(boolean b) {
		mComponent.setEnabled(b);
	}

	public JComponent getComponent() {
		return mComponent;
	}

	@Override
	public DialogEventConsumer getEventConsumer() {
		return mConsumer;
	}
}
