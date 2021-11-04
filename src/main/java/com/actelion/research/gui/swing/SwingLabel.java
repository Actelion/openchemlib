package com.actelion.research.gui.swing;

import com.actelion.research.gui.generic.GenericLabel;

import javax.swing.*;

public class SwingLabel extends SwingComponent implements GenericLabel {
	public SwingLabel(String text) {
		super(new JLabel(text));
	}

	@Override
	public void setText(String text) {
		((JLabel)getComponent()).setText(text);
	}
}
