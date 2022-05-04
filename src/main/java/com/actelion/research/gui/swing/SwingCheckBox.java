package com.actelion.research.gui.swing;

import com.actelion.research.gui.generic.GenericActionEvent;
import com.actelion.research.gui.generic.GenericCheckBox;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class SwingCheckBox extends SwingComponent implements GenericCheckBox,ActionListener {
	private JCheckBox mCheckBox;

	public SwingCheckBox(String text) {
		super(new JCheckBox(text));
		mCheckBox = (JCheckBox)getComponent();
		mCheckBox.addActionListener(this);
	}

	@Override
	public boolean isSelected() {
		return mCheckBox.isSelected();
	}

	@Override
	public void setSelected(boolean b) {
		mCheckBox.setSelected(b);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		fireEvent(new GenericActionEvent(this, GenericActionEvent.WHAT_STATE_TOGGLED, ((JCheckBox)e.getSource()).isSelected() ? 1 : 0));
	}
}
