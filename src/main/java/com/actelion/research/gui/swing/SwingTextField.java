package com.actelion.research.gui.swing;

import com.actelion.research.gui.generic.GenericTextField;

import javax.swing.*;
import javax.swing.text.JTextComponent;

public class SwingTextField extends SwingComponent implements GenericTextField {
	public SwingTextField(int width, int height) {
		super(height == 1 ? new JTextField(width/2) : new JTextArea(height, width));
	}

	public String getText() {
		return ((JTextComponent)getComponent()).getText();
	}

	public void setText(String text) {
		((JTextComponent)getComponent()).setText(text);
	}
}
