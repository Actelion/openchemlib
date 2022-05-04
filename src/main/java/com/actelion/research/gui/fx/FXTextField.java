package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.GenericTextField;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.TextInputControl;

public class FXTextField extends FXComponent implements GenericTextField {
	public FXTextField(int width, int height) {
		super(height == 1 ? new TextField() : new TextArea());
		((TextInputControl)getNode()).setPrefWidth(width);
		((TextInputControl)getNode()).setPrefHeight(height);
	}

	public String getText() {
		return ((TextInputControl)getNode()).getText();
	}

	public void setText(String text) {
		((TextInputControl)getNode()).setText(text);
	}
}
