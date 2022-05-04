package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.GenericLabel;
import javafx.scene.control.Label;

public class FXLabel extends FXComponent implements GenericLabel {
	public FXLabel(String text) {
		super(new Label(text));
	}

	@Override
	public void setText(String text) {
		((Label)getNode()).setText(text);
	}
}
