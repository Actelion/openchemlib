package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.GenericActionEvent;
import com.actelion.research.gui.generic.GenericCheckBox;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.control.CheckBox;


public class FXCheckBox extends FXComponent implements GenericCheckBox,EventHandler<ActionEvent> {
	private CheckBox mCheckBox;

	public FXCheckBox(String text) {
		super(new CheckBox(text));
		mCheckBox = (CheckBox)getNode();
		mCheckBox.addEventHandler(javafx.event.ActionEvent.ACTION, this);
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
	public void handle(ActionEvent event) {
		fireEvent(new GenericActionEvent(this, GenericActionEvent.WHAT_STATE_TOGGLED, isSelected() ? 1 : 0));
	}
}
