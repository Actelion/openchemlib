package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.GenericActionEvent;
import com.actelion.research.gui.generic.GenericComboBox;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.control.ComboBox;

public class FXComboBox extends FXComponent implements GenericComboBox,EventHandler<ActionEvent> {
	private ComboBox mComboBox;

	public FXComboBox() {
		super(new ComboBox<String>());
		mComboBox = (ComboBox)getNode();
		mComboBox.addEventHandler(javafx.event.ActionEvent.ACTION, this);
	}

	@Override
	public void handle(ActionEvent event) {
		fireEvent(new GenericActionEvent(this, GenericActionEvent.WHAT_ITEM_SELECTED, mComboBox.getSelectionModel().getSelectedIndex()));
	}

	@Override
	public void removeAllItems() {
		mComboBox.getItems().removeAll();
	}

	@Override
	public void addItem(String item) {
		mComboBox.getItems().add(item);
	}

	@Override
	public int getSelectedIndex() {
		return mComboBox.getSelectionModel().getSelectedIndex();
	}

	@Override
	public String getSelectedItem() {
		return (String)mComboBox.getSelectionModel().getSelectedItem();
	}

	@Override
	public void setSelectedIndex(int index) {
		mComboBox.getSelectionModel().select(index);
	}

	@Override
	public void setSelectedItem(String item) {
		mComboBox.getSelectionModel().select(item);
	}

	@Override
	public void setEditable(boolean b) {
		mComboBox.setEditable(b);
	}
}
