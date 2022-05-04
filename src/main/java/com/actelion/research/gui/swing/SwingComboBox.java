package com.actelion.research.gui.swing;

import com.actelion.research.gui.generic.GenericActionEvent;
import com.actelion.research.gui.generic.GenericComboBox;

import javax.swing.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

public class SwingComboBox extends SwingComponent implements GenericComboBox,ItemListener {
	private JComboBox mComboBox;

	public SwingComboBox() {
		super(new JComboBox<String>());
		mComboBox = (JComboBox)getComponent();
		mComboBox.addItemListener(this);
	}

	@Override
	public void itemStateChanged(ItemEvent e) {
		if (e.getStateChange() == ItemEvent.SELECTED)
			fireEvent(new GenericActionEvent(this, GenericActionEvent.WHAT_ITEM_SELECTED, mComboBox.getSelectedIndex()));
	}

	@Override
	public void removeAllItems() {
		mComboBox.removeAllItems();
		}

	@Override
	public void addItem(String item) {
		mComboBox.addItem(item);
	}

	@Override
	public int getSelectedIndex() {
		return mComboBox.getSelectedIndex();
	}

	@Override
	public String getSelectedItem() {
		return (String)mComboBox.getSelectedItem();
	}

	@Override
	public void setSelectedIndex(int index) {
		mComboBox.setSelectedIndex(index);
	}

	@Override
	public void setSelectedItem(String item) {
		mComboBox.setSelectedItem(item);
	}

	@Override
	public void setEditable(boolean b) {
		mComboBox.setEditable(b);
	}
}
