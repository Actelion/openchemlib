package com.actelion.research.gui.generic;

public interface GenericComboBox extends GenericComponent {
	void addItem(String item);
	void removeAllItems();
	int getSelectedIndex();
	String getSelectedItem();
	void setSelectedIndex(int index);
	void setSelectedItem(String item);
	void setEditable(boolean b);
}
