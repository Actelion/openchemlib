package com.actelion.research.gui.generic;

public interface GenericPopupMenu extends GenericComponent {
	void addItem(String text, String command, boolean enabled);
	void addRadioButtonItem(String text, String command, int rgb, boolean isSelected);
	void startSubMenu(String text);
	void endSubMenu();
	void addSeparator();
	void show(int x, int y);
}
