package com.actelion.research.gui.generic;

import info.clearthought.layout.TableLayout;

public interface GenericDialog extends GenericComponent {
	int PREFERRED = (int)TableLayout.PREFERRED;
	int FILL = (int)TableLayout.FILL;

	void setLayout(int[] hLayout, int[] vLayout);
	void add(GenericComponent c, int x, int y);
	void add(GenericComponent c, int x1, int y1, int x2, int y2);
	GenericCheckBox createCheckBox(String text);
	GenericComboBox createComboBox();
	GenericLabel createLabel(String text);
	GenericTextField createTextField(int width, int height);
	void showDialog();
	void disposeDialog();
	void showMessage(String message);
}
