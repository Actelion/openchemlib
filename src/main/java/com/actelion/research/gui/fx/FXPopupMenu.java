package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.GenericActionEvent;
import com.actelion.research.gui.generic.GenericEventListener;
import com.actelion.research.gui.generic.GenericPopupMenu;
import javafx.geometry.Point2D;
import javafx.scene.Node;
import javafx.scene.control.*;


public class FXPopupMenu extends FXComponent implements GenericPopupMenu {
	private ContextMenu mPopupMenu;
	private Menu mSubMenu;
	private Node mOwner;

	public FXPopupMenu(Node owner, GenericEventListener<GenericActionEvent> consumer) {
		super(null);
		mPopupMenu = new ContextMenu();
		mOwner = owner;
		addEventConsumer(consumer);
	}

	@Override
	public void addItem(String text, String command, boolean enabled) {
		MenuItem item = new MenuItem(text);
		item.disableProperty().set(!enabled);
		String _command = (command == null) ? text : command;
		item.setOnAction(event -> fireEvent(new GenericActionEvent(this, GenericActionEvent.WHAT_ITEM_SELECTED, _command)));

		if (mSubMenu != null)
			mSubMenu.getItems().add(item);
		else
			mPopupMenu.getItems().add(item);
	}

	@Override
	public void addRadioButtonItem(String text, String command, int colorRGB, boolean isSelected) {
		RadioMenuItem item = new RadioMenuItem(text);
		item.selectedProperty().set(isSelected);
		String _command = (command == null) ? text : command;
		item.setOnAction(event -> fireEvent(new GenericActionEvent(this, GenericActionEvent.WHAT_STATE_TOGGLED, _command)));

		if (colorRGB != 0)
			item.setStyle("-fx-background-color: #"+Integer.toHexString(colorRGB & 0x00FFFFFF)+"; ");
		if (mSubMenu != null)
			mSubMenu.getItems().add(item);
		else
			mPopupMenu.getItems().add(item);
	}

	@Override
	public void startSubMenu(String text) {
		mSubMenu = new Menu(text);
	}

	@Override
	public void endSubMenu() {
		mPopupMenu.getItems().add(mSubMenu);
		mSubMenu = null;
	}

	@Override
	public void addSeparator() {
		if (mSubMenu != null)
			mSubMenu.getItems().add(new SeparatorMenuItem());
		else
			mPopupMenu.getItems().add(new SeparatorMenuItem());
	}

	@Override
	public void show(int x, int y) {
		Point2D p = mOwner.localToScreen(x, y);
		mPopupMenu.show(mOwner.getScene().getWindow(), p.getX(), p.getY());
	}
}
