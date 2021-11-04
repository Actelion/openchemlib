package com.actelion.research.gui.swing;

import com.actelion.research.gui.editor.DialogEvent;
import com.actelion.research.gui.editor.DialogEventConsumer;
import com.actelion.research.gui.generic.GenericPopupMenu;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class SwingPopupMenu extends SwingComponent implements ActionListener,GenericPopupMenu {
	private JPopupMenu mPopupMenu;
	private JMenu mSubMenu;
	private JComponent mOwner;

	public SwingPopupMenu(JComponent owner, DialogEventConsumer consumer) {
		super(new JPopupMenu());
		mPopupMenu = (JPopupMenu)getComponent();
		mOwner = owner;
		setEventConsumer(consumer);
	}

	@Override
	public void addItem(String text, String command, boolean enabled) {
		JMenuItem item = new JMenuItem(text);
		item.setEnabled(enabled);
		if (command != null)
			item.setActionCommand(command);
		item.addActionListener(this);
		if (mSubMenu != null)
			mSubMenu.add(item);
		else
			mPopupMenu.add(item);
		}

	@Override
	public void addRadioButtonItem(String text, String command, Color color, boolean isSelected) {
		JRadioButtonMenuItem item = new JRadioButtonMenuItem(text, isSelected);
		if (command != null)
			item.setActionCommand(command);
		if (color != null)
			item.setBackground(color);
		item.addActionListener(this);
		if (mSubMenu != null)
			mSubMenu.add(item);
		else
			mPopupMenu.add(item);
		}

	@Override
	public void startSubMenu(String text) {
		mSubMenu = new JMenu(text);
		}

	@Override
	public void endSubMenu() {
		mPopupMenu.add(mSubMenu);
		mSubMenu = null;
	}

	@Override
	public void addSeparator() {
		if (mSubMenu != null)
			mSubMenu.addSeparator();
		else
			mPopupMenu.addSeparator();
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		getEventConsumer().dialogEventHappened(new DialogEvent(this, DialogEvent.WHAT_ITEM_SELECTED, e.getActionCommand()));
	}

	@Override
	public void show(int x, int y) {
		mPopupMenu.show(mOwner, x, y);
	}
}
