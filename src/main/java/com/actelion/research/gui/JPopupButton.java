package com.actelion.research.gui;

import com.actelion.research.gui.hidpi.HiDPIIconButton;

import javax.swing.*;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

public class JPopupButton extends HiDPIIconButton {
	private JPopupMenu mPopup;
	private ActionListener mListener;

	public JPopupButton(ActionListener al) {
		super("popup14.png", null, "showPopup");
		mListener = al;
		mPopup = new JPopupMenu();

		addMouseListener(new MouseAdapter() {
			@Override
			public void mousePressed(MouseEvent e) {
				showPopup(e);
			}
			@Override
			public void mouseReleased(MouseEvent e) {
				mPopup.setVisible(false);
			}
		} );
	}

	public void addItem(String text, String command) {
		JMenuItem item = new JMenuItem(text);
		item.setActionCommand(command);
		item.addActionListener(mListener);
		mPopup.add(item);
	}

	private void showPopup(MouseEvent e) {
		mPopup.show(getParent(), getBounds().x, getBounds().y + getBounds().height);
		}

	private void hidePopup(MouseEvent e) {
		mPopup.setVisible(false);
	}
}
