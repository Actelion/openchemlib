package com.actelion.research.gui.hidpi;

import com.actelion.research.gui.VerticalFlowLayout;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

public class JBrowseButtons extends JPanel implements KeyListener {
	private boolean mShiftPressed;
	private final ActionListener mActionListener;

	public JBrowseButtons(boolean isVertical, ActionListener al) {
		this(isVertical, 4, 4, al);
		}

	public JBrowseButtons(boolean isVertical, int hgap, int vgap, ActionListener al) {
		mActionListener = al;
		JButton b1,b2,b3,b4;
		if (isVertical) {
			setLayout(new VerticalFlowLayout(VerticalFlowLayout.LEFT, VerticalFlowLayout.TOP, hgap, vgap, true));
			b1 = new HiDPIIconButton("toLast.png", null, "|<", 270);
			b2 = new HiDPIIconButton("toNext.png", null, "<", 270);
			b3 = new HiDPIIconButton("toNext.png", null, ">", 90);
			b4 = new HiDPIIconButton("toLast.png", null, ">|", 90);
			}
		else {
			setLayout(new FlowLayout(FlowLayout.CENTER, hgap, vgap));
			b1 = new HiDPIIconButton("toLast.png", null, "|<", 180);
			b2 = new HiDPIIconButton("toNext.png", null, "<", 180);
			b3 = new HiDPIIconButton("toNext.png", null, ">", 0);
			b4 = new HiDPIIconButton("toLast.png", null, ">|", 0);
			}
		b1.addActionListener(al);
		b2.addActionListener(al);
		b3.addActionListener(al);
		b4.addActionListener(al);
		add(b1);
		add(b2);
		add(b3);
		add(b4);

		addKeyListener(this);
		setFocusable(true);
		}

	@Override
	public void keyReleased(KeyEvent e) {}

	@Override
	public void keyTyped(KeyEvent e) {}

	@Override
	public void keyPressed(KeyEvent e) {
		int code = e.getKeyCode();
		if (code == KeyEvent.VK_SHIFT) {
			mShiftPressed = true;
			}
		else if (code == KeyEvent.VK_LEFT) {
			fireActionEvent(mShiftPressed ? "|<" : "<");
			}
		else if (code == KeyEvent.VK_RIGHT) {
			fireActionEvent(mShiftPressed ? ">|" : ">");
			}
		}

	private void fireActionEvent(String command) {
		mActionListener.actionPerformed(new ActionEvent(this, ActionEvent.ACTION_FIRST, command));
		}
	}
