/*
 * @(#)MessageBar.java
 *
 * Copyright 1997-2001 Actelion Ltd., Inc. All Rights Reserved.
 * 
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 * 
 * @author Thomas Sander
 */

package com.actelion.research.gui;

import java.awt.*;

import javax.swing.*;

public class JMessageBar extends JPanel {
	private static JMessageBar sMessageBar = null;
	private static Color sBackground = new Color(99,99,99);
	private static Color sTextColor = Color.GREEN;
	private static Color sBusyColor = Color.ORANGE;
	private static Color sErrorColor = Color.RED;
	private JLabel mLabel;
	private String mText;

	public static JMessageBar getBar() {
		if (sMessageBar == null)
			sMessageBar = new JMessageBar();
		return sMessageBar;
		}

	private JMessageBar() {
		setLayout(new BorderLayout());

		mLabel = new JLabel("", SwingConstants.CENTER);
		mLabel.setFont(mLabel.getFont().deriveFont(Font.PLAIN,12));
		mLabel.setForeground(sTextColor);
		mLabel.setBackground(sBackground);
		mLabel.setOpaque(true);
		add(mLabel,BorderLayout.CENTER);
		}

	public static void setTextColor(Color fg) {
		sTextColor = fg;
		}

	public static void setBusyColor(Color fg) {
		sTextColor = fg;
		}

	public static void setErrorColor(Color fg) {
		sTextColor = fg;
		}

	public static void setBackgroundColor(Color bg) {
		sBackground = bg;
		}

	public void setText(String theText) {
		mLabel.setForeground(sTextColor);
		mLabel.setText(theText);
		mText = theText;
		}

	public void showText() {
		mLabel.setForeground(sTextColor);
		mLabel.setText(mText);
		}

	public void showBusyMessage(String msg) {
		mLabel.setForeground(sBusyColor);
		mLabel.setText(msg);
		}

	public void showErrorMessage(String msg) {
		mLabel.setForeground(sErrorColor);
		mLabel.setText(msg);
		}
	}
