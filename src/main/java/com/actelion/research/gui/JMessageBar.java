/*
 * Copyright (c) 1997 - 2016
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the the copyright holder nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Thomas Sander
 */

package com.actelion.research.gui;

import com.actelion.research.gui.hidpi.HiDPIHelper;

import javax.swing.*;
import java.awt.*;

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
		mLabel.setFont(mLabel.getFont().deriveFont(Font.PLAIN, HiDPIHelper.scale(12)));
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
