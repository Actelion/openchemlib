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

import com.actelion.research.util.Prefs;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;

public class JLoginDialog extends JDialog implements WindowListener {
    private static final long serialVersionUID = 0x20100518;

    public static final String cLoginCancel = "loginCancel";
	public static final String cLoginOK = "loginOK";
	public static final String PREFERENCES_KEY_USERID = "databaseUserID";

	private JTextField		mTextFieldUserID;
	private JPasswordField	mTextFieldPassword;

	public JLoginDialog(Dialog owner, ActionListener listener) {
		super(owner, "Database Login", true);
		initialize(listener);
		setLocationRelativeTo(owner);
		}

	public JLoginDialog(Frame owner, ActionListener listener) {
		super(owner, "Database Login", true);
		initialize(listener);
		setLocationRelativeTo(owner);
		}

	private void initialize(ActionListener listener) {
		JPanel p1 = new JPanel();
		p1.setBorder(BorderFactory.createEmptyBorder(8, 8, 8, 8));
		p1.setLayout(new GridLayout(2, 2, 4, 4));
		p1.add(new JLabel("User-ID:", SwingConstants.TRAILING));
		mTextFieldUserID = new JTextField(12);
		p1.add(mTextFieldUserID);
		p1.add(new JLabel("Password:", SwingConstants.TRAILING));
		mTextFieldPassword = new JPasswordField(12);
		p1.add(mTextFieldPassword);

		JPanel p2 = new JPanel();
		p2.setBorder(BorderFactory.createEmptyBorder(12, 8, 8, 8));
		p2.setLayout(new BorderLayout());
		JPanel bp = new JPanel();
		bp.setLayout(new GridLayout(1, 6, 8, 0));
		JButton bcancel = new JButton("Cancel");
		bcancel.setActionCommand(cLoginCancel);
		bcancel.addActionListener(listener);
		bp.add(bcancel);
		JButton bok = new JButton("OK");
		bok.setActionCommand(cLoginOK);
		bok.addActionListener(listener);
		bp.add(bok);
		p2.add(bp, BorderLayout.EAST);

		getContentPane().add(p1, BorderLayout.CENTER);
		getContentPane().add(p2, BorderLayout.SOUTH);
		getRootPane().setDefaultButton(bok);

		String userid = Prefs.getString(PREFERENCES_KEY_USERID, "");
		if (userid != null) {
			mTextFieldUserID.setText(userid);
			addWindowListener(this);
				// a little strange just to set initial focus, but it works
			}

		pack();
		}

	public String getUserID() {
		String userid = mTextFieldUserID.getText();
		if (userid.length() != 0)
			Prefs.setString(PREFERENCES_KEY_USERID, userid);

		return userid;
		}

	public String getPassword() {
		return new String(mTextFieldPassword.getPassword());
		}

    public void windowOpened(WindowEvent e) {
		SwingUtilities.invokeLater(() -> mTextFieldPassword.requestFocus());
		}
    public void windowClosing(WindowEvent e) {}
    public void windowClosed(WindowEvent e) {}
    public void windowIconified(WindowEvent e) {}
    public void windowDeiconified(WindowEvent e) {}
    public void windowActivated(WindowEvent e) {}
    public void windowDeactivated(WindowEvent e) {}
	}