/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 *
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.gui;

import com.actelion.research.calc.ProgressController;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.actelion.research.gui.hidpi.HiDPIIconButton;
import info.clearthought.layout.TableLayout;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class JProgressPanel extends JPanel implements ActionListener,ProgressController {
	private static final long serialVersionUID = 0x20140404;

	private static final int SHOW_ERROR = 1;
	private static final int START_PROGRESS = 2;
	private static final int UPDATE_PROGRESS = 3;
	private static final int STOP_PROGRESS = 4;

	private JProgressBar mProgressBar = new JProgressBar();
	private JLabel mProgressLabel = new JLabel();
	private JButton mCancelButton;
	private volatile boolean mCancelAction;
	private volatile float mUpdateFactor;
	private volatile boolean mUpdateActionPending;
	private volatile int mUpdateStatus;

	public JProgressPanel(boolean showCancelButton) {
		Font oldFont = UIManager.getFont("Label.font");
		Font font = oldFont.deriveFont(Font.BOLD, 0.88f*oldFont.getSize());
		mProgressLabel.setForeground(LookAndFeelHelper.isDarkLookAndFeel() ? Color.RED.brighter() : Color.RED.darker());
		mProgressLabel.setFont(font);

		Dimension dim = new Dimension(HiDPIHelper.scale(80), HiDPIHelper.scale(10));
		mProgressBar.setVisible(false);
		mProgressBar.setPreferredSize(dim);
		mProgressBar.setMaximumSize(dim);
		mProgressBar.setMinimumSize(dim);
		mProgressBar.setSize(dim);

		double[][] size = { {TableLayout.PREFERRED, 4, TableLayout.PREFERRED, 4, TableLayout.FILL},
				{TableLayout.FILL, TableLayout.PREFERRED, TableLayout.FILL} };

		setLayout(new TableLayout(size));
		add(mProgressBar, "0,1");
		add(mProgressLabel, "4,0,4,2");

		if (showCancelButton) {
			mCancelButton = new HiDPIIconButton("closeButton.png", null, "close");
			mCancelButton.setVisible(false);
			mCancelButton.addActionListener(this);
			double[][] bs = { {TableLayout.PREFERRED}, {TableLayout.FILL, TableLayout.PREFERRED, TableLayout.FILL} };
			JPanel bp = new JPanel();
			bp.setLayout(new TableLayout(bs));
			bp.add(mCancelButton, "0,1");
			add(bp, "2,0,2,2");
		}
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("close")) {
			mCancelAction = true;
			return;
		}
	}

	public void cancel() {
		mCancelAction = true;
	}

	/**
	 * Resets the progress cancelled flag, which is set by pressing the cancel button.
	 */
	public void initializeThreadMustDie() {
		mCancelAction = false;
	}

	@Override
	public void startProgress(String text, int min, int max) {
		doActionThreadSafe(START_PROGRESS, text, min, max);
	}

	@Override
	public void updateProgress(int value) {
		doActionThreadSafe(UPDATE_PROGRESS, null, value, 0);
	}

	@Override
	public void updateProgress(int value, String message) {
		doActionThreadSafe(UPDATE_PROGRESS, message, value, 0);
	}

	@Override
	public void stopProgress() {
		doActionThreadSafe(STOP_PROGRESS, null, 0, 0);
	}

	@Override
	public void showErrorMessage(final String message) {
		doActionThreadSafe(SHOW_ERROR, message, 0, 0);
	}

	public void showMessage(final String message) {
		mProgressLabel.setText(message);
	}

	@Override
	public boolean threadMustDie() {
		return mCancelAction;
	}

	private void doActionThreadSafe(final int action, final String text, final int v1, final int v2) {
		if (action == START_PROGRESS) {
			mUpdateFactor = (v2 - v1 <= 10000) ? 1.0f : 10000.0f / (v2 - v1);
			mUpdateStatus = v1;
			}
		if (action == UPDATE_PROGRESS)
			mUpdateStatus = (v1 >= 0) ? v1 : mUpdateStatus - v1;

		if (SwingUtilities.isEventDispatchThread()) {
			doAction(action, text, v1, v2);
			}
		else {
			if (action == UPDATE_PROGRESS) {
				if (!mUpdateActionPending) {
					mUpdateActionPending = true;
					try {
						SwingUtilities.invokeLater(() -> {
							doAction(action, text, v1, v2);
							mUpdateActionPending = false;
							} );
						}
					catch (Exception e) {}
					}
				}
			else {
				try {
					SwingUtilities.invokeLater(() -> doAction(action, text, v1, v2) );
					}
				catch (Exception e) {}
				}
			}
		}

	private void doAction(final int action, final String text, final int v1, final int v2) {
		switch (action) {
			case SHOW_ERROR:
				Component c = this;
				while (!(c instanceof Frame))
					c = c.getParent();
				JOptionPane.showMessageDialog(c, text);
				break;
			case START_PROGRESS:
				mProgressBar.setVisible(true);
				mProgressBar.setMinimum(Math.round(mUpdateFactor*v1));
				mProgressBar.setMaximum(Math.round(mUpdateFactor*v2));
				mProgressBar.setValue(Math.round(mUpdateFactor*v1));
				if (mCancelButton != null)
					mCancelButton.setVisible(true);
				mProgressLabel.setText(text);
				break;
			case UPDATE_PROGRESS:
				mProgressBar.setValue(Math.round(mUpdateFactor*mUpdateStatus));
				if (text != null)
					mProgressLabel.setText(text);
				break;
			case STOP_PROGRESS:
				mProgressLabel.setText("");
				mProgressBar.setValue(mProgressBar.getMinimum());
				mProgressBar.setVisible(false);
				if (mCancelButton != null)
					mCancelButton.setVisible(false);
				break;
		}
	}
}
