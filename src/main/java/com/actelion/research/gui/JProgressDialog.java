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

import com.actelion.research.calc.ProgressController;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import info.clearthought.layout.TableLayout;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class JProgressDialog implements ActionListener,ProgressController {
	private static final int sActionDispose = 2;
	private static final int sActionStart = 4;
	private static final int sActionUpdate = 8;
	private static final int sActionMessage = 16;
	private static final int sActionStop = 32;

	private volatile JDialog mDialog;
	private Frame			mOwner;
	private JProgressBar	mProgressBar;
	private JLabel			mProgressLabel,mLabelRemainingTime;
	private Frame			mNewFrontFrame;
	private volatile String	mBusyText;
	private volatile String	mProgressMessage;
	private volatile boolean mProcessCancelled;
	private volatile int	mProgressMin,mProgressMax,mProgressValue,mAction;
	private volatile long   mProgressStart,mLastUpdate;
	private volatile float  mProgressFactor;

	/**
	 * Creates a facade to a progress dialog, which is not instantiated and shown immediately.
	 * The real dialog is instantiated and shown with the first call of startProgress().
	 * If startProgress() is never called, the dialog is neither created nor shown and none of the
	 * other actions have any effect.
	 * @param owner may be null, if later setOwner() is called before the first call of startProgress()
	 */
	public JProgressDialog(Frame owner) {
		mOwner = owner;
		}

	public void setOwner(Frame owner) {
		mOwner = owner;
	}

	/**
	 * The first call of startProgress() instatiates the
	 * May be called safely from any thread
	 * @param text
	 * @param min
	 * @param max
	 */
	@Override
	public void startProgress(String text, int min, int max) {
		if (mDialog == null && mOwner != null) {
			if (SwingUtilities.isEventDispatchThread())
				initialize();
			else
				SwingUtilities.invokeLater(() -> initialize() );
		}

		mProgressStart = System.currentTimeMillis();
		mBusyText = text;
		mProgressMin = min;
		mProgressValue = min;
		mProgressMax = max;
		mProgressFactor = (max - min <= 10000) ? 1.0f : 10000.0f / (max - min);
		mAction |= sActionStart;
		update();
		}

	/**
	 * Update progress status in an absolute or relative way.
	 * @param value if negative, its abs value is added to current progress.
	 */
	@Override
	public void updateProgress(int value) {
			// may be called safely from any thread
		int newValue = (value >= 0) ? value : mProgressValue-value;
		if (mProgressValue < newValue) {
			mProgressValue = newValue;
			mAction |= sActionUpdate;
			update();
			}
		}

	/**
	 * Update progress status in an absolute or relative way.
	 * This method also updates the message of the progress dialog
	 * @param value if negative, its abs value is added to current progress.
	 * @param message
	 */
	@Override
	public void updateProgress(int value, String message) {
		// may be called safely from any thread
		int newValue = (value >= 0) ? value : mProgressValue-value;
		if (mProgressValue < newValue) {
			mProgressValue = newValue;
			mProgressMessage = message;
			mAction |= sActionUpdate | sActionMessage;
			update();
			}
		}

	/**
	 * Sets progress back to zero, hides the progress bar, but keeps the dialog open.
	 * May be called safely from any thread
	 */
	@Override
	public void stopProgress() {
		mAction |= sActionStop;
		update();
		}

	/**
	 * Closes and disposes the dialog.
	 * If you need to move a newly created window to the front, then use close(Frame newFrontFrame).
	 * May be called safely from any thread
	 */
	public void close() {
		close(null);
		}

	/**
	 * Disposes of the progress dialog and optionally moves the specified
	 * frame to the front. This may be called safely from any thread.
	 * @param newFrontFrame null or frame to be moved to the front
	 */
	public void close(Frame newFrontFrame) {
		mAction |= sActionDispose;
		mNewFrontFrame = newFrontFrame;
		update();
		}

	/**
	 * Shows the given message in an interactive dialog.
	 * This may be called safely from any thread.
	 * @param message
	 */
	@Override
	public void showErrorMessage(final String message) {
		if (SwingUtilities.isEventDispatchThread()) {
			showMessage(message);
			}
		else {
			try {
				SwingUtilities.invokeLater(() -> showMessage(message) );
				}
			catch (Exception e) {}
			}
		}

	private void showMessage(String message) {
		JOptionPane.showMessageDialog(mDialog, message);
		}

	@Override
	public boolean threadMustDie() {
		return mProcessCancelled;
		}

	@Override
	public void actionPerformed(ActionEvent e) {
		mProcessCancelled = true;
		mAction |= sActionStop;
		update();
		}

	private void initialize() {
		mDialog = new JDialog(mOwner, true);

		double[][] size = { {8, HiDPIHelper.scale(200), 8, TableLayout.PREFERRED, 8},
							{8, TableLayout.PREFERRED, 8, TableLayout.PREFERRED, 16, TableLayout.PREFERRED, 8} };
		mDialog.getContentPane().setLayout(new TableLayout(size));

		mProgressLabel = new JLabel(" ");
		mDialog.getContentPane().add(mProgressLabel, "1,1,3,1");

		mProgressBar = new JProgressBar();
		mProgressBar.setVisible(false);
		mDialog.getContentPane().add(mProgressBar, "1,3,3,3");

		mLabelRemainingTime = new JLabel();
		mDialog.getContentPane().add(mLabelRemainingTime, "1,5");

		JButton cancelButton = new JButton("Cancel");
		cancelButton.addActionListener(this);
		mDialog.getContentPane().add(cancelButton, "3,5");
		mDialog.pack();
		mDialog.setLocationRelativeTo(mOwner);
		mDialog.setVisible(true);
		}

	private void update() {
		if (mDialog != null) {
			long millis = System.currentTimeMillis();
			if ((mAction & ~sActionUpdate) != 0 || millis - mLastUpdate > 100) {
				if (SwingUtilities.isEventDispatchThread()) {
					mLastUpdate = millis;
					doUpdate();
					}
				else {
					mLastUpdate = millis;
					SwingUtilities.invokeLater(() -> doUpdate());
					}
				}
			}
		}

	private void doUpdate() {
		if ((mAction & sActionDispose) != 0) {
			mAction &= ~sActionDispose;
			mDialog.setVisible(false);
			mDialog.dispose();
			if (mNewFrontFrame != null) {
				SwingUtilities.invokeLater(() -> mNewFrontFrame.toFront() );
			}
			return;
			}
		if ((mAction & sActionStart) != 0) {
			mAction &= ~sActionStart;
			if (mProgressMin == mProgressMax) {
				mProgressBar.setIndeterminate(true);
				}
			else {
				mProgressBar.setIndeterminate(false);
				mProgressBar.setMinimum(Math.round(mProgressFactor * mProgressMin));
				mProgressBar.setMaximum(Math.round(mProgressFactor * mProgressMax));
				mProgressBar.setValue(Math.round(mProgressFactor * mProgressMin));
				}
			mProgressBar.setVisible(true);
			mProgressLabel.setText(mBusyText);
			}
		if ((mAction & sActionUpdate) != 0) {
			mAction &= ~sActionUpdate;
			mProgressBar.setValue(Math.round(mProgressFactor * mProgressValue));
			if (mProgressValue > mProgressMin) {
				long milliesUsed = System.currentTimeMillis() - mProgressStart;
				long milliesToGo = milliesUsed * (mProgressMax - mProgressValue) / (mProgressValue - mProgressMin);
				long millis = (mProgressMax == mProgressMin) ? milliesUsed : milliesToGo;
				String timeToShow = (millis >= 7200000) ? ""+(millis/3600000)+" hours"
								  : (millis >=  120000) ? ""+(millis/60000)+" minutes"
								  : 					  ""+(millis/1000)+" seconds";
				String message = (mProgressMax == mProgressMin) ? " waiting" : " remaining";
				mLabelRemainingTime.setText(timeToShow+message);
				}
			if ((mAction & sActionMessage) != 0) {
				mAction &= ~sActionMessage;
				mProgressLabel.setText(mProgressMessage);
				}
			}
		if ((mAction & sActionStop) != 0) {
			mAction &= ~sActionStop;
			mProgressLabel.setText(mProcessCancelled ? "Cancelled; cleaning up..." : "");
			if (mProcessCancelled)
				mProgressBar.setIndeterminate(true);
			else {
				mProgressBar.setValue(mProgressBar.getMinimum());
				mProgressBar.setVisible(false);
				}
			}
		}
	}
