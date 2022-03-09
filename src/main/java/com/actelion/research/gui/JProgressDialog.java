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

public class JProgressDialog extends JDialog implements ActionListener,ProgressController,Runnable {
	private static final long serialVersionUID = 0x20070301;

	private static final int sActionDispose = 2;
	private static final int sActionStart = 4;
	private static final int sActionUpdate = 8;
	private static final int sActionMessage = 16;
	private static final int sActionStop = 32;

	private JProgressBar	mProgressBar;
	private JLabel			mProgressLabel,mLabelRemainingTime;
	private Frame			mNewFrontFrame;
	private String			mBusyText;
	private volatile String	mProgressMessage;
	private volatile boolean mProcessCancelled;
	private volatile int	mProgressMin,mProgressMax,mProgressValue,mAction;
	private volatile long   mProgressStart,mLastUpdate;
	private volatile float  mProgressFactor;

	/**
	 * Creates a JProgressDialog and schedule it with invokeLater()
	 * to be set visible without blocking the calling thread.
	 * @param owner
	 */
	public JProgressDialog(Frame owner) {
		this(owner, true);
		}

	/**
	 * Creates a JProgressDialog. If invokeSetVisible is true,
	 * then the dialog is scheduled to be set visible with invokeLater()
	 * without blocking the current thread. Otherwise the caller needs
	 * to call setVisible() manually.
	 * @param owner
	 * @param invokeSetVisible
	 */
	public JProgressDialog(Frame owner, boolean invokeSetVisible) {
			// initialized and shows the modal dialog without blocking the thread
		super(owner, true);
		initialize();
		setLocationRelativeTo(owner);
		if (invokeSetVisible)
			SwingUtilities.invokeLater(() -> setVisible(true) );
		}

	public void startProgress(String text, int min, int max) {
			// may be called safely from any thread
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

	public void showErrorMessage(final String message) {
		if (SwingUtilities.isEventDispatchThread()) {
			showMessage(message);
			}
		else {
			try {
				SwingUtilities.invokeAndWait(() -> showMessage(message) );
				}
			catch (Exception e) {}
			}
		}		

	private void showMessage(String message) {
		JOptionPane.showMessageDialog(this, message);
		}

	public boolean threadMustDie() {
		return mProcessCancelled;
		}

	public void actionPerformed(ActionEvent e) {
		mProcessCancelled = true;
		mAction |= sActionStop;
		update();
		}

	private void initialize() {
		double[][] size = { {8, HiDPIHelper.scale(200), 8, TableLayout.PREFERRED, 8},
							{8, TableLayout.PREFERRED, 8, TableLayout.PREFERRED, 16, TableLayout.PREFERRED, 8} };
		getContentPane().setLayout(new TableLayout(size));

		mProgressLabel = new JLabel(" ");
		getContentPane().add(mProgressLabel, "1,1,3,1");

		mProgressBar = new JProgressBar();
		mProgressBar.setVisible(false);
		getContentPane().add(mProgressBar, "1,3,3,3");

		mLabelRemainingTime = new JLabel();
		getContentPane().add(mLabelRemainingTime, "1,5");

		JButton cancelButton = new JButton("Cancel");
		cancelButton.addActionListener(this);
		getContentPane().add(cancelButton, "3,5");
		pack();
		}

	private void update() {
		long millis = System.currentTimeMillis();
		if ((mAction & ~sActionUpdate) != 0 || millis - mLastUpdate > 100) {
			if (SwingUtilities.isEventDispatchThread()) {
				mLastUpdate = millis;
				run();
				}
			else {
//				try {
					mLastUpdate = millis;
					SwingUtilities.invokeLater(this);
   //				 }
 //			   catch (Exception e) {}
				}
			}		
		}

	public void run() {
		if ((mAction & sActionDispose) != 0) {
			mAction &= ~sActionDispose;
			setVisible(false);
			dispose();
			if (mNewFrontFrame != null)
				mNewFrontFrame.toFront();
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
