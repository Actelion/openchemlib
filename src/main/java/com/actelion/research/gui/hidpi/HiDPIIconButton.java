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
*/

package com.actelion.research.gui.hidpi;

import javax.swing.*;
import java.awt.*;

/**
 * Created by sandert on 04/12/15.
 */
public class HiDPIIconButton extends JButton {
	private final String mImageName;
	private final int mRotation;
	private Icon[] mAnimationIcon;
	private Animator mAnimator;

	/**
	 * Creates a button that, if image2 is given, toggles between two states indicated
	 * by two different button images. The button optimizes its size to for QuaQua
	 * and Substance look&feels and uses adequate higher resolution images on HiDPI monitors.
	 * For Retina displays (Mac) it expects double resulution images named 'originalName@2x.png'.
	 *
	 * @param imageName initial appearance
	 * @param tooltip may be null
	 * @param command action command to be used for action listeners (may be null)
	 */
	public HiDPIIconButton(String imageName, String tooltip, String command) {
		this(imageName, tooltip, command, 0);
		}

		/**
		 * Creates a button that, if image2 is given, toggles between two states indicated
		 * by two different button images. The button optimizes its size to for QuaQua
		 * and Substance look&feels and uses adequate higher resolution images on HiDPI monitors.
		 * For Retina displays (Mac) it expects double resolution images named 'originalName@2x.png'.
		 *
		 * @param imageName initial appearance
		 * @param tooltip may be null
		 * @param command action command to be used for action listeners (may be null)
		 * @param rotation 0, 90, 180, or 270 degrees in clockwise direction
		 */
	public HiDPIIconButton(String imageName, String tooltip, String command, int rotation) {
		super();

		mImageName = imageName;
		mRotation = rotation;
		updateIconSet();

		if (command != null)
			setActionCommand(command);

		setFocusable(false);

		if (tooltip != null)
			setToolTipText(tooltip);
		}

	/**
	 * In order for a button to support an animation, it need additional images to follow after its normal appearance.
	 * Animation image names are built from the button's imageName with appended '_0', '_1', etc.
	 * This methods loads and prepares the additional images needed for the button animation.
	 * @param animationImageCount number of animation images including the base image
	 */
	public void startAnimation(int animationImageCount) {
		if (mAnimationIcon == null) {
			mAnimationIcon = new Icon[animationImageCount];
			mAnimationIcon[0] = getIcon();
			for (int i=1; i<animationImageCount; i++) {
				int index = mImageName.lastIndexOf('.');
				String fileName = mImageName.substring(0, index)+"_"+i+mImageName.substring(index);
				mAnimationIcon[i] = HiDPIIcon.createIcon(fileName, 0, false);
				}
			}

		mAnimator = new Animator(animationImageCount);
		mAnimator.start();
		}

	public void stopAnimation() {
		if (mAnimator != null) {
			mAnimator.stop();
			}
		}

	public boolean isAnimating() {
		return mAnimator != null && mAnimator.isRunning();
		}

	private void updateIconSet() {
		if (mImageName != null) {
			setIcon(HiDPIIcon.createIcon(mImageName, mRotation, false));
			setDisabledIcon(HiDPIIcon.createIcon(mImageName, mRotation, true));

			Icon icon = getIcon();
			int w = Math.round(icon.getIconWidth() / HiDPIHelper.getPixelPerComponentSizeFactor()) + 2;
			int h = Math.round(icon.getIconHeight() / HiDPIHelper.getPixelPerComponentSizeFactor()) + 2;
			setPreferredSize(new Dimension(w, h));
			}
		}

	@Override
	public void updateUI() {
		updateIconSet();
		super.updateUI();
		}

	private class Animator {
		private static final long FRAME_RATE = 80L;
		private final int mFrameCount;
		private volatile boolean mRunning;

		public Animator(int frameCount) {
			mFrameCount = frameCount;
		}

		public void start() {
			if (!mRunning) {
				mRunning = true;
				new Thread(() -> {
					int state = 0;
					while (mRunning) {
						try {
							Thread.sleep(FRAME_RATE);
							final int frameNo = ++state % mFrameCount;
							SwingUtilities.invokeLater(() -> setIcon(mAnimationIcon[frameNo]));
							}
						catch (InterruptedException ie) {}
						}
					}).start();
				}
			}

		public void stop() {
			mRunning = false;
			}

		public boolean isRunning() {
			return mRunning;
		}
		}
	}
