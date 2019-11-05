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

import com.actelion.research.gui.LookAndFeelHelper;

import javax.swing.*;
import java.awt.*;

/**
 * Created by sandert on 04/12/15.
 */
public class HiDPIIconButton extends JButton {
	private String mImageName,mStyle;
	private int mRotation;

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
		this(imageName, tooltip, command, 0, "bevel");
		}

		/**
		 * Creates a button that, if image2 is given, toggles between two states indicated
		 * by two different button images. The button optimizes its size to for QuaQua
		 * and Substance look&feels and uses adequate higher resolution images on HiDPI monitors.
		 * For Retina displays (Mac) it expects double resulution images named 'originalName@2x.png'.
		 *
		 * @param imageName initial appearance
		 * @param tooltip may be null
		 * @param command action command to be used for action listeners (may be null)
		 * @param rotation 0, 90, 180, or 270 degrees in clockwise direction
		 * @param style one of "bevel","square",null (used for Quaqua LaF only)
		 */
	public HiDPIIconButton(String imageName, String tooltip, String command, int rotation, String style) {
		super();

		mImageName = imageName;
		mRotation = rotation;
		mStyle = style;
		updateIconSet();

		if (command != null)
			setActionCommand(command);

		setFocusable(false);

		if (tooltip != null)
			setToolTipText(tooltip);
		}

	private void updateIconSet() {
		if (mImageName != null) {
			setIcon(HiDPIHelper.createIcon(mImageName, mRotation));
			setDisabledIcon(HiDPIHelper.createDisabledIcon(mImageName, mRotation));

			Icon icon = getIcon();
			int w = Math.round(icon.getIconWidth() / HiDPIHelper.getRetinaScaleFactor()) + 2;
			int h = Math.round(icon.getIconHeight() / HiDPIHelper.getRetinaScaleFactor()) + 2;
			setPreferredSize(new Dimension(w, h));
			}
		}

	@Override
	public void updateUI() {
		updateIconSet();
		super.updateUI();
		}
	}
