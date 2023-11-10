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

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.image.AffineTransformOp;
import java.awt.image.BufferedImage;
import java.io.BufferedInputStream;

public class JImagePanelFixedSize extends JPanel {
	private BufferedImage   mImage;

    public JImagePanelFixedSize() {
		super();
		}

	public JImagePanelFixedSize(String fileName) {
		this(fileName, 1.0);
		}

    public JImagePanelFixedSize(String fileName, double scale) {
		super();

		readImage(fileName);
	    int width = (int)(scale * HiDPIHelper.scale(mImage.getWidth()));
	    int height = (int)(scale * HiDPIHelper.scale(mImage.getHeight()));
	    scaleImage(width, height);
		Dimension size = new Dimension(width, height);
		setMinimumSize(size);
		setMaximumSize(size);
		setPreferredSize(size);
		}

	public void setImage(String fileName) {
		readImage(fileName);
		repaint();
		}

	@Override public void paintComponent(Graphics g) {
		g.drawImage(mImage,0,0, this);
		}

	private void readImage(String fileName) {
		try {
			BufferedInputStream in = new BufferedInputStream(getClass().getResourceAsStream(fileName));
			mImage = ImageIO.read(in);
			}
		catch (Exception e) {}
		}

	public void scaleImage(int width, int height) {
		if (width != mImage.getWidth()
		 || height != mImage.getHeight()) {
			AffineTransform scaleTransform = AffineTransform.getScaleInstance(
					(double)width / mImage.getWidth(), (double)height / mImage.getHeight());
			AffineTransformOp bilinearScaleOp = new AffineTransformOp(scaleTransform, AffineTransformOp.TYPE_BILINEAR);
			mImage = bilinearScaleOp.filter(mImage, new BufferedImage(width, height, mImage.getType()));
			}
		}
	}