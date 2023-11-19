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

	/**
	 * Creates an image panel with the supplied image file name.
	 * After reading the image, it is automatically UI-scaled (enlarged on high-res screens).
	 * @param fileName
	 */
	public JImagePanelFixedSize(String fileName) {
		this(fileName, 1.0);
		}

	/**
	 * Creates an image panel with the supplied image file name.
	 * After reading the image, it is automatically UI-scaled (enlarged on high-res screens).
	 * The provided scale factor does not include UI-scaling and is used if source image
	 * has a high resolution to allow lossless UI-scaling.
	 * @param fileName
	 * @param scale typically smaller than 1.0 for high-res images that shall not loose during UI-scaling
	 */
    public JImagePanelFixedSize(String fileName, double scale) {
		super();
		readAndScaleImage(fileName, scale);
		}

	/**
	 * Creates a JImagePanelFixedSize without any scaling.
	 * If you need UI-scaling, the image has to be scaled in advance.
	 * @param image
	 */
	public JImagePanelFixedSize(BufferedImage image) {
		super();
		mImage = image;
		initializeSize(mImage.getWidth(), mImage.getHeight());
		}

	private void initializeSize(int width, int height) {
		Dimension size = new Dimension(width, height);
		setMinimumSize(size);
		setMaximumSize(size);
		setPreferredSize(size);
		}

	public void setImage(String fileName) {
		setImage(fileName, 1.0);
		}

	public void setImage(String fileName, double scale) {
		readAndScaleImage(fileName, scale);
		repaint();
		}

	@Override public void paintComponent(Graphics g) {
		if (mImage != null)
			g.drawImage(mImage,0,0, this);
		}

	private void readAndScaleImage(String fileName, double scale) {
		try {
			BufferedInputStream in = new BufferedInputStream(getClass().getResourceAsStream(fileName));
			mImage = ImageIO.read(in);
			if (mImage != null) {
				int width = (int)(scale * HiDPIHelper.scale(mImage.getWidth()));
				int height = (int)(scale * HiDPIHelper.scale(mImage.getHeight()));
				mImage = scaleImage(mImage, width, height);
				initializeSize(width, height);
				}
			}
		catch (Exception e) {}
		}

	public static BufferedImage scaleImage(BufferedImage image, int width, int height) {
		if (width == image.getWidth()
		 && height == image.getHeight())
			return image;

		AffineTransform scaleTransform = AffineTransform.getScaleInstance(
				(double)width / image.getWidth(), (double)height / image.getHeight());
		AffineTransformOp bilinearScaleOp = new AffineTransformOp(scaleTransform, AffineTransformOp.TYPE_BILINEAR);
		return bilinearScaleOp.filter(image, new BufferedImage(width, height, image.getType()));
		}
	}