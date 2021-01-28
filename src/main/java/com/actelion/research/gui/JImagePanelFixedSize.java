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

import javax.swing.*;
import java.awt.*;
import java.io.BufferedInputStream;

public class JImagePanelFixedSize extends JPanel {
	private Image		mImage;

    public JImagePanelFixedSize() {
		super();
		}

    public JImagePanelFixedSize(String fileName) {
		super();

		readImage(fileName);
		Dimension size = new Dimension(mImage.getWidth(this), mImage.getHeight(this));
		setMinimumSize(size);
		setMaximumSize(size);
		setPreferredSize(size);
		}

	public void setImage(String fileName) {
		readImage(fileName);
		repaint();
		}

	@Override public void paintComponent(Graphics g) {
		g.drawImage(mImage,0,0,this);
		}

	private void readImage(String fileName) {
		Toolkit tk = Toolkit.getDefaultToolkit();
		mImage = null;
		try {
			BufferedInputStream in=new BufferedInputStream(getClass().getResourceAsStream(fileName));
			byte[] imageData = new byte[102400]; // make it as big as necessary
			in.read(imageData);
			mImage = tk.createImage(imageData);
			}
		catch (Exception e) {}

		MediaTracker t = new MediaTracker(this);
		t.addImage(mImage, 0);
		try {
			t.waitForAll();
			}
		catch (InterruptedException e) {}
		}
	}