package com.actelion.research.gui.hidpi;

import com.actelion.research.gui.LookAndFeelHelper;

import javax.swing.*;
import java.awt.*;

/**
 * Created by thomas on 07/12/15.
 */
public class HiDPIIcon extends ImageIcon {

	public HiDPIIcon(final Image image) {
		super(image);
	}

	public synchronized void paintIcon(Component c, Graphics g, int x, int y) {
		if (HiDPIHelper.getRetinaScaleFactor() == 2) {
			Image image = getImage();
			int width = image.getWidth(null) / 2;
			int height = image.getHeight(null) / 2;

			// because of double size image, x and y are too small and need to be corrected
			// (in case of aqua x & y are 1 - size/2, in case of new substance x & y are 0)
			x += width / 2;
			y += height / 2;

			if (LookAndFeelHelper.isNewSubstance()) {
				// for some reason y needs to be adjusted by 1
				y++;
				}

			g.drawImage(image, x, y, width, height, null);
			}
		else {
			Image image = getImage();
			g.drawImage(image, x, y, null);
			}
		}
	}
