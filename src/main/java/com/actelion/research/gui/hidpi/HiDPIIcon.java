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
			x = 1;
			y = 1;
			if (LookAndFeelHelper.isQuaQua()) {
				x++;
				y++;
				}
			else if (LookAndFeelHelper.isNewSubstance()) {
				x--;
				y--;
				}

			Image image = getImage();
			int width = image.getWidth(null);
			int height = image.getHeight(null);
			g.drawImage(image, x, y, width / 2, height / 2, null);
			}
		else {
			x = 1;
			y = 1;
			if (LookAndFeelHelper.isQuaQua()) {
				x++;
				y++;
				}
			else if (LookAndFeelHelper.isNewSubstance()) {
				x--;
				y--;
				}

			Image image = getImage();
			g.drawImage(image, x, y, null);
			}
		}
	}
