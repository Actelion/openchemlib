package com.actelion.research.gui.hidpi;

import com.actelion.research.gui.LookAndFeelHelper;
import com.actelion.research.util.Platform;

import javax.swing.*;
import java.awt.*;

/**
 * Created by thomas on 07/12/15.
 */
public class HiDPIIcon extends ImageIcon {

	public HiDPIIcon(final Image image, final boolean isScaled) {
		super(image);
		}

	public synchronized void paintIcon(Component c, Graphics g, int x, int y) {
		if (HiDPIHelper.getRetinaScaleFactor() != 1f) {
			Image image = getImage();
			int width = Math.round(image.getWidth(null) / HiDPIHelper.getRetinaScaleFactor());
			int height = Math.round(image.getHeight(null) / HiDPIHelper.getRetinaScaleFactor());

//			System.out.println("HiDPIIcon.paintIcon() x:"+x+" y:"+y);

			// because of double size image, x and y are too small and need to be corrected
			// (in case of aqua x & y are 1 - size/2, in case of new substance and Windows x & y are 0)
			if (Platform.isMacintosh()) {
				x += Math.round(width / HiDPIHelper.getRetinaScaleFactor());
				y += Math.round(height / HiDPIHelper.getRetinaScaleFactor());
				}
			else {
				x += (image.getWidth(null) - width) / 2;
				y += (image.getHeight(null) - height) / 2;
				}

			// for some reason y needs to be adjusted by 1 in new substance
			if (LookAndFeelHelper.isNewSubstance())
				y++;

			g.drawImage(image, x, y, width, height, null);
			}
		else {
			// for some reason y needs to be adjusted by 1 in new substance
			if (LookAndFeelHelper.isNewSubstance())
				y++;

			Image image = getImage();
			g.drawImage(image, x, y, null);
			}
		}
	}
