package com.actelion.research.gui.hidpi;

import com.actelion.research.gui.LookAndFeelHelper;

import javax.swing.*;
import java.awt.*;

/**
 * Created by thomas on 07/12/15.
 */
public class HiDPIToggleButton extends JToggleButton {
	private String mImageName1,mImageName2;

	/**
	 * Creates a button that, if image2 is given, toggles between two states indicated
	 * by two different button images. The button optimizes its size to for QuaQua
	 * and Substance look&feels and uses adequate higher resolution images on HiDPI monitors.
	 * For Retina displays (Mac) it expects double resolution images named 'originalName@2x.png'.
	 *
	 * @param imageName1 initial appearance
	 * @param imageName2 toggled appearance
	 * @param tooltip may be null
	 * @param command action command to be used for action listeners (may be null)
	 */
	public HiDPIToggleButton(String imageName1, String imageName2, String tooltip, String command) {
		super();

		mImageName1 = imageName1;
		mImageName2 = imageName2;

		updateIconSet();

		if (command != null)
			setActionCommand(command);

		setFocusable(false);

		if (tooltip != null)
			setToolTipText(tooltip);
		}

	private void updateIconSet() {
		if (mImageName1 != null) {
			setIcon(HiDPIHelper.createIcon(mImageName1, 0));
			setSelectedIcon(HiDPIHelper.createIcon(mImageName2, 0));
			setDisabledIcon(HiDPIHelper.createDisabledIcon(mImageName1, 0));
			setDisabledSelectedIcon(HiDPIHelper.createDisabledIcon(mImageName2, 0));

			Icon icon = getIcon();
			int w = icon.getIconWidth() / (int)HiDPIHelper.getRetinaScaleFactor() + 2;
			int h = icon.getIconHeight() / (int)HiDPIHelper.getRetinaScaleFactor() + 2;
			if (LookAndFeelHelper.isQuaQua()) {
				w += 2;
				h += 2;
				putClientProperty("Quaqua.Component.visualMargin", new Insets(1, 1, 1, 1));
				putClientProperty("Quaqua.Button.style", "bevel");
				}
			setPreferredSize(new Dimension(w, h));
			}
		}

	@Override
	public void updateUI() {
		updateIconSet();
		super.updateUI();
		}
	}
