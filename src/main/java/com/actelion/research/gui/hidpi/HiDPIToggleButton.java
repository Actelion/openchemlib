package com.actelion.research.gui.hidpi;

import javax.swing.*;
import java.awt.*;

/**
 * Created by thomas on 07/12/15.
 */
public class HiDPIToggleButton extends JToggleButton {
	private final String mImageName1,mImageName2;

	/**
	 * Creates a button that, if image2 is given, toggles between two states indicated
	 * by two different button images. The button optimizes its size to for QuaQua
	 * and Substance look&feels and uses adequate higher resolution images on HiDPI monitors.
	 * For Retina displays (Mac) it expects double resolution images named 'originalName@2x.png'.
	 *
	 * @param imageName1 initial appearance
	 * @param imageName2 toggled appearance (may be null)
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
			setIcon(HiDPIIcon.createIcon(mImageName1, 0, false));
			setSelectedIcon(HiDPIIcon.createIcon(mImageName2 != null ? mImageName2 : mImageName1, 0, false));
			setDisabledIcon(HiDPIIcon.createIcon(mImageName1, 0, true));
			setDisabledSelectedIcon(HiDPIIcon.createIcon(mImageName2 != null ? mImageName2 : mImageName1, 0, true));

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
	}
