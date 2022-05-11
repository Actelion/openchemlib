package com.actelion.research.gui.hidpi;

import com.actelion.research.gui.LookAndFeelHelper;
import com.actelion.research.gui.generic.GenericImage;
import com.actelion.research.gui.swing.SwingImage;
import com.actelion.research.util.Platform;

import javax.swing.*;
import java.awt.*;

/**
 * Created by thomas on 07/12/15.
 */
public class HiDPIIcon extends ImageIcon {
	private static final float ICON_SCALE_LIMIT_1 = 1.1f; // custom dpi scale factors smaller than this will be neglected
	private static final float ICON_SCALE_LIMIT_2 = 1.9f; // custom dpi scale factors between ICON_SCALE_LIMIT_2 and ICON_SCALE_LIMIT_3
	private static final float ICON_SCALE_LIMIT_3 = 2.1f; // use larger image, but won't be scaled

	public HiDPIIcon(final Image image) {
		super(image);
		}

	public static Icon createIcon(String fileName, int rotation, boolean isDisabled) {
		GenericImage image = createIconImage(fileName);
		HiDPIHelper.adaptForLookAndFeel(image);
		rotate(image, rotation);
		if (isDisabled)
			HiDPIHelper.disableImage(image);
		capCorners(image);
		return new HiDPIIcon((Image)(mustScale() ? scale(image) : image).get());
		}

	/**
	 * Creates an image from the fileName. On HiDPI devices this is a high
	 * resolution image. If the current look&feel is dark, then colors are adapted
	 * for optimal contrast.
	 * @param fileName
	 * @return
	 */
	public static SwingImage createIconImage(String fileName) {
		return new SwingImage(useDoubleImage() ? getDoubleResolutionFileName(fileName) : fileName);
		}

	private static String getDoubleResolutionFileName(String fileName) {
		int index = fileName.lastIndexOf('.');
		return fileName.substring(0, index).concat("@2x").concat(fileName.substring(index));
		}

	public static float getIconScaleFactor() {
		if (!mustScale())
			return 1f;

		float scale = HiDPIHelper.getUIScaleFactor() * HiDPIHelper.getRetinaScaleFactor();
		if (useDoubleImage())
			scale *= 0.5f;

		return scale;
		}

	public static GenericImage scale(GenericImage image) {
		float scale = getIconScaleFactor();
		image.scale(Math.round(scale * image.getWidth()), Math.round(scale * image.getHeight()));
		return image;
		}

	private static boolean mustScale() {
		return (HiDPIHelper.getUIScaleFactor() * HiDPIHelper.getRetinaScaleFactor() > ICON_SCALE_LIMIT_1
				&& HiDPIHelper.getUIScaleFactor() * HiDPIHelper.getRetinaScaleFactor() < ICON_SCALE_LIMIT_2)
				|| HiDPIHelper.getUIScaleFactor() * HiDPIHelper.getRetinaScaleFactor() < ICON_SCALE_LIMIT_3;
		}


	private static boolean useDoubleImage() {
		return HiDPIHelper.getUIScaleFactor() * HiDPIHelper.getRetinaScaleFactor() > ICON_SCALE_LIMIT_1;
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

	private static GenericImage capCorners(GenericImage image) {
		image.setRGB(0, 0, 0x00000000);
		image.setRGB(1, 0, 0x00000000);
		image.setRGB(0, 1, 0x00000000);
		image.setRGB(image.getWidth()-2, 0, 0x00000000);
		image.setRGB(image.getWidth()-1, 0, 0x00000000);
		image.setRGB(image.getWidth()-1, 1, 0x00000000);
		image.setRGB(0, image.getHeight()-1, 0x00000000);
		image.setRGB(1, image.getHeight()-1, 0x00000000);
		image.setRGB(0, image.getHeight()-2, 0x00000000);
		image.setRGB(image.getWidth()-2, image.getHeight()-1, 0x00000000);
		image.setRGB(image.getWidth()-1, image.getHeight()-1, 0x00000000);
		image.setRGB(image.getWidth()-1, image.getHeight()-2, 0x00000000);
		return image;
		}

	private static GenericImage rotate(GenericImage image, int rotation) {
		int size = image.getHeight();
		int max = size-1;
		if (rotation == 90) {
			for (int x=0; x<size/2; x++) {
				for (int y=0; y<size/2; y++) {
					int argb = image.getRGB(x, y);
					image.setRGB(x, y, image.getRGB(y, max-x));
					image.setRGB(y, max-x, image.getRGB(max-x, max-y));
					image.setRGB(max-x, max-y, image.getRGB(max-y, x));
					image.setRGB(max-y, x, argb);
				}
			}
		}
		if (rotation == 180) {
			for (int x=0; x<size/2; x++) {
				for (int y=0; y<size; y++) {
					int argb = image.getRGB(x, y);
					image.setRGB(x, y, image.getRGB(max-x, max-y));
					image.setRGB(max-x, max-y, argb);
				}
			}
		}
		if (rotation == 270) {
			for (int x=0; x<size/2; x++) {
				for (int y=0; y<size/2; y++) {
					int argb = image.getRGB(x, y);
					image.setRGB(x, y, image.getRGB(max-y, x));
					image.setRGB(max-y, x, image.getRGB(max-x, max-y));
					image.setRGB(max-x, max-y, image.getRGB(y, max-x));
					image.setRGB(y, max-x, argb);
				}
			}
		}
		return image;
	}
}
