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

	private static final int[] ICON_SPOT_COLOR = {   // original spot colors used in icon images (bright L&F)
			0x00503CB4, 0x00000000 };
	private static final int[] DARK_LAF_SPOT_COLOR = {   // default replacement spot colors for dark L&F)
			0x00B4A0FF, 0x00E0E0E0 };

	private final Image mUnscaledImage;
	private static int[] sSpotColor = null;

	public HiDPIIcon(Image scaled, Image unscaled) {
		super(scaled);
		mUnscaledImage = unscaled;
		}

	public static Icon createIcon(String fileName, int rotation, boolean isDisabled) {
		GenericImage image = createIconImage(fileName);
		adaptForLookAndFeel(image);
		rotate(image, rotation);
		if (isDisabled)
			HiDPIHelper.disableImage(image);
		capCorners(image);
		Image unscaled = (Image)image.get();
		return new HiDPIIcon((Image)scale(image).get(), unscaled);
		}

	public static void setIconSpotColors(int[] rgb) {
		sSpotColor = rgb;
	}

	public static int[] getThemeSpotRGBs() {
		return (sSpotColor != null) ? sSpotColor
				: LookAndFeelHelper.isDarkLookAndFeel() ? DARK_LAF_SPOT_COLOR
				: ICON_SPOT_COLOR;
	}

	public static void adaptForLookAndFeel(GenericImage image) {
		if (sSpotColor != null)
			replaceSpotColors(image, sSpotColor);
		else if (LookAndFeelHelper.isDarkLookAndFeel())
			replaceSpotColors(image, DARK_LAF_SPOT_COLOR);
	}

	private static void replaceSpotColors(GenericImage image, int[] altRGB) {
		for (int x=0; x<image.getWidth(); x++) {
			for (int y=0; y<image.getHeight(); y++) {
				int argb = image.getRGB(x, y);
				int rgb = argb & 0x00FFFFFF;
				for (int i=0; i<ICON_SPOT_COLOR.length; i++) {
					if (rgb == ICON_SPOT_COLOR[i]) {
						image.setRGB(x, y, (0xFF000000 & argb) + altRGB[i]);
						break;
					}
				}
			}
		}
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

	public static GenericImage scale(GenericImage image) {
		if (!mustScale())
			return image;

		float scale = HiDPIHelper.getIconScaleFactor();
		if (useDoubleImage())
			scale *= 0.5f;

		image.scale(Math.round(scale * image.getWidth()), Math.round(scale * image.getHeight()));
		return image;
		}

	private static boolean mustScale() {
		return (HiDPIHelper.getIconScaleFactor() > ICON_SCALE_LIMIT_1
			 && HiDPIHelper.getIconScaleFactor() < ICON_SCALE_LIMIT_2)
			|| HiDPIHelper.getIconScaleFactor() > ICON_SCALE_LIMIT_3;
		}


	private static boolean useDoubleImage() {
		return HiDPIHelper.getIconScaleFactor() > ICON_SCALE_LIMIT_1;
		}

	@Override
	public synchronized void paintIcon(Component c, Graphics g, int x, int y) {
		if (HiDPIHelper.getPixelPerComponentSizeFactor() != 1f) {
			if (Platform.isMacintosh()) {
				Image image = getImage();
				int width = Math.round(image.getWidth(null) / HiDPIHelper.getPixelPerComponentSizeFactor());
				int height = Math.round(image.getHeight(null) / HiDPIHelper.getPixelPerComponentSizeFactor());

				// because of larger size image, x and y are too small and need to be corrected
				// (in case of aqua x & y are 1 - size/2, in case of new substance and Windows x & y are 0)
				x += Math.round(width / HiDPIHelper.getPixelPerComponentSizeFactor());
				y += Math.round(height / HiDPIHelper.getPixelPerComponentSizeFactor());

				// for some reason y needs to be adjusted by 1 in new substance
				if (LookAndFeelHelper.isNewSubstance()
				 || LookAndFeelHelper.isRadiance())
					y++;

				g.drawImage(image, x, y, width, height, null);
				}
			else if (Platform.isWindows()) {
				// we have a centered viewport on a rectangle of the size of the image
				Image image = getImage();
				float width = image.getWidth(null) / HiDPIHelper.getPixelPerComponentSizeFactor();
				float height = image.getHeight(null) / HiDPIHelper.getPixelPerComponentSizeFactor();

//				System.out.println("HiDPIIcon.paintIcon() x:"+x+" icw:"+getIconWidth()+" w:"+width+" imw:"+image.getWidth(null)+" pixFac:"+HiDPIHelper.getPixelPerComponentSizeFactor());

				// because of larger size image, x and y are too small and need to be corrected
				// (in case of aqua x & y are 1 - size/2, in case of new substance and Windows x & y are 0)
				int dx = Math.round((image.getWidth(null) - width) / 2f);
				int dy = Math.round((image.getHeight(null) - height) / 2f);

				// for some reason y needs to be adjusted by 1 in new substance
				if (LookAndFeelHelper.isNewSubstance()
				 || LookAndFeelHelper.isRadiance())
					dy++;

				g.drawImage(image, dx, dy, Math.round(width), Math.round(height), null);
//				g.drawImage(mUnscaledImage, dx, dy, Math.round(width), Math.round(height), null);
				}
			}
		else {
			// for some reason y needs to be adjusted by 1 in new substance
			if (LookAndFeelHelper.isNewSubstance()
			 || LookAndFeelHelper.isRadiance())
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
