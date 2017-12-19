package com.actelion.research.gui.hidpi;

import com.actelion.research.gui.LookAndFeelHelper;
import com.actelion.research.util.ColorHelper;
import com.actelion.research.util.Platform;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.lang.reflect.Field;
import java.net.URL;

public class HiDPIHelper {
	private static final int COLOR1_BRIGHT_LAF = 0x00503CB4;	// main button color in supplied images
	private static final int COLOR2_BRIGHT_LAF = 0x00000000;	// second button color in supplied images

	private static final int COLOR1_DARK_LAF = 0x00B4A0FF;	// main button color in supplied images
	private static final int COLOR2_DARK_LAF = 0x00E0E0E0;	// second button color in supplied images

	private static final float ICON_SCALE_LIMIT = 1.2f; // custom dpi scale factors smaller than this will be neglected

	// This is an Apple only solution and needs to be adapted to support high-res displays of other vendors
	private static float sRetinaFactor = -1f;
	private static float sUIScaleFactor = -1f;

	/**
	 * Macintosh retina display support for Java 7 and newer.
	 *
	 * @return 1.0 on standard resolution devices and 2.0 for retina screens
	 */
	public static float getRetinaScaleFactor() {
		/* with Apple-Java-6 this was:
		Object sContentScaleFactorObject = Toolkit.getDefaultToolkit().getDesktopProperty("apple.awt.contentScaleFactor");
		private static final float sRetinaFactor = (sContentScaleFactorObject == null) ? 1f : ((Float)sContentScaleFactorObject).floatValue();
		*/
		if (sRetinaFactor != -1f)
			return sRetinaFactor;

		sRetinaFactor = 1f;

		GraphicsEnvironment env = GraphicsEnvironment.getLocalGraphicsEnvironment();
		final GraphicsDevice device = env.getDefaultScreenDevice();

		try {
			Field field = device.getClass().getDeclaredField("scale");
			if (field != null) {
				field.setAccessible(true);
				Object scale = field.get(device);

				if (scale instanceof Integer)
					sRetinaFactor = (Integer) scale;
				else
					System.out.println("Unexpected content scale (not 1 nor 2): "+scale.toString());
				}
			}
		catch (Throwable e) {}

/*	the above code gives WARNING under Java 9:
 			WARNING: An illegal reflective access operation has occurred
 			WARNING: All illegal access operations will be denied in a future release

			If we know, we are on a Mac, we could do something like:

		if (device instanceof CGraphicsDevice) {	// apple.awt.CGraphicsDevice
			final CGraphicsDevice cgd = (CGraphicsDevice)device;

			// this is the missing correction factor, it's equal to 2 on HiDPI a.k.a. Retina displays
			final int scaleFactor = cgd.getScaleFactor();

			// now we can compute the real DPI of the screen
			final double realDPI = scaleFactor * (cgd.getXResolution() + cgd.getYResolution()) / 2;
			}*/

		return sRetinaFactor;
		}

	/**
	 * For Windows and Linux this method returns the user defined UI scaling factor.
	 * This is done by judging from the size of the UIManager's Label.font
	 * and comparing it to the unscaled default (13). Typically this factor is larger
	 * than 1.0 on HiDPI devices. For this method to work the Look&Feel must consider
	 * the OS provided setting and scale its fonts accordingly (Substance LaF does).<br>
	 * On the Macintosh this factor is usually 1.0, because HiDPI device support uses
	 * a different mechanism (see getRetinaScaleFactor()).
	 * @return typically 1.0 or 1.25, 1.5, ...
	 */
	public static float getUIScaleFactor() {
		if (sUIScaleFactor == -1)
			sUIScaleFactor = Platform.isMacintosh() ? 1f : (float) UIManager.getFont("Label.font").getSize() / 12f;

		return sUIScaleFactor;
		}

	/**
	 * This is a convenience method that scales the passed int value with getUIScaleFactor()
	 * and returns the rounded result.
	 * @param value
	 * @return
	 */
	public static int scale(int value) {
		return Math.round(getUIScaleFactor() * value);
		}

	/**
	 * This is a convenience method that scales the passed int value with getUIScaleFactor()
	 * and with getRetinaScaleFactor() and returns the rounded result.
	 * @param value
	 * @return
	 */
	public static int scaleRetinaAndUI(int value) {
		return Math.round(getUIScaleFactor() * getRetinaScaleFactor() * value);
		}

	public static Icon createIcon(String fileName, int rotation) {
		return new HiDPIIcon(scale(capCorners(rotate(createLaFCompatibleImage(fileName), rotation))));
		}

	public static Icon createDisabledIcon(String fileName, int rotation) {
		return new HiDPIIcon(scale(capCorners(rotate(createDisabledImage(fileName), rotation))));
		}

	private static String getDoubleResolutionFileName(String fileName) {
		int index = fileName.lastIndexOf('.');
		return fileName.substring(0, index).concat("@2x").concat(fileName.substring(index));
		}

	private static BufferedImage rotate(BufferedImage image, int rotation) {
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

	/**
	 * Creates an image from the fileName. On HiDPI devices this is a high
	 * resolution image. If the current look&feel is dark, then colors are adapted
	 * for optimal contrast.
	 * @param fileName
	 * @return
	 */
	public static BufferedImage createImage(String fileName) {
		if (getRetinaScaleFactor() == 2 || getUIScaleFactor() > ICON_SCALE_LIMIT)
			fileName = getDoubleResolutionFileName(fileName);

		URL url = HiDPIIconButton.class.getResource("/images/" + fileName);
		if (url == null)
			throw new RuntimeException("Could not find: " + fileName);

		try {
			BufferedImage image = ImageIO.read(url);
			return image;
		}
		catch (IOException ioe) { return null; }
	}

	/**
	 * Creates an image from the fileName. On HiDPI devices this is a high
	 * resolution image. If the current look&feel is dark, then colors are adapted
	 * for optimal contrast.
	 * @param fileName
	 * @return
	 */
	public static BufferedImage createLaFCompatibleImage(String fileName) {
		BufferedImage image = createImage(fileName);
		if (image != null && LookAndFeelHelper.isDarkLookAndFeel())
			brightenImage(image);
		return image;
	}

	private static BufferedImage createDisabledImage(String fileName) {
		BufferedImage image = createLaFCompatibleImage(fileName);

		Color gray = LookAndFeelHelper.isDarkLookAndFeel() ?
				ColorHelper.brighter(UIManager.getColor("Panel.background"), 0.8f)
			  : ColorHelper.darker(UIManager.getColor("Panel.background"), 0.8f);
		int grayRGB = 0x00FFFFFF & gray.getRGB();

		for (int x=0; x<image.getWidth(); x++) {
			for (int y=0; y<image.getHeight(); y++) {
				int argb = image.getRGB(x, y);
				image.setRGB(x, y, (0xFF000000 & argb) + grayRGB);
			}
		}
		return image;
	}

	private static Image scale(BufferedImage image) {
		float scale = getUIScaleFactor();
		if (scale > ICON_SCALE_LIMIT)   // in this case we have double size images
			return image.getScaledInstance(Math.round(0.5f * scale * image.getWidth()),
										   Math.round(0.5f * scale * image.getHeight()), Image.SCALE_SMOOTH);
		else
			return image;
	}

	private static BufferedImage capCorners(BufferedImage image) {
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

	private static void brightenImage(BufferedImage image) {
		for (int x=0; x<image.getWidth(); x++) {
			for (int y=0; y<image.getHeight(); y++) {
				int argb = image.getRGB(x, y);
				int rgb = argb & 0x00FFFFFF;
				int color = (rgb == COLOR1_BRIGHT_LAF) ? COLOR1_DARK_LAF
						: (rgb == COLOR2_BRIGHT_LAF) ? COLOR2_DARK_LAF : rgb;
				image.setRGB(x, y, (0xFF000000 & argb) + color);
			}
		}
	}
}
