package com.actelion.research.gui.hidpi;

import com.actelion.research.gui.LookAndFeelHelper;
import com.actelion.research.gui.generic.GenericImage;
import com.actelion.research.util.ColorHelper;
import com.actelion.research.util.Platform;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.lang.reflect.Field;

public class HiDPIHelper {
	// This is an Apple only solution and needs to be adapted to support high-res displays of other vendors
	private static float sRetinaFactor = -1f;
	private static float sUIScaleFactor = -1f;

	/**
	 * Macintosh retina display support for Java 7 and newer.
	 *
	 * @return 1.0 on standard resolution devices and 2.0 for retina screens
	 */
	public static float getRetinaScaleFactor() {
		if (!Platform.isMacintosh())
			return 1f;

		if (sRetinaFactor != -1f)
			return sRetinaFactor;

		sRetinaFactor = 1f;

		GraphicsEnvironment env = GraphicsEnvironment.getLocalGraphicsEnvironment();
		final GraphicsDevice device = env.getDefaultScreenDevice();

		if (System.getProperty("java.version").startsWith("1.")) {  // for JRE8 and earlier
			try {
				Field field = device.getClass().getDeclaredField("scale");
				if (field != null) {
					field.setAccessible(true);
					Object scale = field.get(device);

					if (scale instanceof Integer)
						sRetinaFactor = (Integer) scale;
					else
						System.out.println("Unexpected content scale (not 1 nor 2): " + scale.toString());
					}
				}
			catch (Throwable e) {}
			}
		else {
			GraphicsDevice sd = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice();
			sRetinaFactor = (float) sd.getDefaultConfiguration().getDefaultTransform().getScaleX();
		}
		return sRetinaFactor;
	}

	/**
	 * For Windows and Linux this method returns a (hopefully) reasonable UI scaling factor.
	 * If DataWarrior was launched with a user supplied system property 'dpifactor', then we use
	 * this value. Otherwise, we take the screen resolution devided by 96, which finally results
	 * in a scaling factor of 1.0 at 96 dpi monitor resolution. If we use factor*12 point system
	 * fonts, then these are about 3 mm high on monitors independent of their resolution.
	 * On the Macintosh this factor is 1.0, because HiDPI device support uses
	 * a different mechanism on OSX (see getRetinaScaleFactor()).
	 * @return typically 1.0 or 1.25, 1.5, ...
	 */
	public static float getUIScaleFactor() {
		if (sUIScaleFactor == -1) {
			float f = 0;
			String dpiFactor = System.getProperty("dpifactor");
			if (dpiFactor != null)
				try { f = Float.parseFloat(dpiFactor); } catch (NumberFormatException nfe) {}
			if (f != 0)
				sUIScaleFactor = f;
			else if (Platform.isMacintosh()) {
				sUIScaleFactor = 1.0f;
				}
			else if (Platform.isWindows()) {
				try {
					// with JRE8 we used (float)UIManager.getFont("Label.font").getSize() / 12f
					sUIScaleFactor = Toolkit.getDefaultToolkit().getScreenResolution() / 96f;
					}
				catch (HeadlessException hle) {
					sUIScaleFactor = 1.0f;
					}
				}
			else {  // Linux; Toolkit.getDefaultToolkit().getScreenResolution() always returns 1.0
				try {
					sUIScaleFactor = 1.0f;  // default in case of error
					Process process = Runtime.getRuntime().exec("xrdb -q");
					process.waitFor();
					BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()));
					String line;
					while ((line = br.readLine()) != null) {
						if (line.startsWith("Xft.dpi:")) {
							int dpi = Integer.parseInt(line.substring(8).trim());
							sUIScaleFactor = dpi / 96f;
							break;
							}
						}
					br.close();
					}
				catch (Exception ioe) {}
				}
			if (sUIScaleFactor < 1.1f)
				sUIScaleFactor = 1.0f;
			}
		return sUIScaleFactor;
		}

	public static void setUIScaleFactor(float factor) {
		sUIScaleFactor = factor;
	}

	/**
	 * This is a convenience method that scales the passed int value with getUIScaleFactor()
	 * and returns the rounded result.
	 * @param value
	 * @return
	 */
	public static int scale(float value) {
		return Math.round(getUIScaleFactor() * value);
		}

	/**
	 * This is a convenience method that scales the passed int value with getUIScaleFactor()
	 * and with getRetinaScaleFactor() and returns the rounded result.
	 * @param value
	 * @return
	 */
	public static int scaleRetinaAndUI(float value) {
		return Math.round(getUIScaleFactor() * getRetinaScaleFactor() * value);
		}

	/**
	 * If the current look&feel is dark, then colors are adapted for optimal contrast.
	 * @param image
	 * @return
	 */

	public static void disableImage(GenericImage image) {
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
	}
}
