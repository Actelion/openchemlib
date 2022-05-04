package com.actelion.research.gui.swing;

import com.actelion.research.gui.generic.GenericImage;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.net.URL;

public class SwingImage implements GenericImage {
	private Image mImage;

	public SwingImage(String name) {
		// Once double resolution is available use HiDPIHelper.createImage() !!!

		URL url = SwingImage.class.getResource("/images/" + name);
		if (url == null)
			throw new RuntimeException("Could not find: " + name);

		try {
			mImage = ImageIO.read(url);
			}
		catch (IOException ioe) {}
		}

	public SwingImage(int width, int height) {
		mImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		}

	@Override
	public void setRGB(int x, int y, int argb) {
		((BufferedImage)mImage).setRGB(x, y, argb);
		}

	@Override
	public Image get() {
		return mImage;
		}

	@Override
	public int getWidth() {
		return mImage.getWidth(null);
		}

	@Override
	public int getHeight() {
		return mImage.getHeight(null);
		}

	@Override
	public void scale(int width, int height) {
		mImage = mImage.getScaledInstance(width, height, Image.SCALE_SMOOTH);
		}
	}
