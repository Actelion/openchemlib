package com.actelion.research.gui.swing;

import com.actelion.research.gui.generic.GenericImage;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.net.URL;

public class SwingImage implements GenericImage {
	private BufferedImage mImage;

	public SwingImage(BufferedImage image) {
		mImage = image;
		}

	public SwingImage(String name) {
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
	public int getRGB(int x, int y) {
		return mImage.getRGB(x, y);
		}

	@Override
	public void setRGB(int x, int y, int argb) {
		mImage.setRGB(x, y, argb);
		}

	@Override
	public BufferedImage get() {
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
		java.awt.Image scaledImage = mImage.getScaledInstance(width, height, java.awt.Image.SCALE_SMOOTH);
		mImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		Graphics g = mImage.createGraphics();
		g.drawImage(scaledImage, 0, 0, null);
		g.dispose();
		}
	}
