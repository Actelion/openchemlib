package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.GenericImage;
import com.actelion.research.gui.swing.SwingImage;
import javafx.embed.swing.SwingFXUtils;
import javafx.scene.image.Image;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.net.URL;

public class FXImage implements GenericImage {
	private Image mFXImage;
	private BufferedImage mAWTImage;

	public FXImage(String name) {
		URL url = SwingImage.class.getResource("/images/" + name);
		if (url == null)
			throw new RuntimeException("Could not find: " + name);

		try {
			mAWTImage = ImageIO.read(url);
		}
		catch (IOException ioe) {}

// We don't use an FX Image or WritableImage, because there is no easy way to scale them with retaining transparency
//		Image image = new Image("images/" + name);
//		mImage = new WritableImage(image.getPixelReader(), (int)image.getWidth(), (int)image.getHeight());
	}

	public FXImage(int width, int height) {
		mAWTImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
	}

	@Override
	public int getWidth() {
		return Math.round((float)mAWTImage.getWidth());
	}

	@Override
	public int getHeight() {
		return Math.round((float)mAWTImage.getHeight());
	}

	@Override
	public int getRGB(int x, int y) {
		return mAWTImage.getRGB(x, y);
	}

	@Override
	public void setRGB(int x, int y, int argb) {
		mAWTImage.setRGB(x, y, argb);
		mFXImage = null;
	}

	@Override
	public Image get() {
		if (mFXImage == null) {
			mFXImage = SwingFXUtils.toFXImage(mAWTImage, null);
			}
		return mFXImage;
		}

	@Override
	public void scale(int width, int height) {
		java.awt.Image scaledImage = mAWTImage.getScaledInstance(width, height, java.awt.Image.SCALE_SMOOTH);
		mAWTImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		Graphics g = mAWTImage.createGraphics();
		g.drawImage(scaledImage, 0, 0, null);
		g.dispose();
		mFXImage = null;

// This method of scaling looses the transparency
//		ImageView imageView = new ImageView(mImage);
//		imageView.setPreserveRatio(false);
//		imageView.setFitWidth(width);
//		imageView.setFitHeight(height);
//		mImage = imageView.snapshot(null, null);
		}
	}
