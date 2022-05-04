package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.GenericImage;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.image.WritableImage;

public class FXImage implements GenericImage {
	private Image mImage;

	public FXImage(String name) {
//		URL url = FXImage.class.getResource("/images/" + name);
//		if (url == null)
//			throw new RuntimeException("Could not find: " + name);
//
//		mImage = new Image(url.getPath());
		mImage = new Image("images/" + name);
	}

	public FXImage(int width, int height) {
		mImage = new WritableImage(width, height);
	}

	@Override
	public int getWidth() {
		return Math.round((float)mImage.getWidth());
	}

	@Override
	public int getHeight() {
		return Math.round((float)mImage.getHeight());
	}

	@Override
	public void setRGB(int x, int y, int argb) {
		((WritableImage)mImage).getPixelWriter().setArgb(x, y, argb);
	}

	@Override
	public Image get() {
		return mImage;
		}

	@Override
	public void scale(int width, int height) {
		ImageView imageView = new ImageView(mImage);
		imageView.setPreserveRatio(false);
		imageView.setFitWidth(width);
		imageView.setFitHeight(height);
		mImage = imageView.snapshot(null, null);
		}
	}
