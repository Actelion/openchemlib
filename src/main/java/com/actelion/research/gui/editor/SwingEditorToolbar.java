package com.actelion.research.gui.editor;

import com.actelion.research.gui.generic.GenericCanvas;
import com.actelion.research.gui.generic.GenericDrawContext;
import com.actelion.research.gui.swing.SwingDrawContext;
import com.actelion.research.gui.swing.SwingMouseHandler;

import javax.swing.*;
import java.awt.*;

public class SwingEditorToolbar extends JPanel implements GenericCanvas {
	private GenericEditorToolbar mGenericToolbar;

	public SwingEditorToolbar(SwingEditorArea swingEditorArea) {
		mGenericToolbar = new GenericEditorToolbar(this, swingEditorArea.getGenericDrawArea());

		SwingMouseHandler mouseHandler = new SwingMouseHandler(mGenericToolbar);
		addMouseListener(mouseHandler);
		addMouseMotionListener(mouseHandler);
		mouseHandler.addListener(mGenericToolbar);

		int w = mGenericToolbar.getWidth();
		int h = mGenericToolbar.getHeight();
		setMinimumSize(new Dimension(w,h));
		setPreferredSize(new Dimension(w,h));
		setSize(w,h);
		}

	@Override
	public double getCanvasWidth() {
		return getWidth();
		}

	@Override
	public double getCanvasHeight() {
		return getHeight();
		}

	@Override
	public int getBackgroundRGB() {
		return getBackground().getRGB();
		}

	@Override
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		mGenericToolbar.paintContent(new SwingDrawContext((Graphics2D)g));
		}

	@Override
	public GenericDrawContext getDrawContext() {
		return new SwingDrawContext((Graphics2D)getGraphics());
		}
	}
