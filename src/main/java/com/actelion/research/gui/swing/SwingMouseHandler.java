package com.actelion.research.gui.swing;

import com.actelion.research.gui.generic.GenericMouseEvent;
import com.actelion.research.gui.generic.GenericMouseListener;

import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;

public class SwingMouseHandler implements MouseListener,MouseMotionListener {
	private ArrayList<GenericMouseListener> mListeners;

	public SwingMouseHandler() {
		mListeners = new ArrayList<>();
		}

	public void addListener(GenericMouseListener gml) {
		mListeners.add(gml);
		}

	public void removeListener(GenericMouseListener gml) {
		mListeners.remove(gml);
		}

	@Override
	public void mousePressed(MouseEvent e) {
		fireMouseEvent(createEvent(GenericMouseEvent.MOUSE_PRESSED, e));
	}

	@Override
	public void mouseReleased(MouseEvent e) {
		fireMouseEvent(createEvent(GenericMouseEvent.MOUSE_RELEASED, e));
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		fireMouseEvent(createEvent(GenericMouseEvent.MOUSE_CLICKED, e));
	}

	@Override
	public void mouseEntered(MouseEvent e) {
		fireMouseEvent(createEvent(GenericMouseEvent.MOUSE_ENTERED, e));
	}

	@Override
	public void mouseExited(MouseEvent e) {
		fireMouseEvent(createEvent(GenericMouseEvent.MOUSE_EXITED, e));
	}

	@Override
	public void mouseMoved(MouseEvent e) {
		fireMouseEvent(createEvent(GenericMouseEvent.MOUSE_MOVED, e));
	}

	@Override
	public void mouseDragged(MouseEvent e) {
		fireMouseEvent(createEvent(GenericMouseEvent.MOUSE_DRAGGED, e));
	}

	private void fireMouseEvent(GenericMouseEvent gme) {
		for (GenericMouseListener gml:mListeners)
			gml.mouseActionHappened(gme);
	}

	private GenericMouseEvent createEvent(int what, MouseEvent e) {
		int button = (e.getModifiers() & InputEvent.BUTTON1_MASK) != 0 ? 1
				   : (e.getModifiers() & InputEvent.BUTTON2_MASK) != 0 ? 2 : 3;
		return new GenericMouseEvent(what, button, e.getClickCount(), e.getX(), e.getY(),
				e.isShiftDown(), e.isControlDown(), e.isAltDown(), e.isPopupTrigger());
	}
}
