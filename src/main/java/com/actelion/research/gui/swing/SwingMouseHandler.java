package com.actelion.research.gui.swing;

import com.actelion.research.gui.generic.GenericEventHandler;
import com.actelion.research.gui.generic.GenericMouseEvent;

import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

public class SwingMouseHandler extends GenericEventHandler<GenericMouseEvent> implements MouseListener,MouseMotionListener {

	public SwingMouseHandler(Object source) {
		super(source);
		}

	@Override
	public void mousePressed(MouseEvent e) {
		fireEvent(createEvent(GenericMouseEvent.MOUSE_PRESSED, e));
	}

	@Override
	public void mouseReleased(MouseEvent e) {
		fireEvent(createEvent(GenericMouseEvent.MOUSE_RELEASED, e));
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		fireEvent(createEvent(GenericMouseEvent.MOUSE_CLICKED, e));
	}

	@Override
	public void mouseEntered(MouseEvent e) {
		fireEvent(createEvent(GenericMouseEvent.MOUSE_ENTERED, e));
	}

	@Override
	public void mouseExited(MouseEvent e) {
		fireEvent(createEvent(GenericMouseEvent.MOUSE_EXITED, e));
	}

	@Override
	public void mouseMoved(MouseEvent e) {
		fireEvent(createEvent(GenericMouseEvent.MOUSE_MOVED, e));
	}

	@Override
	public void mouseDragged(MouseEvent e) {
		fireEvent(createEvent(GenericMouseEvent.MOUSE_DRAGGED, e));
	}

	private GenericMouseEvent createEvent(int what, MouseEvent e) {
		int button = (e.getModifiers() & InputEvent.BUTTON1_MASK) != 0 ? 1
				   : (e.getModifiers() & InputEvent.BUTTON2_MASK) != 0 ? 2 : 3;
		return new GenericMouseEvent(what, button, e.getClickCount(), e.getX(), e.getY(),
				e.isShiftDown(), e.isControlDown(), e.isAltDown(), e.isPopupTrigger(), getSource());
	}
}
