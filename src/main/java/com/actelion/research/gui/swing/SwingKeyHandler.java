package com.actelion.research.gui.swing;

import com.actelion.research.gui.generic.GenericEventHandler;
import com.actelion.research.gui.generic.GenericKeyEvent;
import javafx.scene.input.KeyCode;

import java.awt.*;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

public class SwingKeyHandler extends GenericEventHandler<GenericKeyEvent> implements KeyListener {

	public SwingKeyHandler(Object source) {
		super(source);
	}

	@Override
	public void keyTyped(KeyEvent e) {
//		fireKeyEvent(createEvent(GenericKeyEvent.KEY_TYPED, e));
		}

	@Override
	public void keyPressed(KeyEvent e) {
		fireEvent(createEvent(GenericKeyEvent.KEY_PRESSED, e));
		}

	@Override
	public void keyReleased(KeyEvent e) {
		fireEvent(createEvent(GenericKeyEvent.KEY_RELEASED, e));
		}

	private GenericKeyEvent createEvent(int what, KeyEvent e) {
		int key = e.getKeyCode() == KeyEvent.VK_ALT ? GenericKeyEvent.KEY_ALT
				: e.getKeyCode() == KeyEvent.VK_CONTROL ? GenericKeyEvent.KEY_CTRL
				: e.getKeyCode() == KeyEvent.VK_SHIFT ? GenericKeyEvent.KEY_SHIFT
				: e.getKeyCode() == KeyEvent.VK_DELETE ? GenericKeyEvent.KEY_DELETE
				: e.getKeyCode() == KeyEvent.VK_BACK_SPACE ? GenericKeyEvent.KEY_BACK_SPACE
				: e.getKeyCode() == KeyEvent.VK_HELP ? GenericKeyEvent.KEY_HELP
				: e.getKeyCode() == KeyEvent.VK_ESCAPE ? GenericKeyEvent.KEY_ESCAPE
				: e.getKeyCode() == KeyEvent.VK_ENTER ? GenericKeyEvent.KEY_ENTER
				: e.getKeyCode() == KeyEvent.VK_LEFT ? GenericKeyEvent.KEY_LEFT
				: e.getKeyCode() == KeyEvent.VK_RIGHT ? GenericKeyEvent.KEY_RIGHT
				: e.getKeyCode() == KeyEvent.VK_UP ? GenericKeyEvent.KEY_UP
				: e.getKeyCode() == KeyEvent.VK_DOWN ? GenericKeyEvent.KEY_DOWN
				: e.getKeyChar();
		boolean isCtrlDown = (e.getModifiers() & Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()) != 0;
		boolean isAltDown = (e.getModifiers() & Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()) != 0;
		boolean isShiftDown = (e.getModifiers() & Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()) != 0;
		boolean isShortcut = (e.getModifiers() & Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()) != 0;
		// BH needs to check isCtrlDown here otherwise VK_ENTER becomes 'j'
		// BH the adjustment is only applicable if e.getKeyChar() != e.getKeyCode() 
		if (isCtrlDown 
				&& key != e.getKeyCode()
				&& key >= 1 && key <= 26)  // strangely, if Ctrl is pressed, letters are encoded from 1-26
			key = 'a' + key - 1;
		return new GenericKeyEvent(what, key, isAltDown, isCtrlDown, isShiftDown, isShortcut, getSource());
	}
}
