package com.actelion.research.gui.swing;

import com.actelion.research.gui.generic.GenericKeyEvent;
import com.actelion.research.gui.generic.GenericKeyListener;

import java.awt.*;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.util.ArrayList;

public class SwingKeyHandler implements KeyListener {
	private ArrayList<GenericKeyListener> mListeners;

	public SwingKeyHandler() {
		mListeners = new ArrayList<>();
	}

	public void addListener(GenericKeyListener gkl) {
		mListeners.add(gkl);
	}

	public void removeListener(GenericKeyListener gkl) {
		mListeners.remove(gkl);
	}

	@Override
	public void keyTyped(KeyEvent e) {
		fireKeyEvent(createEvent(GenericKeyEvent.KEY_TYPED, e));
		}

	@Override
	public void keyPressed(KeyEvent e) {
		fireKeyEvent(createEvent(GenericKeyEvent.KEY_PRESSED, e));
		}

	@Override
	public void keyReleased(KeyEvent e) {
		fireKeyEvent(createEvent(GenericKeyEvent.KEY_RELEASED, e));
		}

	private void fireKeyEvent(GenericKeyEvent gke) {
		for (GenericKeyListener gkl:mListeners)
			gkl.keyActionHappened(gke);
	}

	private GenericKeyEvent createEvent(int what, KeyEvent e) {
		int button = (e.getModifiers() & InputEvent.BUTTON1_MASK) != 0 ? 1
				: (e.getModifiers() & InputEvent.BUTTON2_MASK) != 0 ? 2 : 3;
		return new GenericKeyEvent(what, e.getKeyCode(), e.getKeyChar(),
				(e.getModifiers() & Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()) != 0);
	}
}
