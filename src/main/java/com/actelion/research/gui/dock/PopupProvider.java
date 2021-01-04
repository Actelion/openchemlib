package com.actelion.research.gui.dock;

import javax.swing.*;

public interface PopupProvider {
    public JPopupMenu createPopupMenu(String title, boolean isMaximized);
    }
