package com.actelion.research.gui.dock;

import java.awt.event.MouseEvent;
import javax.swing.JPopupMenu;
import javax.swing.event.MouseInputAdapter;

public class HeaderMouseAdapter extends MouseInputAdapter {
    private PopupProvider mPopupProvider = null;
    private DockableHeader mHeader;

    public HeaderMouseAdapter(DockableHeader header) {
        mHeader = header;
        }

    public PopupProvider getPopupProvider() {
        return mPopupProvider;
        }

    public void setPopupProvider(PopupProvider p) {
        mPopupProvider = p;
        }

    public void mousePressed(MouseEvent e) {
        handlePopupTrigger(e);
        }

    public void mouseReleased(MouseEvent e) {
        handlePopupTrigger(e);
        }

    private void handlePopupTrigger(MouseEvent e) {
        if (mPopupProvider != null && e.isPopupTrigger()) {
            JPopupMenu popup = mPopupProvider.createPopupMenu(mHeader.getTitle(), mHeader.getDockable().isMaximized());
            popup.show(mHeader, e.getX(), e.getY());
            }
        }
    }

