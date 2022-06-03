package com.actelion.research.gui.dock;

import com.actelion.research.gui.swing.SwingCursorHelper;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import java.awt.event.MouseEvent;

public class HeaderMouseAdapter extends MouseInputAdapter {
    private PopupProvider mPopupProvider = null;
    private JLabel mTitleLabel;
    private Dockable mDockable;
    private boolean mIsMouseDown;

    public HeaderMouseAdapter(JLabel titleLabel, Dockable dockable) {
        mTitleLabel = titleLabel;
        mDockable = dockable;
        }

    public PopupProvider getPopupProvider() {
        return mPopupProvider;
        }

    public void setPopupProvider(PopupProvider p) {
        mPopupProvider = p;
        }

    @Override
    public void mouseEntered(MouseEvent e) {
        mTitleLabel.setCursor(SwingCursorHelper.getCursor(SwingCursorHelper.cHandCursor));
        mIsMouseDown = false;
        }

    @Override
    public void mouseExited(MouseEvent e) {
        if (!mIsMouseDown)
            mTitleLabel.setCursor(SwingCursorHelper.getCursor(SwingCursorHelper.cPointerCursor));
        mIsMouseDown = false;
        }

    @Override
    public void mousePressed(MouseEvent e) {
        if (!handlePopupTrigger(e)) {
            mTitleLabel.setCursor(SwingCursorHelper.getCursor(e.getButton() == MouseEvent.BUTTON1 ? SwingCursorHelper.cFistCursor : SwingCursorHelper.cPointerCursor));
            mIsMouseDown = true;
            mDockable.getDockingPanel().selectDockable(mDockable);
            }
        }

    @Override
    public void mouseReleased(MouseEvent e) {
        handlePopupTrigger(e);
        mTitleLabel.setCursor(SwingCursorHelper.getCursor(SwingCursorHelper.cHandCursor));
        mIsMouseDown = false;
        }

    private boolean handlePopupTrigger(MouseEvent e) {
        if (mPopupProvider != null && e.isPopupTrigger()) {
            mTitleLabel.setCursor(SwingCursorHelper.getCursor(SwingCursorHelper.cPointerCursor));
            JPopupMenu popup = mPopupProvider.createPopupMenu(mTitleLabel.getText(), mDockable.isMaximized());
            popup.show(mTitleLabel, e.getX(), e.getY());
            return true;
            }

        return false;
        }
    }
