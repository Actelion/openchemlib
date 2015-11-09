/*
 * Project: DD_jfx
 * @(#)MoleculeDragAdapter.java
 *
 * Copyright (c) 1997- 2015
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All Rights Reserved.
 *
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 *
 * Author: Christian Rufener
 */

package com.actelion.research.gui.dnd;

import com.actelion.research.chem.AbstractDepictor;
import com.actelion.research.chem.Depictor2D;
import com.actelion.research.chem.StereoMolecule;

import java.awt.*;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.dnd.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.IOException;


public abstract class MoleculeDragAdapter implements DragSourceListener//,DragSourceMotionListener
{
    int mAllowedDragAction;
    Cursor cursor = null;

    public MoleculeDragAdapter(Component c)
    {
//        DragSource ds = new DragSource();
//        ds.addDragSourceListener(new WindowMovingDragSourceListener());

        DragSource.getDefaultDragSource().createDefaultDragGestureRecognizer(
            c, DnDConstants.ACTION_COPY_OR_MOVE, new DGListener());
        mAllowedDragAction = DnDConstants.ACTION_COPY_OR_MOVE;
    }

    public void allowDrag(boolean allow)
    {
        if (allow)
            mAllowedDragAction = DnDConstants.ACTION_COPY_OR_MOVE;
        else
            mAllowedDragAction = DnDConstants.ACTION_NONE;
    }

    public abstract Transferable getTransferable(Point origin);

    public void onDragEnter()
    {

    }

    public void onDragOver()
    {

    }

    public void onDragExit()
    {

    }

    public final void dragDropEnd(DragSourceDropEvent e)
    {
    }

    public final void dragEnter(DragSourceDragEvent e)
    {
        DragSourceContext context = e.getDragSourceContext();
        int dropAction = e.getDropAction();
        if ((dropAction & mAllowedDragAction) != 0) {
//            if (cursor != null)
//                context.setCursor(cursor);
//            else
                context.setCursor(DragSource.DefaultCopyDrop);
        } else {
            context.setCursor(DragSource.DefaultCopyNoDrop);
        }
        onDragEnter();
    }

    public final void dragOver(DragSourceDragEvent e)
    {
        onDragOver();
    }

    public final void dragExit(DragSourceEvent e)
    {
        onDragExit();
    }

    public final void dropActionChanged(DragSourceDragEvent e)
    {
        DragSourceContext context = e.getDragSourceContext();
        int dropAction = e.getDropAction();
        if ((dropAction & mAllowedDragAction) != 0) {
//            if (cursor != null)
//                context.setCursor(cursor);
//            else
                context.setCursor(DragSource.DefaultCopyDrop);
        } else {
            context.setCursor(DragSource.DefaultCopyNoDrop);
        }
    }


//    class WindowMovingDragSourceListener extends DragSourceAdapter implements DragSourceMotionListener {
//
//        private Window win;
//        private boolean firstCall = true;
//
//        public WindowMovingDragSourceListener(Window w) {
//            win = w;
//        }
//
//        public void dragDropEnd(DragSourceDropEvent dsde) {
//            System.err.println("dragDropEnd(): " + dsde);
//            win.setVisible(false);
//            firstCall = true;
//        }
//
//        public void dragMouseMoved(DragSourceDragEvent dsde) {
//            Point p = dsde.getLocation();
//            p.translate(32, 16);
//            win.setLocation(p);
//            if (firstCall) {
//                win.setVisible(true);
//                firstCall = false;
//            }
//        }
//    }

    class DGListener implements DragGestureListener
    {
        public void dragGestureRecognized(DragGestureEvent e)
        {
            onDragBegin(e);
        }

    }


    protected void onDragBegin(DragGestureEvent e)
    {
        cursor = null;
        if ((e.getDragAction() & mAllowedDragAction) != 0) {
            Transferable transferable = getTransferable(e.getDragOrigin());
            if (transferable != null) {
                try {
                    Image img = null;
                    if (DragSource.isDragImageSupported() && (img = drawDragImage(transferable,400,400)) != null) {
                        e.startDrag(DragSource.DefaultCopyDrop, img, new Point(200, 200), transferable, this);
                    } else {
                        e.startDrag(cursor, transferable, this);
                    }
                } catch (InvalidDnDOperationException idoe) {
                    System.err.println(idoe);
                }
            }
        }
    }

    public Image drawDragImage(Transferable transferable,int width,int height)
    {
        if (transferable instanceof MoleculeTransferable) {
            try {
                MoleculeTransferable t = (MoleculeTransferable) transferable;
                Object o = t.getTransferData(MoleculeFlavors.DF_SERIALIZEDOBJECT);
                if (o instanceof StereoMolecule) {
                    StereoMolecule mol = (StereoMolecule) o;
                    Depictor2D depict = new Depictor2D(mol);
                    BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
                    Graphics g = img.getGraphics();
                    ((Graphics2D)g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
    		        ((Graphics2D)g).setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);
                    depict.validateView(g, new Rectangle2D.Float(0, 0, width, height), AbstractDepictor.cModeInflateToMaxAVBL);
                    depict.paint(g);
                    return img;
                }
            } catch (IOException e1) {
                System.err.println(e1);
            } catch (UnsupportedFlavorException e1) {
                System.err.println(e1);
            }
        }
        return null;
    }


//    final Window w = new Window() {
//                public void paint(Graphics g) {
//                    g.drawString("image", 5, 20); // or render an image via g.drawImage()
//                }
//            };

}
