package com.actelion.research.gui.dock;

import java.awt.*;
import java.awt.image.BufferedImage;

/**
 * @author Thomas Sander, idea taken from flexdock GhostPreview by Christopher Butler
 */
public class GhostPreview {
    private Rectangle mPreviewBounds;
    private Dockable mPreviewDockable;
    private BufferedImage mImage;

	public void createPreview(Dockable draggedDockable, Dockable targetDockable, int position, JDockingPanel port) {
        mPreviewBounds = getPreviewBounds(port.getAbsoluteBounds(targetDockable), position);
		if (mPreviewDockable != draggedDockable) {
			mImage = createImage(draggedDockable);
            mPreviewDockable = draggedDockable;
    		}
    	}
	
	public void drawPreview(Graphics2D g) {
		if (mImage == null)
			return;

        // create a solid preview outline
		g.setColor(Color.BLACK);
		g.drawRect(mPreviewBounds.x, mPreviewBounds.y, mPreviewBounds.width, mPreviewBounds.height);
		
		// make the graphics 20% translucent
		Composite original = g.getComposite();
		Composite composite = AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.8f);
		g.setComposite(composite);
		// now draw the preview image
		g.drawImage(mImage, mPreviewBounds.x, mPreviewBounds.y, mPreviewBounds.width, mPreviewBounds.height, null);
		g.setComposite(original);
		}
	
    private Rectangle getPreviewBounds(Rectangle targetBounds, int position) {
        switch (position) {
        case JDockingPanel.DOCK_TOP:
            targetBounds.height /= 2;
            return targetBounds;
        case JDockingPanel.DOCK_LEFT:
            targetBounds.width /= 2;
            return targetBounds;
        case JDockingPanel.DOCK_BOTTOM:
            targetBounds.y += targetBounds.height;
            targetBounds.height /= 2;
            targetBounds.y -= targetBounds.height;
            return targetBounds;
        case JDockingPanel.DOCK_RIGHT:
            targetBounds.x += targetBounds.width;
            targetBounds.width /= 2;
            targetBounds.x -= targetBounds.width;
            return targetBounds;
        default:
            return targetBounds;
            }
        }

    private BufferedImage createImage(Component comp) {
        BufferedImage image = (BufferedImage)comp.createImage(comp.getWidth(), comp.getHeight());
        Graphics g = image.createGraphics();
        comp.paintAll(g);
        return image;
    }
}
