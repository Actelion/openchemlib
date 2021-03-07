/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.gui;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.hidpi.HiDPIHelper;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;

public class JEditableStructureView extends JStructureView {
    static final long serialVersionUID = 0x20090727;

    private static final String EDIT_MESSAGE = "<double click or drag & drop>";
    private boolean mAllowQueryFeatures;

    public JEditableStructureView() {
        super(null);
	    setEditable(true);
		mAllowQueryFeatures = true;
		}

	public JEditableStructureView(StereoMolecule mol) {
        super(mol);
		setEditable(true);
		mAllowQueryFeatures = true;
	    }

	public JEditableStructureView(int dragAction, int dropAction) {
        super(null, dragAction, dropAction);
		setEditable(true);
		mAllowQueryFeatures = true;
	    }

	public JEditableStructureView(StereoMolecule mol, int dragAction, int dropAction) {
        super(mol, dragAction, dropAction);
		setEditable(true);
		mAllowQueryFeatures = true;
	    }

	@Override
	public void paintComponent(Graphics g) {
        Dimension theSize = getSize();
        Insets insets = getInsets();
        theSize.width -= insets.left + insets.right;
        theSize.height -= insets.top + insets.bottom;

		super.paintComponent(g);

		if (isEnabled() && isEditable() && getMolecule().getAllAtoms() == 0) {
	        g.setFont(g.getFont().deriveFont(Font.PLAIN, HiDPIHelper.scale(10)));
	        FontMetrics metrics = g.getFontMetrics();
	        Rectangle2D bounds = metrics.getStringBounds(EDIT_MESSAGE, g);
	        g.drawString(EDIT_MESSAGE, (int)(insets.left+theSize.width-bounds.getWidth())/2,
	                                   (insets.top+theSize.height-metrics.getHeight())/2+metrics.getAscent());
	        }
	    }

    public void mouseClicked(MouseEvent e) {
        if (e.getClickCount() == 2 && isEnabled() && isEditable()) {
            JDrawDialog theDialog = createDrawDialog();
            theDialog.getDrawArea().setAllowQueryFeatures(mAllowQueryFeatures);
            theDialog.addStructureListener(this);
            theDialog.setVisible(true);
            }
        }

    protected JDrawDialog createDrawDialog() {
		Component c = this;
		while (!(c instanceof Frame || c instanceof Dialog))
			c = c.getParent();
		return (c instanceof Frame) ? new JDrawDialog((Frame) c, getMolecule(), Dialog.ModalityType.DOCUMENT_MODAL) : new JDrawDialog((Dialog) c, getMolecule(), Dialog.ModalityType.DOCUMENT_MODAL);
		}

	public void setAllowQueryFeatures(boolean allow) {
		if (mAllowQueryFeatures != allow) {
			mAllowQueryFeatures = allow;
			if (!allow && getMolecule().removeQueryFeatures())
				structureChanged();
			}
		}
    }
