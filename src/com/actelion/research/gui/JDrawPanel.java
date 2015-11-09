/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
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

import java.awt.*;
import javax.swing.*;
import com.actelion.research.chem.*;
import com.actelion.research.chem.reaction.IReactionMapper;

public class JDrawPanel extends JPanel {
    static final long serialVersionUID = 0x20061019;

    protected JDrawToolbar mToolBar;
    
	protected JDrawArea mArea;

	public JDrawPanel(Frame parent, StereoMolecule mol) {
		this(parent, mol, 0);
	}

	public JDrawPanel(Frame parent, StereoMolecule mol, boolean isReaction) {
		this(parent, mol, JDrawArea.MODE_MULTIPLE_FRAGMENTS
				| JDrawArea.MODE_REACTION);
	}

	public JDrawPanel(Frame parent, StereoMolecule mol, int mode) {
		setLayout(new BorderLayout());

		mArea = new JDrawArea(mol, mode);
		add(mArea, BorderLayout.CENTER);

		mToolBar = new JDrawToolbar(mArea, mode);
		add(mToolBar, BorderLayout.WEST);
	}

    public void setMapper(IReactionMapper mapper)
    {
        mArea.setMapper(mapper);
    }
	public JDrawArea getDrawArea() {
		return mArea;
	}
	
	public void cleanStructure(){
		mArea.toolChanged(1);
	}
}


