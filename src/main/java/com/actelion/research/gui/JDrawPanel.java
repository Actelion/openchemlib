/*
* Copyright (c) 1997 - 2016
* Actelion Pharmaceuticals Ltd.
* Gewerbestrasse 16
* CH-4123 Allschwil, Switzerland
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 3. Neither the name of the the copyright holder nor the
*    names of its contributors may be used to endorse or promote products
*    derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
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

	public JDrawPanel(StereoMolecule mol) {
		this(mol, 0);
	}

	public JDrawPanel(StereoMolecule mol, boolean isReaction) {
		this(mol, JDrawArea.MODE_MULTIPLE_FRAGMENTS
				| JDrawArea.MODE_REACTION);
	}

	public JDrawPanel(StereoMolecule mol, int mode) {
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
