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

import javax.swing.*;
import java.io.File;

public class JFileChooserOverwrite extends JFileChooser {
	private static final long serialVersionUID = 20150101L;

	private File	mFile;
	private String[] mExtensions = null;

	public JFileChooserOverwrite() {
//		super(System.getProperty("user.dir"));
		}

	public void setExtensions(String[] extensions) {
		mExtensions = extensions;
		}

	public File getFile() {
		return mFile;
		}

	public void approveSelection() {
		if (getSelectedFile() != null) {
			String filename = getSelectedFile().getPath();
			if (mExtensions != null) {
				int dotIndex = filename.lastIndexOf('.');
				int slashIndex = filename.lastIndexOf(File.separator);
				if (dotIndex == -1
				 || dotIndex < slashIndex) {
					filename = filename.concat(mExtensions[0]);
					}
				else {
					boolean found = false;
					for (String extension:mExtensions)
				        if (filename.substring(dotIndex).equalsIgnoreCase(extension))
				        	found = true;

				    if (!found) {
						JOptionPane.showMessageDialog(this, "Incompatible file name extension.");
					    return;
						}
					}
				}

			mFile = new File(filename);
			if (mFile != null) {
				if (mFile.exists()) {
					int answer = JOptionPane.showConfirmDialog(this,
						"This file already exists. Do you want to replace it?", "Warning",
						JOptionPane.OK_CANCEL_OPTION, JOptionPane.WARNING_MESSAGE);

					if (answer == JOptionPane.OK_OPTION)
						super.approveSelection();
					}
				else
		        	super.approveSelection();
				}
			}
		}
	}
