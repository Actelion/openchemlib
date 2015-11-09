/*
 * Project: DD_jfx
 * @(#)IClipboardHandler.java
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
package com.actelion.research.gui.clipboard;

import com.actelion.research.chem.*;
import com.actelion.research.chem.reaction.Reaction;

public interface IClipboardHandler {
	public StereoMolecule pasteMolecule();
	public Reaction pasteReaction();
	public boolean copyMolecule(String molfile);
	public boolean copyMolecule(StereoMolecule mol);
	public boolean copyReaction(Reaction r);
	public boolean copyReaction(String ctab);
	public boolean copyImage(java.awt.Image img);
	public java.awt.Image pasteImage();
	}
