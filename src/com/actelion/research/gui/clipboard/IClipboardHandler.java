/*
* @(#)IClipboardHandler.java
*
* Copyright 2005 Actelion Ltd. All Rights Reserved.
*
* This software is the proprietary information of Actelion Ltd.
* Use is subject to license terms.
*
 */
package com.actelion.research.gui.clipboard;

import com.actelion.research.chem.*;
import com.actelion.research.chem.reaction.Reaction;

/**
 *
 * <p>Title: Actelion Library</p>
 * <p>Description: Actelion Java Library</p>
 * <p>Copyright: Copyright (c) 2002-2003</p>
 * <p>Company: Actelion Ltd</p>
 * @author Thomas Sander, Christian Rufener
 * @version 1.0
 */
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
