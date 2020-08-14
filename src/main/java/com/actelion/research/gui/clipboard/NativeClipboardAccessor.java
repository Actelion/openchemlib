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

package com.actelion.research.gui.clipboard;

public class NativeClipboardAccessor {
	public static final String NC_SKETCH	= "MDLSK";
	public static final String NC_CTAB		= "MDLCT";
	public static final String NC_MOLFILE	= "MDL_MOL";
	public static final String NC_METAFILE	= "CF_METAFILEPICT";
	public static final String NC_DIB		= "CF_DIB";
	public static final String NC_BITMAP	= "CF_BITMAP";
	public static final String NC_SERIALIZEMOLECULE = "ACT_MOLECULE";
	public static final String NC_SERIALIZEREACTION = "ACT_REACTION";
	public static final String NC_ALDUSMETAFILE	= "ALDUS_METAFILE";;
	public static final String NC_EMBEDDEDSKETCH = "MDLSK_EMBEDDED";
	public static final String NC_CHEMDRAWINTERCHANGE = "ChemDraw Interchange Format";
	public static final String NC_IDCODE		= "IDCODE";

    public static native boolean copyMoleculeToClipboard(String filname, byte[] sketch, byte[] serializedObject);
	public static native boolean copyReactionToClipboard(byte[] ctab, byte[] sketch, byte[] serializedObject);
	public static native boolean copySizedMoleculeToClipboard(String filname, byte[] sketch, byte[] serializedObject, int cx,int cy);
	/* Formats are "MDLSK","MDLCT","MDL_MOL","CF_METAFILEPICT","CF_DIB","ACT_MOLECULE","ACT_REACTION" */
    public static native byte[] getClipboardData(String format);
    public static native boolean setClipBoardData(String format, byte[] buffer);

    static {
        try {
			System.loadLibrary("actelionclip");
            System.out.println("actelionclip loaded");
        } catch (UnsatisfiedLinkError e) {
        	// added to retain compatibility with DataWarrior installations; TLS 11Jan2018
			e.printStackTrace();
        } catch (SecurityException e) {
        	e.printStackTrace();
        }
    }
}
