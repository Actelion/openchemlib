/*
 * Project: DD_jfx
 * @(#)LinuxNativeClipboardAccessor.java
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

import java.awt.*;
import java.awt.datatransfer.*;

public class LinuxNativeClipboardAccessor implements ClipboardOwner
{
	static public DataFlavor MOLFLAVOUR = new DataFlavor(com.actelion.research.chem.StereoMolecule.class,"Actelion Molecule");
	static public DataFlavor REACTIONFLAVOUR = new DataFlavor(com.actelion.research.chem.reaction.Reaction.class,"Actelion Reaction");


	public static boolean copyMoleculeToClipboard(String filname,byte[] sketch,byte[] serializedObject)
	{
		Toolkit.getDefaultToolkit().getSystemClipboard().setContents(new SerializedObjectSelection(serializedObject,MOLFLAVOUR),new LinuxNativeClipboardAccessor());
		return true;
	}

	public static boolean copyReactionToClipboard(String filname,byte[] sketch, byte[] serializedObject)
	{
		Toolkit.getDefaultToolkit().getSystemClipboard().setContents(new SerializedObjectSelection(serializedObject,REACTIONFLAVOUR),new LinuxNativeClipboardAccessor());
		return true;
	}

	public static byte[] getClipboardData(String format)
	{
        System.out.println("GetClipboardData " + format);
		try {
			Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
			DataFlavor df[] = {MOLFLAVOUR,REACTIONFLAVOUR};
			for (int i = 0; i < df.length; i++) {
				try{
					Object o = t.getTransferData(df[i]);
                    System.out.println("GetClipboardData Data: " + o);
					if(o instanceof byte[]){
						return(byte[])o;
					}
				} catch(Exception e){
                    System.out.println("Exception in getClipboardData: " + e);
				}
			}
		} catch (Exception e) {
			System.err.println("error getting clipboard data "+ e);

		}
		return null;
	}

	public static boolean setClipBoardData(String format,byte[] buffer)
	{
		throw new RuntimeException("NOT IMPLEMENTED YET");
	}


	public void lostOwnership(Clipboard parClipboard,Transferable parTransferable)
	{

		System.out.println("Lost ownership");
	}

	/**
	 * <p>Copyright: Copyright (c) 2003-2005</p>
	 * <p>Company: Actelion Ltd.</p>
	 *
	 * @author Christian Rufener
	 * @version 1.0
	 */
	static class SerializedObjectSelection implements Transferable
	{

		private DataFlavor flavor_ = null;
		private Object mol_;

		public SerializedObjectSelection(Object mol,DataFlavor f)
		{
			mol_ = mol;
			flavor_ = f;
		}

		public synchronized DataFlavor[] getTransferDataFlavors()
		{

			return new DataFlavor[] {flavor_} ;

		}

		public boolean isDataFlavorSupported(DataFlavor parFlavor)
		{

			return parFlavor.equals(flavor_);

		}

		public synchronized Object getTransferData(DataFlavor parFlavor) throws UnsupportedFlavorException
		{
			if(parFlavor.equals(flavor_)) {
				return mol_ ;
			}
			throw new UnsupportedFlavorException(flavor_);
		}

	}
}
