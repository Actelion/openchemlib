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
//			DataFlavor[] flavors = t.getTransferDataFlavors();
//			for (DataFlavor df: flavors) {
//				System.out.printf("Flavor: %s\n",df);
//				Object o = t.getTransferData(df);
//				System.out.printf("Returned object is %s\n",o);
//			}

			if (format.equals(NativeClipboardHandler.NC_IDCODE)) {
				Object o = t.getTransferData(DataFlavor.stringFlavor);
				if (o != null) {
					return o.toString().getBytes();
				}
			}
			DataFlavor df = format.equalsIgnoreCase(NativeClipboardHandler.NC_SERIALIZEREACTION) ? REACTIONFLAVOUR : MOLFLAVOUR;
			try{
				Object o = t.getTransferData(df);
				System.out.println("GetClipboardData Data: " + o);
				if(o instanceof byte[]){
					return(byte[])o;
				}
			} catch(Exception e){
				System.out.println("Exception in getClipboardData: " + e);
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
