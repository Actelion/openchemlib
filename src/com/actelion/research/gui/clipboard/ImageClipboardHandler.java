/*
 * Project: DD_jfx
 * @(#)ImageClipboardHandler.java
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

import java.awt.geom.Rectangle2D;

import java.io.*;

import com.actelion.research.chem.*;

import com.actelion.research.util.*;



public class ImageClipboardHandler {

    public static Image pasteImage() {

        Clipboard systemClipboard = Toolkit.getDefaultToolkit().getSystemClipboard();

        Transferable clipboardContents = systemClipboard.getContents(null);

  

        if (clipboardContents == null)

            return null;

        else {

            try {

                if (clipboardContents.isDataFlavorSupported(DataFlavor.imageFlavor)) {

                    Image image = (Image)clipboardContents.getTransferData(DataFlavor.imageFlavor);

                    return image;

               	}

            } catch (UnsupportedFlavorException ufe) {

                ufe.printStackTrace();

            } catch (IOException ioe) {

                ioe.printStackTrace();

            }

            return null;

        }

    }



    public static boolean copyImage(Image image) {

        ImageSelection imageSelection = new ImageSelection(image);

        Toolkit.getDefaultToolkit().getSystemClipboard().setContents(imageSelection, imageSelection);

        return true;

    }



    public static class ImageSelection implements Transferable,ClipboardOwner {


	    // the Image object which will be housed by the ImageSelection

        private Image mImage;



        public ImageSelection(Image image) {

            mImage = image;

        }

		public void lostOwnership(Clipboard clipboard, Transferable contents)
		{
			System.out.println("Lost ownership...");	
		}

        public DataFlavor[] getTransferDataFlavors() {

            return new DataFlavor[] { DataFlavor.imageFlavor };

        }

    

        public boolean isDataFlavorSupported(DataFlavor flavor) {

            return DataFlavor.imageFlavor.equals(flavor);

        }



        public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException,IOException {

            if (!DataFlavor.imageFlavor.equals(flavor)) 

                throw new UnsupportedFlavorException(flavor);



            return mImage;

        }

    }
	
	public static void main(String args[])
	{
		Image img = pasteImage();
		System.out.println("PAsted image is " + img);
		copyImage(img);
		System.out.println("Copied image...");
	}
}

