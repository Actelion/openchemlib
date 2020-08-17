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
import java.io.IOException;



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

		public void lostOwnership(Clipboard clipboard, Transferable contents){
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
	
	public static void main(String args[]) {
		Image img = pasteImage();
		System.out.println("Pasted image is " + img);
		copyImage(img);
		System.out.println("Copied image...");
	}
}

