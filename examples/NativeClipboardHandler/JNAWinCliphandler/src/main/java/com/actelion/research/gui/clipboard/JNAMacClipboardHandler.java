/*
 * Copyright (c) 2025,
 * Christian Rufener
 * Alipheron AG
 * www.alipheron.com
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

import com.actelion.research.gui.clipboard.mac.foundation.*;
import com.sun.jna.Pointer;

public class JNAMacClipboardHandler
{
    private static boolean isInitOK = true; // dummy to match NativeClipboardAccessor
    private static NSPasteboard pasteboard = NSPasteboard.generalPasteboard();

    public static boolean isInitOK()
    {
        return isInitOK;
    }

    public static byte[] getClipboardData(String format)
    {
        byte[] data = null;
        NSData dataForType = pasteboard.getDataForType(new NSString(format));
        if (dataForType != null) {
            int size = dataForType.length();
            data = new Pointer(dataForType.getData().longValue()).getByteArray(0, size);
        }
        return data;
    }

/*
   static void test() {
        // Example Codde
        Pointer classId = Foundation.INSTANCE.objc_getClass("NSPasteboard");
        Pointer clearContents = Foundation.INSTANCE.sel_registerName("clearContents");
        Pointer selector = Foundation.INSTANCE.sel_registerName("generalPasteboard");
        System.out.printf("class id %s Selector %s\n", classId, selector);
        NativeLong pasteBoard = Foundation.INSTANCE.objc_msgSend(classId, selector);
        System.out.printf("Pointer %s \n", pasteBoard);
        NativeLong ptr = Foundation.INSTANCE.objc_msgSend(pasteBoard, clearContents);
        System.out.printf("clearContent returned %s\n", ptr.longValue());
        Pointer name = Foundation.INSTANCE.sel_registerName("name");
        NativeLong namePtr = Foundation.INSTANCE.objc_msgSend(pasteBoard, name);
        //new ID(namePtr);
        System.out.printf("name returned '%s'\n", toStringViaUTF8(namePtr));

        try {
            System.out.printf("PasteBoard %s\n",pasteboard.getId());
            //pasteboard.clearContents();
            NSString ret = pasteboard.getString(new NSString("public.rtf"));
            NSData dataForType = pasteboard.getDataForType(new NSString("com.mdli.sketchfile"));
            if (dataForType != null) {
                int size = dataForType.length();
                byte[] byteArray = new Pointer(dataForType.getData().longValue()).getByteArray(0, size);

                System.out.printf("Data %s\n", new String(byteArray));
                StereoMolecule mol = new StereoMolecule();
                Sketch.createMolFromSketchBuffer(mol, byteArray);
                mol.ensureHelperArrays(Molecule.cHelperAll);
                System.out.printf("Mol Atoms %d Bonds %d MW %s\n", mol.getAllAtoms(), mol.getBonds(), mol.getMolweight());
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        Pointer pasteboardType = Foundation.INSTANCE.sel_registerName("PasteboardType");
        Pointer init = Foundation.INSTANCE.sel_registerName("init");
        NativeLong foo = Foundation.INSTANCE.objc_msgSend(pasteBoard, pasteboardType
                //,init
        );
        System.out.printf("PasteboardType returned %s\n", foo);
    }

    static String toStringViaUTF8(NativeLong cfString)
    {
        //if (ID.NIL.equals(cfString)) return null;

        int lengthInChars = Foundation.INSTANCE.CFStringGetLength(cfString);
        int potentialLengthInBytes = 3 * lengthInChars + 1; // UTF8 fully escaped 16 bit chars, plus nul

        byte[] buffer = new byte[potentialLengthInBytes];
        byte ok = Foundation.INSTANCE.CFStringGetCString(cfString, buffer, buffer.length, Foundation.INSTANCE.kCFStringEncodingUTF8);
        if (ok == 0) throw new RuntimeException("Could not convert string");
        return Native.toString(buffer);
    }
*/
}
