package com.actelion.research.gui.clipboard;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.clipboard.mac.foundation.*;
import com.actelion.research.util.Sketch;
import com.sun.jna.*;
import com.sun.jna.platform.mac.CoreFoundation;
import com.sun.jna.platform.win32.User32;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Map;

/**
 * Copyright (c) 2009-2025
 * Christian Rufener async.ch
 * All rights reserved
 * User: christian
 * Date: 12.10.25
 * Time: 19:17
 */
// The interface of the runtime library of Objective-C.
/*
interface Foundation extends Library
{
    int kCFStringEncodingUTF8 = 0x08000100;
    Map<String, Integer> loadoptions = new HashMap()
    {
        {
            put(Library.OPTION_STRING_ENCODING, StandardCharsets.UTF_8.name());
        }
    };

    Foundation INSTANCE = Native.load("Foundation", Foundation.class, loadoptions);

    // https://developer.apple.com/documentation/objectivec/1418952-objc_getclass?language=objc
    Pointer objc_getClass(String className);

    // https://developer.apple.com/documentation/objectivec/1418557-sel_registername?language=objc
    Pointer sel_registerName(String selectorName);

    // https://developer.apple.com/documentation/objectivec/1456712-objc_msgsend?language=objc
    // The return type is actually "generic". You might need to declare this function
    // multiple times with different return types if you need them.
    Pointer objc_msgSend(Pointer receiver, Pointer selector, Object... args);

    int CFStringGetLength(Pointer theString);

    byte CFStringGetCString(Pointer theString, byte[] buffer, int bufferSize, int encoding);

}
*/

final class ID extends NativeLong
{

    public ID()
    {
    }

    public ID(long peer)
    {
        super(peer);
    }

    public static final ID NIL = new ID(0L);

    public boolean booleanValue()
    {
        return intValue() != 0;
    }
}

public class JNAMacClipboardHandler
{

    public static void foo() throws IOException
    {
/*
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
*/
        NSPasteboard pasteboard = NSPasteboard.generalPasteboard();
        System.out.printf("PasteBoard %s\n",pasteboard.getId());
        //pasteboard.clearContents();
        NSString ret = pasteboard.getString(new NSString("public.rtf"));
//        System.out.printf("Data %s\n",ret);
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

        /*
        Pointer pasteboardType = Foundation.INSTANCE.sel_registerName("PasteboardType");
        Pointer init = Foundation.INSTANCE.sel_registerName("init");
        NativeLong foo = Foundation.INSTANCE.objc_msgSend(pasteBoard, pasteboardType
                //,init
        );
//                CoreFoundation.CFStringRef .createCFString("ch.async.cr"));
        System.out.printf("PasteboardType returned %s\n", foo);
*/

        //let NSFilenamesPboardTypeTemp = NSPasteboard.PasteboardType("NSFilenamesPboardType")
    }

    public static String toStringViaUTF8(NativeLong cfString)
    {
        //if (ID.NIL.equals(cfString)) return null;

        int lengthInChars = Foundation.INSTANCE.CFStringGetLength(cfString);
        int potentialLengthInBytes = 3 * lengthInChars + 1; // UTF8 fully escaped 16 bit chars, plus nul

        byte[] buffer = new byte[potentialLengthInBytes];
        byte ok = Foundation.INSTANCE.CFStringGetCString(cfString, buffer, buffer.length, Foundation.INSTANCE.kCFStringEncodingUTF8);
        if (ok == 0) throw new RuntimeException("Could not convert string");
        return Native.toString(buffer);
    }

}
