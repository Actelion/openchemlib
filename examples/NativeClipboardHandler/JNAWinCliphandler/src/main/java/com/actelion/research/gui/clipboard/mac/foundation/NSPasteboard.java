package com.actelion.research.gui.clipboard.mac.foundation;

import com.sun.jna.NativeLong;
import com.sun.jna.Pointer;

public class NSPasteboard extends NSObject {
    // The equivalent of NSPasteboardTypeString
    public static final NSString TypeString = new NSString("public.utf8-plain-text");

    // The equivalent of NSPasteboardTypePNG
    public static final NSString TypePNG = new NSString("public.png");

    // NSPasteboard -> https://developer.apple.com/documentation/foundation/nspasteboard?language=objc
    private static final Pointer nativeClass = Foundation.INSTANCE.objc_getClass("NSPasteboard");

    // [NSPasteboard generalPasteboard]; -> https://developer.apple.com/documentation/appkit/nspasteboard/1530091-generalpasteboard?language=objc
    private static final Pointer generalPasteboardSelector = Foundation.INSTANCE.sel_registerName("generalPasteboard");

    // [NSPasteboard setString:forType:];
    private static final Pointer setStringForTypeSelector = Foundation.INSTANCE.sel_registerName("setString:forType:");

    // [NSPasteboard setData:forType:];
    private static final Pointer setDataForTypeSelector = Foundation.INSTANCE.sel_registerName("setData:forType:");

    // [NSPasteboard stringForType:];
    private static final Pointer stringForTypeSelector = Foundation.INSTANCE.sel_registerName("stringForType:");

    // [NSPasteboard stringForType:];
    private static final Pointer dataForTypeSelector = Foundation.INSTANCE.sel_registerName("dataForType:");

    // [NSPasteboard clearContents]; -> https://developer.apple.com/documentation/appkit/nspasteboard/1533599-clearcontents?language=objc
    private static final Pointer clearContentsSelector = Foundation.INSTANCE.sel_registerName("clearContents");

    public NSPasteboard(NativeLong id) {
        super(id);
    }

    public static NSPasteboard generalPasteboard() {
        return new NSPasteboard(Foundation.INSTANCE.objc_msgSend(nativeClass, generalPasteboardSelector));
    }

    // https://developer.apple.com/documentation/appkit/nspasteboard/1528225-setstring?language=objc
    public void setString(String string, NSString type) {
        NSString nativeString = new NSString(string);
        Foundation.INSTANCE.objc_msgSend(getId(), setStringForTypeSelector, nativeString.getId(), type.getId());
    }

    // https://developer.apple.com/documentation/appkit/nspasteboard/1531214-setdata?language=objc
    public void setData(NSData data, NSString type) {
        Foundation.INSTANCE.objc_msgSend(getId(), setDataForTypeSelector, data.getId(), type.getId());
    }

    public NSString getString(NSString type) {
        NativeLong nativeResult = Foundation.INSTANCE.objc_msgSend(getId(), stringForTypeSelector, type.getId());

        long result = nativeResult.longValue();
        return result != 0 ? new NSString(nativeResult) : null;
    }

    public NSData getDataForType(NSString type) {
        NativeLong nativeResult = Foundation.INSTANCE.objc_msgSend(getId(), dataForTypeSelector, type.getId());

        long result = nativeResult.longValue();
        return result != 0 ? new NSData(nativeResult) : null;
    }

    public void clearContents() {
        Foundation.INSTANCE.objc_msgSend(getId(), clearContentsSelector);
    }
}