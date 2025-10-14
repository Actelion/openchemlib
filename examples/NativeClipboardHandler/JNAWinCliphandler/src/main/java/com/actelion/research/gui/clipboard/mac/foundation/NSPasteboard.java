/*
MIT License

Copyright (c) 2024 Caoimhe Byrne, 2025 Alipheron AG

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */
package com.actelion.research.gui.clipboard.mac.foundation;

import com.sun.jna.NativeLong;
import com.sun.jna.Pointer;

public class NSPasteboard extends NSObject
{
    // The equivalent of NSPasteboardTypeString
    public static final NSString TypeString = new NSString("public.utfplain-text");

    // The equivalent of NSPasteboardTypePNG
    public static final NSString TypePNG = new NSString("public.png");

    // NSPasteboard -> https://developer.apple.com/documentation/foundation/nspasteboard?language=objc
    private static final Pointer nativeClass = Foundation.INSTANCE.objc_getClass("NSPasteboard");

    // [NSPasteboard generalPasteboard]; -> https://developer.apple.com/documentation/appkit/nspasteboard/generalpasteboard?language=objc
    private static final Pointer generalPasteboardSelector = Foundation.INSTANCE.sel_registerName("generalPasteboard");

    // [NSPasteboard setString:forType:];
    private static final Pointer setStringForTypeSelector = Foundation.INSTANCE.sel_registerName("setString:forType:");

    // [NSPasteboard setData:forType:];
    private static final Pointer setDataForTypeSelector = Foundation.INSTANCE.sel_registerName("setData:forType:");

    // [NSPasteboard stringForType:];
    private static final Pointer stringForTypeSelector = Foundation.INSTANCE.sel_registerName("stringForType:");

    // [NSPasteboard stringForType:];
    private static final Pointer dataForTypeSelector = Foundation.INSTANCE.sel_registerName("dataForType:");

    // [NSPasteboard clearContents]; -> https://developer.apple.com/documentation/appkit/nspasteboard/clearcontents?language=objc
    private static final Pointer clearContentsSelector = Foundation.INSTANCE.sel_registerName("clearContents");

    public NSPasteboard(NativeLong id)
    {
        super(id);
    }

    public static NSPasteboard generalPasteboard()
    {
        return new NSPasteboard(Foundation.INSTANCE.objc_msgSend(nativeClass, generalPasteboardSelector));
    }

    // https://developer.apple.com/documentation/appkit/nspasteboard/setstring?language=objc
    public void setString(String string, NSString type)
    {
        NSString nativeString = new NSString(string);
        Foundation.INSTANCE.objc_msgSend(getId(), setStringForTypeSelector, nativeString.getId(), type.getId());
    }

    // https://developer.apple.com/documentation/appkit/nspasteboard/setdata?language=objc
    public void setData(NSData data, NSString type)
    {
        Foundation.INSTANCE.objc_msgSend(getId(), setDataForTypeSelector, data.getId(), type.getId());
    }

    public NSString getString(NSString type)
    {
        NativeLong nativeResult = Foundation.INSTANCE.objc_msgSend(getId(), stringForTypeSelector, type.getId());

        long result = nativeResult.longValue();
        return result != 0 ? new NSString(nativeResult) : null;
    }

    public NSData getDataForType(NSString type)
    {
        NativeLong nativeResult = Foundation.INSTANCE.objc_msgSend(getId(), dataForTypeSelector, type.getId());

        long result = nativeResult.longValue();
        return result != 0 ? new NSData(nativeResult) : null;
    }

    public void clearContents()
    {
        Foundation.INSTANCE.objc_msgSend(getId(), clearContentsSelector);
    }
}