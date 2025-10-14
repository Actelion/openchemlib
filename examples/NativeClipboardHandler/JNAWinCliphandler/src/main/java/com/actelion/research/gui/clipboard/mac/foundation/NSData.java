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

public class NSData extends NSObject.Releasable
{
    // NSData -> https://developer.apple.com/documentation/foundation/NSData?language=objc
    private static final Pointer nativeClass = Foundation.INSTANCE.objc_getClass("NSData");

    // [NSData dataWithBytes:length:]; -> https://developer.apple.com/documentation/foundation/nsdata/datawithbytes?language=objc
    private static final Pointer dataWithBytesSelector = Foundation.INSTANCE.sel_registerName("dataWithBytes:length:");

    // [NSData data]; -> https://developer.apple.com/documentation/foundation/nsdata/data?language=objc
    private static final Pointer dataSelector = Foundation.INSTANCE.sel_registerName("data");

    private static final Pointer bytesSelector = Foundation.INSTANCE.sel_registerName("bytes");

    // [NSData length]; -> https://developer.apple.com/documentation/foundation/nsdata/length?language=objc
    private static final Pointer lengthSelector = Foundation.INSTANCE.sel_registerName("length");

    public NSData(NativeLong id)
    {
        super(id);
    }

    /**
     * Creates a new NSData object from the given byte array.
     *
     * @param bytes The byte array to create the NSData object from.
     * @return if bytes is null or 0, return an empty NSData object, otherwise return a new NSData object.
     */
    public static NSData initWithBytes(byte[] bytes)
    {
        if (bytes == null || bytes.length == 0)
            return new NSData(Foundation.INSTANCE.objc_msgSend(nativeClass, dataSelector));

        return new NSData(Foundation.INSTANCE.objc_msgSend(nativeClass, dataWithBytesSelector, bytes, new NativeLong(bytes.length, true)));
    }

    public int length()
    {
        return Foundation.INSTANCE.objc_msgSend(getId(), lengthSelector).intValue();
    }

    public NativeLong getData()
    {
        return Foundation.INSTANCE.objc_msgSend(getId(), bytesSelector);
    }

}