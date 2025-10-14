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

import com.sun.jna.*;

public interface Foundation extends Library
{
    int kCFStringEncodingUTF8 = 0x08000100;

    Foundation INSTANCE = Native.load("Foundation", Foundation.class);

    // void objc_msgSend(*); -> https://developer.apple.com/documentation/objectivec/1456712-objc_msgsend?language=objc
    NativeLong objc_msgSend(Pointer receiver, Pointer selector);

    // void objc_msgSend(*); -> https://developer.apple.com/documentation/objectivec/1456712-objc_msgsend?language=objc
    NativeLong objc_msgSend(NativeLong receiver, Pointer selector);

    // void objc_msgSend(*); -> https://developer.apple.com/documentation/objectivec/1456712-objc_msgsend?language=objc
    NativeLong objc_msgSend(NativeLong receiver, Pointer selector, NativeLong arg1, NativeLong arg2);

    // void objc_msgSend(*); -> https://developer.apple.com/documentation/objectivec/1456712-objc_msgsend?language=objc
    NativeLong objc_msgSend(Pointer receiver, Pointer selector, byte[] arg1, NativeLong arg2);

    // void objc_msgSend(*); -> https://developer.apple.com/documentation/objectivec/1456712-objc_msgsend?language=objc
    NativeLong objc_msgSend(Pointer receiver, Pointer selector, String arg1);

    // void objc_msgSend(*); -> https://developer.apple.com/documentation/objectivec/1456712-objc_msgsend?language=objc
    NativeLong objc_msgSend(NativeLong receiver, Pointer selector, NativeLong arg1);

    // void objc_msgSend(*); -> https://developer.apple.com/documentation/objectivec/1456712-objc_msgsend?language=objc
    NativeLong objc_msgSend(NativeLong receiver, Pointer selector, boolean arg1);

    // void objc_msgSend(*); -> https://developer.apple.com/documentation/objectivec/1456712-objc_msgsend?language=objc
    NativeLong objc_msgSend(NativeLong receiver, Pointer selector, int arg1);

    // void objc_msgSend(*); -> https://developer.apple.com/documentation/objectivec/1456712-objc_msgsend?language=objc
    NativeLong objc_msgSend(Pointer receiver, Pointer selector, NativeLong arg1);

    // id objc_getClass(const char* name); -> https://developer.apple.com/documentation/objectivec/1418952-objc_getclass?language=objc
    Pointer objc_getClass(String className);

    // SEL sel_registerName(const char* name); -> https://developer.apple.com/documentation/objectivec/1418557-sel_registername?language=objc
    Pointer sel_registerName(String selectorName);

    int CFStringGetLength(NativeLong theString);

    byte CFStringGetCString(NativeLong theString, byte[] buffer, int bufferSize, int encoding);

}