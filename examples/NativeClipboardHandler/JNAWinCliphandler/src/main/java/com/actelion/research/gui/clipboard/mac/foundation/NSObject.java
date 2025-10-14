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

/**
 * A reference to an NSObject, which is holds an 'identifier' (a reference) to a specific object
 */
public class NSObject
{
    private static final NativeLong nullPointer = new NativeLong(0L);

    // [NSObject description]; -> https://developer.apple.com/documentation/objectivec/nsobject/description
    private static final Pointer descriptionSelector = Foundation.INSTANCE.sel_registerName("description");

    private final NativeLong id;

    public NSObject(NativeLong id)
    {
        this.id = id;
    }

    /**
     * Returns the identifier for this NSObject
     *
     * @return null pointer if the object is null, otherwise the identifier
     */
    public NativeLong getId()
    {
        return id;
    }

    /**
     * A textual representation of the receiver
     *
     * @return A string that describes the object
     */
    public NSString getDescription()
    {
        return new NSString(Foundation.INSTANCE.objc_msgSend(id, descriptionSelector));
    }

    public boolean isNull()
    {
        return id.equals(nullPointer);
    }

    @Override
    public String toString()
    {
        return getDescription().getJvmString();
    }

    /**
     * A reference to an NSObject, which is holds an 'identifier' (a reference) to a specific object which is
     * automatically released when the JVM invokes a garbage collection.
     * <p>
     * This *only* be used when you allocate the NSObject, for example: NSData or NSString.
     */
    public static class Releasable extends NSObject implements AutoCloseable
    {
        // [NSObject release]; -> https://developer.apple.com/documentation/objectivec/nsobject/release
        private static final Pointer releaseSelector = Foundation.INSTANCE.sel_registerName("release");

        public Releasable(NativeLong id)
        {
            super(id);
        }

        @Override
        public void close()
        {
            if (System.getProperty("jnapple.debug", "false").equals("true"))
                System.out.println("[d] Cleaning up NSObject: " + getId());

            Foundation.INSTANCE.objc_msgSend(getId(), releaseSelector);
        }
    }
}