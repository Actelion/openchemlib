package com.actelion.research.gui.clipboard.mac.foundation;

import com.sun.jna.NativeLong;
import com.sun.jna.Pointer;

/**
 * A reference to an NSObject, which is holds an 'identifier' (a reference) to a specific object
 */
public class NSObject {
    private static final NativeLong nullPointer = new NativeLong(0L);

    // [NSObject description]; -> https://developer.apple.com/documentation/objectivec/1418956-nsobject/1418746-description
    private static final Pointer descriptionSelector = Foundation.INSTANCE.sel_registerName("description");

    private final NativeLong id;

    public NSObject(NativeLong id) {
        this.id = id;
    }

    /**
     * Returns the identifier for this NSObject
     *
     * @return null pointer if the object is null, otherwise the identifier
     */
    public NativeLong getId() {
        return id;
    }

    /**
     * A textual representation of the receiver
     *
     * @return A string that describes the object
     */
    public NSString getDescription() {
        return new NSString(Foundation.INSTANCE.objc_msgSend(id, descriptionSelector));
    }

    public boolean isNull() {
        return id.equals(nullPointer);
    }

    @Override
    public String toString() {
        return getDescription().getJvmString();
    }

    /**
     * A reference to an NSObject, which is holds an 'identifier' (a reference) to a specific object which is
     * automatically released when the JVM invokes a garbage collection.
     * <p>
     * This *only* be used when you allocate the NSObject, for example: NSData or NSString.
     */
    public static class Releasable extends NSObject implements AutoCloseable {
        // [NSObject release]; -> https://developer.apple.com/documentation/objectivec/1418956-nsobject/1571957-release
        private static final Pointer releaseSelector = Foundation.INSTANCE.sel_registerName("release");

        public Releasable(NativeLong id) {
            super(id);
        }

        @Override
        public void close() {
            if (System.getProperty("jnapple.debug", "false").equals("true"))
                System.out.println("[d] Cleaning up NSObject: " + getId());

            Foundation.INSTANCE.objc_msgSend(getId(), releaseSelector);
        }
    }
}