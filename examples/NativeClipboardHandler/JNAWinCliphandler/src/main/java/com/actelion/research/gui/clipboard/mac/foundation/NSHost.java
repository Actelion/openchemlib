package com.actelion.research.gui.clipboard.mac.foundation;

import com.sun.jna.NativeLong;
import com.sun.jna.Pointer;

public class NSHost extends NSObject {
    // NSHost -> https://developer.apple.com/documentation/foundation/nshost?language=objc
    private static final Pointer nativeClass = Foundation.INSTANCE.objc_getClass("NSHost");

    // [NSHost currentHost]; -> https://developer.apple.com/documentation/foundation/nshost/1408946-currenthost
    private static final Pointer currentHostSelector = Foundation.INSTANCE.sel_registerName("currentHost");

    // [NSHost name]; -> https://developer.apple.com/documentation/foundation/nshost/1416949-name
    private static final Pointer nameSelector = Foundation.INSTANCE.sel_registerName("name");

    public NSHost(NativeLong id) {
        super(id);
    }

    /**
     * Returns an NSHost instance for this machine
     * <p>
     * [NSHost currentHost];
     */
    public static NSHost currentHost() {
        return new NSHost(Foundation.INSTANCE.objc_msgSend(nativeClass, currentHostSelector));
    }

    /**
     * Returns the computer's hostname
     *
     * @return null if objc_msgSend returns 0, otherwise an NSString instance
     */
    public NSString getName() {
        NativeLong nativeResult = Foundation.INSTANCE.objc_msgSend(getId(), nameSelector);

        long result = nativeResult.longValue();
        return result != 0 ? new NSString(nativeResult) : null;
    }
}