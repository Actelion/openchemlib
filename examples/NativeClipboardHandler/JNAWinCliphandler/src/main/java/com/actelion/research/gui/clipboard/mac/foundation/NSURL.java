package com.actelion.research.gui.clipboard.mac.foundation;

import com.sun.jna.NativeLong;
import com.sun.jna.Pointer;

import java.net.MalformedURLException;
import java.net.URL;

/**
 * A reference to an NSURL, which is an object that represents the location of a resource, such as an item on a remote
 * server or the path to a local file.
 */
public class NSURL extends NSObject.Releasable {
    // NSURL -> https://developer.apple.com/documentation/foundation/nsurl?language=objc
    private static final Pointer nativeClass = Foundation.INSTANCE.objc_getClass("NSURL");

    // [NSURL URLWithString:]; -> https://developer.apple.com/documentation/foundation/nsurl/1572047-urlwithstring?language=objc
    private static final Pointer urlWithStringSelector = Foundation.INSTANCE.sel_registerName("URLWithString:");

    // [NSURL absoluteString]; -> https://developer.apple.com/documentation/foundation/nsurl/1409868-absolutestring?language=objc
    private static final Pointer absoluteStringSelector = Foundation.INSTANCE.sel_registerName("absoluteString");

    public NSURL(NativeLong id) {
        super(id);
    }

    /**
     * Creates and returns an NSURL object initialized with a provided URL string.
     *
     * @param url The URL string with which to initialize the NSURL object. Must be a URL that conforms to RFC 2396.
     *            This method parses URLString according to RFCs 1738 and 1808.
     */
    public NSURL(String url) {
        super(Foundation.INSTANCE.objc_msgSend(nativeClass, urlWithStringSelector, new NSString(url).getId()));
    }

    /**
     * Returns a Java {@link URL} from an NSURL instance
     *
     * @return the NSURL converted to URL
     */
    public URL getJvmURL() throws MalformedURLException {
        return new URL(getAbsoluteString().getJvmString());
    }

    public NSString getAbsoluteString() {
        return new NSString(Foundation.INSTANCE.objc_msgSend(getId(), absoluteStringSelector));
    }
}