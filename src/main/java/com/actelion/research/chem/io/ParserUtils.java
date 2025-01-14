package com.actelion.research.chem.io;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URL;
import java.util.Arrays;

/**
 * 
 * @author hanson@stolaf.edu
 *
 */
public class ParserUtils {

	private ParserUtils() {
	}

	public static byte[] getResourceBytes(Class<?> c, String fileName) throws FileNotFoundException, IOException {
		return getLimitedStreamBytes(c.getResourceAsStream(fileName), -1, null, true, true);
	}

	public static byte[] getLimitedStreamBytes(InputStream is, long n, OutputStream out, boolean andCloseInput,
			boolean andCloseOutput) throws IOException {

		// Note: You cannot use InputStream.available() to reliably read
		// zip data from the web.

		boolean toOut = (out != null);
		int buflen = (n > 0 && n < 1024 ? (int) n : 1024);
		byte[] buf = new byte[buflen];
		byte[] bytes = (out == null ? new byte[n < 0 ? 4096 : (int) n] : null);
		int len = 0;
		int totalLen = 0;
		if (n < 0)
			n = Integer.MAX_VALUE;
		while (totalLen < n && (len = is.read(buf, 0, buflen)) > 0) {
			totalLen += len;
			if (toOut) {
				out.write(buf, 0, len);
			} else {
				if (totalLen > bytes.length)
					bytes = Arrays.copyOf(bytes, totalLen * 2);
				System.arraycopy(buf, 0, bytes, totalLen - len, len);
				if (n != Integer.MAX_VALUE && totalLen + buflen > bytes.length)
					buflen = bytes.length - totalLen;
			}
		}
		if (andCloseInput) {
			try {
				is.close();
			} catch (IOException e) {
				// ignore
			}
		}
		if (toOut) {
			if (andCloseOutput)
				try {
					out.close();
				} catch (IOException e) {
					// ignore
				}
			return null;
		}
		if (totalLen == bytes.length)
			return bytes;
		buf = new byte[totalLen];
		System.arraycopy(bytes, 0, buf, 0, totalLen);
		return buf;
	}

	public static String getURLContentsAsString(String url) {
		byte[] b = getURLContentsAsBytes(url);
		return (b == null ? null : new String(b));
	}

	public static byte[] getURLContentsAsBytes(String url) {
		try {
			url = ensureURLPath(url);
			return getLimitedStreamBytes(new URL(url).openStream(), -1, null, true, true);
		} catch (Exception e) {
			return null;
		}
	}

	private static String ensureURLPath(String url) {
		return (url.indexOf("://") < 0 ? "file://" + (url.startsWith("/") ? "" : "/") + url : url);
	}

}
