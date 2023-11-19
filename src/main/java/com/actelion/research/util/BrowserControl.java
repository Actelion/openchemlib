package com.actelion.research.util;

import java.awt.*;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

public class BrowserControl {
	/**
	 * Display a file in the system browser.  If you want to display a
	 * file, you must include the absolute path name.
	 *
	 * @param url the file's url (the url must start with either "http://"
	 * or "file://").
	 */
	public static void displayURL(String url) {
		new Thread(() -> {
			try {
				if (Desktop.isDesktopSupported() && Desktop.getDesktop().isSupported(Desktop.Action.BROWSE)) {
					try {
						Desktop.getDesktop().browse(new URI(url));
					} catch(URISyntaxException use) {
						use.printStackTrace();
						}
					}
				else if (Platform.isWindows()) {
					Runtime.getRuntime().exec("rundll32 url.dll,FileProtocolHandler " + url);
				} else if (Platform.isMacintosh()) {
					Runtime.getRuntime().exec("open " + url);
				}
				else {
					Runtime rt = Runtime.getRuntime();
					String[] browsers = { "google-chrome", "firefox", "mozilla", "epiphany", "konqueror",
							"netscape", "opera", "links", "lynx" };

					StringBuffer cmd = new StringBuffer();
					for (int i = 0; i < browsers.length; i++)
						if(i == 0)
							cmd.append(String.format(    "%s \"%s\"", browsers[i], url));
						else
							cmd.append(String.format(" || %s \"%s\"", browsers[i], url));
					// If the first didn't work, try the next browser and so on

					try {
						rt.exec(new String[]{"sh", "-c", cmd.toString()});
					} catch (IOException ioe) {
						ioe.printStackTrace();
					}
				}
			} catch (IOException x) {
				System.err.println("Could not invoke browser: "+x);
			}
		} ).start();
	}
}