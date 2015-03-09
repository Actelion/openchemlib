/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Christian Rufener
 */

package com.actelion.research.util;


import java.awt.*;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

/**
 * Requires JRE 1.6
 */
// Think about renaming this class
public class Platform
{
    private static boolean isWin = (System.getProperty("os.name").toLowerCase().indexOf("win") >= 0);
    private static boolean isMac = (System.getProperty("os.name").toLowerCase().indexOf("mac") >= 0);
    private static boolean isLinux = (System.getProperty("os.name").toLowerCase().indexOf("nux") >= 0);
    private static boolean isUnix = (System.getProperty("os.name").toLowerCase().indexOf("nix") >= 0);

    private static Properties PATH_PROPERTIES = null;

    private static String WINPATH = "\\\\actelch02\\pgm\\ActelionResearch";
    private static String WINPROPERTIES = WINPATH + "\\ApplicationDatabase.db";

    private static final String[][] MACINTOSH_APPLICATION_NAME = { { "orbit" , "Orbit Image Analysis" },
    															   { "datawarrior", "DataWarrior" },
    															   { "actelion3d", "Actelion3D" } };

    public static boolean isWindows()
    {
        return isWin;
    }

    public static boolean isMacintosh()
    {
        return isMac;
    }

    public static boolean isLinux()
    {
        return isLinux;
    }

    public static boolean isUnix()
    {
        return isUnix;
    }


    /**
     * Start an executable with parameters. In non-windows environment this is searching the PATH
     * In the windows environment @Actelion this first checks the standard locations where applications
     * are installed usually \\actelch02\pgm. This is achieved by consulting the ApplicationDatabase.db file
     * located under \\actelch02\pgm\ActelionResearch. The file format is a standard Java properties file
     * where the key contains the name of the application and the value the respective absolute path of the executable
     * e.g.
     * datawarrior=//actelch02/pgm/Datawarrior/DataWarrior.exe
     * @param program
     * @param args
     * @throws IOException
     */
    public static void execute(String program, String... args) throws IOException
    {
        String executable = findExecutable(program);
        List<String> arguments = new ArrayList<String>();
        arguments.add(executable);
        if (args != null && args.length > 0) {
            for (String a : args)
                arguments.add(a);
        }
        Runtime.getRuntime().exec(arguments.toArray(new String[0]));
    }

    /**
     * Given a filename, open this file with the default application
     * @param doc File name to open
     * @throws IOException if the file is unavailable
     * @throws UnsupportedOperationException, if Desktop integration of the platform is unavailable
     */
    public static void openDocument(String doc) throws IOException, UnsupportedOperationException
    {
        if (!Desktop.isDesktopSupported()) {
            if (isLinux())
                throw new UnsupportedOperationException("Please check your OS installation: Is the libgnome2 library missing?");
            else
                throw new UnsupportedOperationException("Desktop integration on this OS/Machine not supported");
        } else {
            Desktop.getDesktop().open(new File(doc));
        }
    }

    /**
     * In a platform specific way tries to find the path to the application.
     * @param name
     * @return valid full path to application or null
     */
    private static String findExecutable(String name)
    {
        String res = name;
        if (isWindows()) {
            try {
                //load a properties file
                if (PATH_PROPERTIES == null) {
                    PATH_PROPERTIES = new Properties();
                    PATH_PROPERTIES.load(new FileInputStream(WINPROPERTIES));
                }
                String value = PATH_PROPERTIES.getProperty(name);
                if (value != null)
                    res = value;
            } catch (Throwable e) {
                System.err.println("Error reading Application Database file: " + e);
            }
        }
        if (isMacintosh()) {
        	for (String[] appKeyAndName:MACINTOSH_APPLICATION_NAME) {
        		if (appKeyAndName[0].equals(name)) {
        			String path = "/Applications/"+appKeyAndName[1]+".app/Contents/MacOS/JavaApplicationStub";
        			return new File(path).exists() ? path : null;
        		}
        	}
        }
        if (isLinux()) {
        	String path = "/opt/actelion/"+name+"/"+name;
        	if (new File(path).exists())
        		return path;
        	path = "/opt/"+name+"/"+name;
        	if (new File(path).exists())
        		return path;
        	return null;
        }
        return res;
    }

/*
    public static void main(String args[]) throws Exception
    {
        if (isWindows()) {
            openDocument("C:\\temp\\supplier.txt");
            openDocument("C:\\temp\\out.elb");
            openDocument("C:\\temp\\foo.sdf");
            execute("Notepad");
            execute("cmd", "/c", "start");
            execute("nemo");
            execute("datawarrior","c:\\temp\\foo.sdf");
        }
        if (isWindows()) {
        	
        }
    }
*/
}
