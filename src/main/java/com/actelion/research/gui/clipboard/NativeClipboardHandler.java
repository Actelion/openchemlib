/*
* Copyright (c) 1997 - 2016
* Actelion Pharmaceuticals Ltd.
* Gewerbestrasse 16
* CH-4123 Allschwil, Switzerland
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 3. Neither the name of the the copyright holder nor the
*    names of its contributors may be used to endorse or promote products
*    derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/

package com.actelion.research.gui.clipboard;


public class NativeClipboardHandler
{

    private static boolean isWin = (System.getProperty("os.name").toLowerCase().indexOf("win") >= 0);
    private static boolean isMac = (System.getProperty("os.name").toLowerCase().indexOf("mac") >= 0);
    private static boolean isLinux = (System.getProperty("os.name").toLowerCase().indexOf("nux") >= 0);
    private static boolean isUnix = (System.getProperty("os.name").toLowerCase().indexOf("nix") >= 0);



    public static final String NC_SKETCH	= "MDLSK";
    public static final String NC_CTAB		= "MDLCT";
    public static final String NC_MOLFILE	= "MDL_MOL";
    public static final String NC_METAFILE	= "CF_METAFILEPICT";
    public static final String NC_DIB		= "CF_DIB";
    public static final String NC_BITMAP	= "CF_BITMAP";
    public static final String NC_SERIALIZEMOLECULE = "ACT_MOLECULE";
    public static final String NC_SERIALIZEREACTION = "ACT_REACTION";
    public static final String NC_ALDUSMETAFILE	= "ALDUS_METAFILE";;
    public static final String NC_EMBEDDEDSKETCH = "MDLSK_EMBEDDED";
    public static final String NC_CHEMDRAWINTERCHANGE = "ChemDraw Interchange Format";
    public static final String NC_IDCODE		= "IDCODE";



    static boolean isWindows()
    {
        return isWin;
    }

    static boolean isMacintosh()
    {
        return isMac;
    }

    static boolean isLinux()
    {
        return isLinux;
    }

    static boolean isUnix()
    {
        return isUnix;
    }



    public static boolean copyMoleculeToClipboard(String filename,byte[] sketch, byte[] serializedObject)
    {
        if (isWindows()) {
            return NativeClipboardAccessor.copyMoleculeToClipboard(filename,sketch,serializedObject);
        } else if (isLinux() || isMacintosh()) {
            return LinuxNativeClipboardAccessor.copyMoleculeToClipboard(filename,sketch,serializedObject);
        } else
            return false;
    }

    public static boolean copyReactionToClipboard(String filename,byte[] sketch, byte[] serializedObject)
    {
        if (isWindows()) {
            return NativeClipboardAccessor.copyReactionToClipboard(filename,sketch,serializedObject);
        } else if (isLinux() || isMacintosh()) {
            return LinuxNativeClipboardAccessor.copyReactionToClipboard(filename,sketch,serializedObject);
        } else
            return false;

    }

    // public static native boolean copyMoleculeToClipboard(String filname,byte[] sketch, byte[] serializedObject);
        /* Formats are "MDLSK","MDLCT","MDL_MOL","CF_METAFILEPICT","CF_DIB" "ACT_MOLECULE" */
    public static byte[] getClipboardData(String format)
    {
        if (isWindows()) {
            return NativeClipboardAccessor.getClipboardData(format);
        } else if (isLinux() || isMacintosh()) {
            return LinuxNativeClipboardAccessor.getClipboardData(format);
        } else
            return null;
    }

    public static boolean setClipBoardData(String format, byte[] buffer)
    {
        if (isWindows()) {
            return NativeClipboardAccessor.setClipBoardData(format,buffer);
        } else if (isLinux() || isMacintosh()) {
            return LinuxNativeClipboardAccessor.setClipBoardData(format,buffer);
        } else
            return false;


    }

    /**
     * Copy a windows enhance metafile to the clipboard
     * @param data byte[]
     * @return boolean
     */
    public static boolean copyMetaFile(byte []data)
    {
        java.awt.datatransfer.Clipboard clip = java.awt.Toolkit.getDefaultToolkit().getSystemClipboard();
        return setClipBoardData(NC_METAFILE,data);
    }

    /**
     * Copies an Image as CF_DIB onto the clipboard. Currently it is using the JIMI framework to perform the DIB conversion
     * @param img   Image to be copies
     * @param width Image Width
     * @param height Image Height
     * @return
     */
/*
    public static boolean copyImage(Image img)
    {
        boolean ok = false;
        int FILE_HEADER_SIZE = 14;
        String clsname = "com.sun.jimi.core.Jimi";
        String methodname = "putImage";
        String id = "image/bmp";
        try {
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            Class cls = Class.forName(clsname);
            java.lang.reflect.Method m = cls.getMethod(methodname,new Class[] { String.class,Image.class,OutputStream.class});
            m.invoke(null,new Object[]{id,img,baos});
            ByteArrayOutputStream b = new ByteArrayOutputStream();
            int size = baos.size();
            if (size > FILE_HEADER_SIZE) {
                b.write(baos.toByteArray(),FILE_HEADER_SIZE,size-FILE_HEADER_SIZE);
                baos = null;
                ok = setClipBoardData(NC_DIB,b.toByteArray());
                b = null;
            }
        } catch (java.lang.LinkageError le) {
            System.err.println("NativeClipboard.copyImage(): Unable to load class: " + clsname);
        } catch (ClassNotFoundException ce) {
            System.err.println("NativeClipboard.copyImage(): Unable to load class: " + clsname);
        } catch (NoSuchMethodException me) {
            System.err.println("NativeClipboard.copyImage(): Method Not Found : " + clsname);
        } catch (java.lang.reflect.InvocationTargetException ie) {
            System.err.println("NativeClipboard.copyImage(): " + ie);
        } catch (IllegalAccessException ae) {
            System.err.println("NativeClipboard.copyImage(): " + ae);
        } catch (Exception te) {
            System.err.println("NativeClipboard.copyImage(): " + te);
        }
        return ok;
    }
*/
    /**
     * Paste a DIB from the clipboard. Currently it is using the JIMI framework to perform the DIB conversion
     * @return Image found on the clipboard or null if there was not DIB or the conversion failed.
     */
/*
    public static java.awt.Image pasteImage()
    {
        byte[] buffer = com.actelion.research.gui.NativeClipboardHandler.pasteDIB();
        if (buffer != null) {
            return pasteDIB(buffer);
        } else {
            buffer = getClipboardData(NC_METAFILE);
            com.actelion.research.util.WMFDecoder dec = new com.actelion.research.util.WMFDecoder(new ByteArrayInputStream(buffer));
            return Toolkit.getDefaultToolkit().createImage(dec);
        }
//        return null;
    }

    private static java.awt.Image pasteDIB(byte[] buffer)
    {
        java.awt.Image img = null;
        if (buffer != null) {
            String clsname = "com.sun.jimi.core.Jimi";
            String methodname = "getImage";
            String id = "image/bmp";
            try{
                Class cls = Class.forName(clsname);
                java.lang.reflect.Method m = cls.getMethod(methodname,
                        new Class[]{InputStream.class,String.class});
                if(buffer != null){
                    Object o = m.invoke(null,
                                        new Object[]{new
                                        ByteArrayInputStream(buffer),id});
                    if(o instanceof java.awt.Image){
                        img = (java.awt.Image)o;
                    }
                }
            } catch(java.lang.LinkageError le){
                System.err.println(
                        "NativeClipboard.pasteImage(): Unable to load class: " +
                        clsname);
            } catch(ClassNotFoundException ce){
                System.err.println(
                        "NativeClipboard.pasteImage(): Unable to load class: " +
                        clsname);
            } catch(NoSuchMethodException me){
                System.err.println(
                        "NativeClipboard.pasteImage(): Method Not Found : " +
                        clsname);
            } catch(java.lang.reflect.InvocationTargetException ie){
                System.err.println("NativeClipboard.pasteImage(): " + ie);
            } catch(IllegalAccessException ae){
                System.err.println("NativeClipboard.pasteImage(): " + ae);
            } catch(Exception te){
                System.err.println("NativeClipboard.pasteImage(): " + te);
            }
        }
        return img;
    }

    private static byte[] pasteDIB()
    {
        return getClipboardData(NC_DIB);
    }

*/
}
