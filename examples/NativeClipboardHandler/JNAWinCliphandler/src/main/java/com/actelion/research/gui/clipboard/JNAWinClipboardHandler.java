package com.actelion.research.gui.clipboard;

import com.sun.jna.Native;
import com.sun.jna.Pointer;
import com.sun.jna.Structure;
import com.sun.jna.platform.win32.*;
import com.sun.jna.win32.W32APIOptions;

import java.io.*;
import java.util.HashMap;
import java.util.Map;


public class JNAWinClipboardHandler {

    private static boolean isInitOK = true; // dummy to match NativeClipboardAccessor
    // Numeric windows default formats are defined in jna. E.g. User32.CF_BITMAP or User32.CF_DIB
    private static Map<String, Integer> stdWinCFValues = new HashMap() {
        {
            put(ClipboardHandler.NC_METAFILE, User32.CF_ENHMETAFILE);
            put(ClipboardHandler.NC_DIB, User32.CF_DIB);
        }
    };

    private static int getFormatNo(String format) {
        return stdWinCFValues.get(format) == null ? USER32.RegisterClipboardFormat(format) : stdWinCFValues.get(format);
    }

    public static byte[] getClipboardData(String format) {
        byte[] data = null;
        int x = getFormatNo(format);
        data = getRawData(x);
        return data;
    }

    private static byte[] getRawData(int clipFormatNo) {
        byte[] buffer = null;
        if (USER32.IsClipboardFormatAvailable(clipFormatNo)) {
            if (USER32.OpenClipboard((WinDef.HWND)null)) {
                try {
                    WinNT.HANDLE handle = USER32.GetClipboardData(clipFormatNo);
                    if (handle != null) {
                        int size = KERNEL32.GlobalSize(handle);
                        Pointer pt = KERNEL32.GlobalLock(handle);
                        buffer = pt.getByteArray(0, size);
                        KERNEL32.GlobalUnlock(handle);
                    }
                } finally {
                    USER32.CloseClipboard();
                }
            }
        }
        return buffer;
    }

    private static WinNT.HANDLE createEnhMetaFile(WinNT.HANDLE hmeta/*, int x, int y*/) {
        WinNT.HANDLE henhmeta = null;
        /*if (x > 0 && y > 0) { // code block in dll, needed to copy sized molecule, but could not find a call, so maybe not required anymore?
            struct POINT would be required.
        } */
        int size = GDI32.GetMetaFileBitsEx(hmeta, 0, null);
        if (size > 0) {
            WinNT.HANDLE handle = KERNEL32.GlobalAlloc(0x0042, size);
            Pointer ptr = KERNEL32.GlobalLock(handle);
            GDI32.GetMetaFileBitsEx(hmeta, size, ptr);
            henhmeta = GDI32.SetWinMetaFileBits(size, ptr, null, null);
            KERNEL32.GlobalUnlock(handle);
            Kernel32 krn = Kernel32.INSTANCE;
            krn.GlobalFree(ptr);
        }
        return henhmeta;
    }

    public static boolean setClipBoardData(String format, byte[] buffer) {
        int fno = getFormatNo(format);
        return setClipBoardData(fno, buffer, true);
    }
    public static boolean setClipBoardData(String format, byte[] buffer, boolean emptyClipboard) {
        int fno = getFormatNo(format);
        return setClipBoardData(fno, buffer, emptyClipboard);
    }

    public static boolean setClipBoardData(int format, byte[] buffer, boolean emptyClipboard) {
        boolean ok = false;
        if (USER32.OpenClipboard(null)) {
            try {
                if (emptyClipboard) {
                    USER32.EmptyClipboard();
                }
                WinNT.HANDLE handleData = formatHandling(format, buffer);
                if (handleData != null) {
                    USER32.SetClipboardData(format, handleData);
                    ok = true;
                }
            } finally {
                USER32.CloseClipboard();
            }
        }
        return ok;
    }

    public static boolean setClipBoardData(int format, byte[] buffer) {
        return setClipBoardData(format, buffer, true);
    }

    private static WinNT.HANDLE formatHandling(int format, byte[] buffer) {
        WinNT.HANDLE handleData = null;
        // TODO add other handling from dll if required, like DIB
        if (format == USER32.CF_ENHMETAFILE) {
            try {
                handleData = wmfToEmf(buffer);
            } catch (IOException e) {
                e.printStackTrace();
            }
        } else {
            handleData = copyToGlobalBuffer(buffer);
        }
        return handleData;
    }

    private static WinNT.HANDLE copyToGlobalBuffer(byte[] b) {
        WinNT.HANDLE handle = KERNEL32.GlobalAlloc(0x0042, b.length);
        Pointer pt = KERNEL32.GlobalLock(handle);
        pt.write(0, b, 0, b.length);
        KERNEL32.GlobalUnlock(handle);
        return handle;
    }

    public static boolean copyMoleculeToClipboard(String filename, byte[] sketch, byte[] serializedObject) {
        boolean ok = false;
        if (USER32.OpenClipboard(null)) {
            try {
                if (USER32.EmptyClipboard()) {
                    WinNT.HANDLE hmeta = wmfToEmf(filename);
                    int fmol = USER32.RegisterClipboardFormat(ClipboardHandler.NC_SERIALIZEMOLECULE);
                    int fcd = USER32.RegisterClipboardFormat(ClipboardHandler.NC_CHEMDRAWINTERCHANGE);

                    WinNT.HANDLE handleMol = copyToGlobalBuffer(serializedObject);
                    WinNT.HANDLE handleCD = copyToGlobalBuffer(sketch);


                    if (hmeta != null) {
                        USER32.SetClipboardData(USER32.CF_ENHMETAFILE, hmeta);
                    }
                    if (handleMol != null) {
                        USER32.SetClipboardData(fmol, handleMol);
                    }
                    if (handleCD != null) {
                        USER32.SetClipboardData(fcd, handleCD);
                    }

                    ok = true;
                }
            } finally {
                USER32.CloseClipboard();
            }
        }
        return ok;
    }

    public static boolean copyReactionToClipboard(byte[] ctab, byte[] sketch, byte[] serializedObject) {
        boolean ok = false;
        if (USER32.OpenClipboard(null)) {
            try {
                if (USER32.EmptyClipboard()) {

                    int fctab = USER32.RegisterClipboardFormat(ClipboardHandler.NC_CTAB);
                    int frxn = USER32.RegisterClipboardFormat(ClipboardHandler.NC_SERIALIZEREACTION);
                    int fcd = USER32.RegisterClipboardFormat(ClipboardHandler.NC_CHEMDRAWINTERCHANGE);

                    WinNT.HANDLE handlectab = copyToGlobalBuffer(ctab);
                    WinNT.HANDLE handleRxn = copyToGlobalBuffer(serializedObject);
                    WinNT.HANDLE handleCD = copyToGlobalBuffer(sketch);

                    USER32.SetClipboardData(fctab, handlectab );
                    USER32.SetClipboardData(fcd, handleCD);
                    USER32.SetClipboardData(frxn, handleRxn);

                    ok = true;

                }
            } finally {
                USER32.CloseClipboard();
            }
        }
        return ok;
    }

    public static WinNT.HANDLE wmfToEmf(String path) {
        WinNT.HANDLE hmeta = GDI32.GetEnhMetaFile(path);
        if (hmeta == null) {
            hmeta = GDI32.GetMetaFile(path);
            hmeta = createEnhMetaFile(hmeta);
        }
        return hmeta;
    }

    public static WinNT.HANDLE wmfToEmf(byte[] buffer) throws IOException {
        File f = File.createTempFile("wmf2emf", "tmp");
        FileOutputStream out = new FileOutputStream(f);
        out.write(buffer);
        out.close();
        WinNT.HANDLE handle = wmfToEmf(f.getPath());
        f.delete();
        return handle;
    }

    public static boolean isInitOK() {return isInitOK;}
    public interface MyUser32 extends User32 {
        MyUser32 INSTANCE = (MyUser32) Native.load("user32", MyUser32.class, W32APIOptions.UNICODE_OPTIONS);
        boolean EmptyClipboard();

        int RegisterClipboardFormat(String var1);

        boolean OpenClipboard(WinDef.HWND var1);

        boolean CloseClipboard();

        WinNT.HANDLE SetClipboardData(int var1, WinNT.HANDLE var2);

        WinNT.HANDLE GetClipboardData(int var1);

        boolean IsClipboardFormatAvailable(int var1);
    }

    public interface MyKernel32 extends Kernel32 {

        MyKernel32 INSTANCE = (MyKernel32) Native.load("kernel32", MyKernel32.class, W32APIOptions.UNICODE_OPTIONS);
        WinNT.HANDLE GlobalAlloc(int var1, int var2);

        Pointer GlobalLock(WinNT.HANDLE var1);

        void GlobalUnlock(WinNT.HANDLE var1);

        int GlobalSize(WinNT.HANDLE var1);
    }

    public interface MyGDI32 extends GDI32 {
        MyGDI32 INSTANCE = (MyGDI32) Native.load("gdi32", MyGDI32.class, W32APIOptions.UNICODE_OPTIONS);
        WinNT.HANDLE GetEnhMetaFile(String path);
        WinNT.HANDLE GetMetaFile(String path);
        int GetMetaFileBitsEx(WinNT.HANDLE hmf, int buffsize, Pointer buff);
        WinNT.HANDLE SetWinMetaFileBits(int size, Pointer bytearrayptr, WinDef.HDC hdc, Structure struct);
    }

    private static MyUser32 USER32 = MyUser32.INSTANCE;
    private static MyKernel32 KERNEL32 = MyKernel32.INSTANCE;
    private static MyGDI32 GDI32 = MyGDI32.INSTANCE;

}
