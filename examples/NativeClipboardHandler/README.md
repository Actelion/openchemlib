## Native Clipboard Handler for OpenChemLib
*OpenChemLib* provides Clipboard functionality like Copy & Paste within its Editors by default. However, to access Data from third parties, like ChemDraw, a native Clipboard Handler is required.
Its best to look at com.actelion.research.gui.clipboard.ClipboardHandler in *OpenChemLib* and these examples to see how things work.
If you build a one of the cliphandler projects and put the resulting jar and JNA dependencies on the classpath, *OpenChemLib* will try to use the native implementation, if executed on the Operating System for which it is implemented.

*OpenChemLib* with a version later then *CHANGEME* will load the native Clipboardhandler by the fully qualified class name defined in the Clipboardhandler class (nativeCliphandler list).
If more than one are defined it will use the first one available. If none is available, *OpenChemLib* will fall back to its java implementation.
A native Clipboardhandler implementation should implement the getClipboardData(String format), setClipBoardData(String format, byte[] buffer, boolean emptyClipboard) and isInitOK() methods as shown in the examples.  
On windows for compatibility with a clipboardhandler dll the methods copyReactionToClipboard and copyMoleculeToClipboard can be implemented, they will not offer more functionality

### Dependencies
These examples use JNA to provide native Clipboard Access. It should work from JNA Version 4.0.0, however a more recent version, preferably the latest release (5.18.0 at time of writing), is recommended.

