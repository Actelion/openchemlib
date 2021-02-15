package com.actelion.research.gui.dock;

import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;

public class TransferableDockable implements Transferable {
	public static final DataFlavor DF_DOCKABLE_DEF = new DataFlavor("application/x-openmolecules-dockable;class=java.lang.String", "Dockable view setting");
	public static final DataFlavor[] DOCKABLE_FLAVORS = { DF_DOCKABLE_DEF };

	private Dockable mDockable;

	public TransferableDockable(Dockable dockable) {
		mDockable = dockable;
		}

	@Override
	public DataFlavor[] getTransferDataFlavors() {
		return DOCKABLE_FLAVORS;
		}

	@Override
	public boolean isDataFlavorSupported(DataFlavor flavor) {
		for (DataFlavor df:DOCKABLE_FLAVORS)
			if (df.equals(flavor))
				return true;
		return false;
		}

	@Override
	public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException {
		if (flavor.equals(DF_DOCKABLE_DEF))
			return mDockable.getTitle();    // for inter app drag&drop we will need the full configuration

		throw new UnsupportedFlavorException(null);
		}
	}
