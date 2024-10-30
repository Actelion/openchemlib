package com.actelion.research.gui.table;

import com.actelion.research.gui.LookAndFeelHelper;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.actelion.research.util.ColorHelper;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.TableColumnModel;
import javax.swing.table.TableModel;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.Vector;

public class JTableWithRowNumbers extends JTable implements TableModelListener {
	protected static Cursor	sResizeCursor,sDefaultCursor;

	private static final long serialVersionUID = 0x20060906;

	private RowNumberPanel	mRowNumberPanel;
	private JScrollPane		mScrollPane = null;
	private int				mRowHeaderMinWidth;

	public JTableWithRowNumbers() {
		super();
		initialize();
		}

    public JTableWithRowNumbers(int numRows, int numColumns) {
		super(numRows, numColumns);
		initialize();
		}

    public JTableWithRowNumbers(final Object[][] rowData, final Object[] columnNames) {
		super(rowData, columnNames);
		initialize();
		}

    public JTableWithRowNumbers(TableModel dm) {
		super(dm);
		initialize();
		}

    public JTableWithRowNumbers(TableModel dm, TableColumnModel cm) {
		super(dm, cm);
		initialize();
		}

    public JTableWithRowNumbers(TableModel dm, TableColumnModel cm, ListSelectionModel sm) {
		super(dm, cm, sm);
		initialize();
		}

	public JTableWithRowNumbers(final Vector rowData, final Vector columnNames) {
		super(rowData, columnNames);
		initialize();
		}

	private void initialize() {
		mRowNumberPanel = new RowNumberPanel(this);

		setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		sDefaultCursor = Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR);
		sResizeCursor = Cursor.getPredefinedCursor(Cursor.S_RESIZE_CURSOR);
		}

	public void addHighlightListener(HighlightListener l) {
		mRowNumberPanel.addHighlightListener(l);
	}

	public void removeHighlightListener(HighlightListener l) {
		mRowNumberPanel.removeHighlightListener(l);
	}

	public int getHighlightedRow() {
		return mRowNumberPanel.getHighlightedRow();
		}

	public void setHighlightedRow(int row, boolean informListeners) {
		mRowNumberPanel.setHighlightedRow(row, informListeners);
		}

	public JScrollPane getScrollPane() {
		return mScrollPane;
		}

	/**
	 * @return font size of parent table, which includes UI scaling factor
	 */
	public int getFontSize() {
		return getFont().getSize();
		}

	/**
	 * @param fontSize including UI scaling factor
	 */
	public void setFontSize(int fontSize) {
		Font font = new Font(Font.SANS_SERIF, Font.PLAIN, fontSize);
		super.setFont(font);
		getTableHeader().setFont(font);
		if (mRowNumberPanel != null) {
			mRowNumberPanel.setFont(font);
			validateRowNumberColumnWidth();
			}
		}

	/**
	 * @param height
	 */
	@Override
	public void setRowHeight(int height) {
		super.setRowHeight(height);
		if (mRowNumberPanel != null)
			mRowNumberPanel.repaint();
		}

	/**
	 * @param row
	 * @param height
	 */
	@Override
	public void setRowHeight(int row, int height) {
		super.setRowHeight(row, height);
		if (mRowNumberPanel != null)
			mRowNumberPanel.repaint();
		}

	@Override
	public void addNotify() {
		super.addNotify();

		if (getParent().getParent() instanceof JScrollPane) {
			mScrollPane = (JScrollPane) getParent().getParent();
			mScrollPane.setRowHeaderView(mRowNumberPanel);
			mScrollPane.getVerticalScrollBar().getModel().addChangeListener(e -> mRowNumberPanel.repaint() );

			validateRowNumberColumnWidth();

// make the background of the scrollpane match that of the table.
//			pane.getViewport().setBackground(getBackground());
//			pane.getColumnHeader().setBackground(getBackground());
//			pane.getRowHeader().setBackground(getBackground());
			}
		}

	public void setRowHeaderMinWidth(int width) {
		mRowHeaderMinWidth = width;
		validateRowNumberColumnWidth();
		}

	private void validateRowNumberColumnWidth() {
		if (mScrollPane != null) {
			JViewport viewport = mScrollPane.getRowHeader();
			Dimension size = viewport.getPreferredSize();
			int maxWidth = mRowNumberPanel.getFontMetrics(mRowNumberPanel.getFont()).stringWidth("  "+getRowCount());
			size.width = Math.max(mRowHeaderMinWidth, maxWidth);
			viewport.setPreferredSize(size);
			}
		}

	@Override
	public void tableChanged(TableModelEvent e) {
		super.tableChanged(e);

		if (mRowNumberPanel != null) {
			validateRowNumberColumnWidth();
			mRowNumberPanel.invalidate();
			mRowNumberPanel.repaint();
			}
		}

	public void rowNumberClicked(int row) {	// overwrite this if you want specific effects
		}

	public void setSelectionModel(ListSelectionModel selectionModel) {
		super.setSelectionModel(selectionModel);
		}
	}

	// Inner class used to display the row numbers on the left side of the table. Order
	// is considered important enough on this screen to re-enforce it with a visual.

class RowNumberPanel extends JPanel implements ListSelectionListener,MouseListener,MouseMotionListener {
	private static final long serialVersionUID = 0x20060906;

	private static final int cResizeTolerance = 2;

	private final JTableWithRowNumbers mTable;
	private final Vector<HighlightListener> mHighlightListeners;
	private boolean mFullRowSelection,mIsResizing;
	private int mDragStartY,mDragStartRow,mResizingStartRow,mCurrentDragRow,mDragStartRowHeight,mAchorSelectionRow;
	private int mResizingRowHeight,mHighlightedRow;

	public RowNumberPanel(JTableWithRowNumbers table) {
		mTable = table;
		mAchorSelectionRow = -1;
		mHighlightedRow = -1;
		mHighlightListeners = new Vector<>();
		addMouseListener(this);
		addMouseMotionListener(this);
		table.getSelectionModel().addListSelectionListener(this);
	}

	@Override
	public void setFont(Font font) {
		super.setFont(font.deriveFont(Font.BOLD));
	}

	private int rowFromMouseY(MouseEvent e) {
		return (e.getY() + mTable.getScrollPane().getViewport().getViewRect().y) / mTable.getRowHeight();
	}

	public void mouseClicked(MouseEvent e) {
		if (mResizingStartRow == -1)
			mTable.rowNumberClicked(rowFromMouseY(e));
	}

	public void mousePressed(MouseEvent e) {
		mIsResizing = (mResizingStartRow != -1);
		mFullRowSelection = !mIsResizing;
		mCurrentDragRow = mDragStartRow = rowFromMouseY(e);
		if (mFullRowSelection) {
			if (e.isShiftDown()) {
				if (mAchorSelectionRow != -1) {
					mTable.getSelectionModel().setSelectionInterval(mAchorSelectionRow, mCurrentDragRow);
					mAchorSelectionRow = -1;
					}
				else {
					mTable.getSelectionModel().addSelectionInterval(mCurrentDragRow, mCurrentDragRow);
					}
				}
			else if (e.isControlDown()) {
				mTable.getSelectionModel().removeSelectionInterval(mCurrentDragRow, mCurrentDragRow);
				mAchorSelectionRow = -1;
				}
			else {
				mTable.getSelectionModel().setSelectionInterval(mCurrentDragRow, mCurrentDragRow);
				mAchorSelectionRow = mCurrentDragRow;
				}
			repaint();
			}
		if (mIsResizing) {
			mDragStartY = e.getY();
			mResizingRowHeight = mDragStartRowHeight;
		}
	}

	public void mouseReleased(MouseEvent e) {
		if (mIsResizing) {
			mIsResizing = false;
			if (mDragStartRowHeight != mResizingRowHeight) {
				int scrollValue = mTable.getScrollPane().getVerticalScrollBar().getValue();
				int newScrollValue = scrollValue + mResizingStartRow * (mResizingRowHeight - mDragStartRowHeight);
				mTable.setRowHeight(mResizingRowHeight);
				mTable.getScrollPane().getVerticalScrollBar().setValue(newScrollValue);
			}
		}
	}

	public void mouseEntered(MouseEvent e) {
	}

	public void mouseExited(MouseEvent e) {
		if (!mIsResizing)
			setCursor(JTableWithRowNumbers.sDefaultCursor);
	}

	public void mouseMoved(MouseEvent e) {
		Point p = e.getPoint();

		int y = p.y + mTable.getScrollPane().getViewport().getViewRect().y;
		mDragStartRowHeight = mTable.getRowHeight();
		int tolerance = HiDPIHelper.scale(cResizeTolerance);
		int row = y / mDragStartRowHeight;
		int fracY = y % mDragStartRowHeight;

		mResizingStartRow = (fracY <= tolerance && row > 0) ? row - 1 : (fracY >= mDragStartRowHeight - tolerance) ? row : -1;
		setCursor(mResizingStartRow == -1 ? JTableWithRowNumbers.sDefaultCursor : JTableWithRowNumbers.sResizeCursor);

		if (mHighlightedRow != row) {
			mHighlightedRow = row;
			informListeners();
			repaint();
		}
	}

	public void mouseDragged(MouseEvent e) {
		if (mFullRowSelection) {
			int row = rowFromMouseY(e);
			if (mCurrentDragRow != row) {
				mCurrentDragRow = row;
				if (e.isShiftDown())
					mTable.getSelectionModel().addSelectionInterval(mDragStartRow, mCurrentDragRow);
				else if (e.isControlDown())
					mTable.getSelectionModel().removeSelectionInterval(mDragStartRow, mCurrentDragRow);
				else
					mTable.getSelectionModel().setSelectionInterval(mDragStartRow, mCurrentDragRow);
				repaint();
			}
		}
		if (mIsResizing) {
			int rowHeight = Math.min(getHeight()*4/5, Math.max(HiDPIHelper.scale(16), mDragStartRowHeight + e.getY() - mDragStartY));
			if (mResizingRowHeight != rowHeight) {
				mResizingRowHeight = rowHeight;
				mTable.setRowHeight(mResizingStartRow, mResizingRowHeight);
			}
			repaint();
		}
	}

	public void valueChanged(ListSelectionEvent e) {
		repaint();
	}

	public void addHighlightListener(HighlightListener l) {
		mHighlightListeners.add(l);
	}

	public void removeHighlightListener(HighlightListener l) {
		mHighlightListeners.remove(l);
	}

	public int getHighlightedRow() {
		return mHighlightedRow;
		}

	public void setHighlightedRow(int row, boolean informListeners) {
		mHighlightedRow = row;
		if (informListeners)
			informListeners();

		repaint();
		}

	private void informListeners() {
		for (HighlightListener l : mHighlightListeners)
			l.highlightedRowChanged(mHighlightedRow);
	}

	@Override public void paintComponent(Graphics g) {
		super.paintComponent(g);

		JScrollPane scrollPane = mTable.getScrollPane();
		if (scrollPane == null)
			return;

		Color bg = ColorHelper.darker(UIManager.getColor("Table.background"), 0.75f);
		Color fg = LookAndFeelHelper.isDarkLookAndFeel() ?
				ColorHelper.darker(UIManager.getColor("Table.foreground"), 0.8f)
				: ColorHelper.brighter(UIManager.getColor("Table.foreground"), 0.6f);
		Color lc = UIManager.getColor("Table.background");
		Color sbg = ColorHelper.darker(UIManager.getColor("Table.selectionBackground"), 0.75f);
		Dimension size = getSize();

		Rectangle viewRect = scrollPane.getViewport().getViewRect();
		int rowAtTop = mTable.rowAtPoint(new Point(0, viewRect.y));
		int rowAtBottom = mTable.rowAtPoint(new Point(0, viewRect.y + viewRect.height));
		int firstRow = Math.max(0, rowAtTop);
		int lastRow = (rowAtBottom == -1) ? mTable.getRowCount()-1 : rowAtBottom;
		int tableRowHeight = mTable.getRowHeight();

		int y = firstRow * tableRowHeight - viewRect.y;
		for (int row=firstRow; row<=lastRow; row++) {
			int rowHeight = (mIsResizing && row == mResizingStartRow) ? mResizingRowHeight : tableRowHeight;
			String rowHeader = Integer.toString(row+1);
			g.setColor(mTable.getSelectionModel().isSelectedIndex(row) ? sbg : bg);
			g.fillRect(0, y, size.width, rowHeight);
			g.setColor(fg);
			g.setFont(getFont());
			int x = (getWidth() - g.getFontMetrics().stringWidth(rowHeader)) / 2;
			g.drawString(rowHeader, x, y + (rowHeight + getFont().getSize()) / 2);
			if (row == mHighlightedRow) {
				g.setColor(Color.BLUE);
				for (int i=0; i<HiDPIHelper.scale(2); i++)
					g.drawRect(i, y+i, size.width-(2*i+1), rowHeight-(2*i+1));
			}
			else {
				g.setColor(lc);
				g.drawLine(0, y+rowHeight-1, size.width-1, y+rowHeight-1);
			}
			y += rowHeight;
		}
	}
}
