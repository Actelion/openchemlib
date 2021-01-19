package com.actelion.research.gui.table;

import com.actelion.research.gui.LookAndFeelHelper;
import com.actelion.research.util.ColorHelper;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.Vector;

public class JTableWithRowNumbers extends JTable implements TableModelListener {
    private static final long serialVersionUID = 0x20060906;

	private RowNumberTable	mRowNumberTable = null;
	private JScrollPane		mScrollPane = null;
	private Cursor			mResizeCursor,mDefaultCursor;
	private boolean			mIsResizing;
	private int				mRowHeaderMinWidth,mResizingRowY,mResizingRowHeight, mHeaderLineCount;

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
		mRowNumberTable.setSelectionModel(sm);
		}

    @SuppressWarnings("unchecked")
	public JTableWithRowNumbers(final Vector rowData, final Vector columnNames) {
		super(rowData, columnNames);
		initialize();
		}

	private void initialize() {
		mRowNumberTable = new RowNumberTable();
		mRowNumberTable.setRowHeight(getRowHeight());
		setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		mDefaultCursor = Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR);
		mResizeCursor = Cursor.getPredefinedCursor(Cursor.S_RESIZE_CURSOR);
		mHeaderLineCount = 1;
		}

	public JTable getRowNumberTable() {
		return mRowNumberTable;
		}

	@Override
	public void paintComponent(Graphics g) {
		super.paintComponent(g);

		if (mIsResizing) {
			int tableWidth = getWidth();
			g.setColor(LookAndFeelHelper.isDarkLookAndFeel() ? Color.WHITE : Color.BLACK);
			g.drawLine(0, mResizingRowY, tableWidth, mResizingRowY);
			g.drawLine(0, mResizingRowY+mResizingRowHeight, tableWidth, mResizingRowY+mResizingRowHeight);
			}
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
		if (mRowNumberTable != null) {
			mRowNumberTable.setFont(font);
			validateRowNumberColumnWidth();
			}
		}

	/**
	 * @param height
	 */
	@Override
	public void setRowHeight(int height) {
		super.setRowHeight(height);
		if (mRowNumberTable != null)
			mRowNumberTable.setRowHeight(height);
		}

	/**
	 * @param row
	 * @param height
	 */
	@Override
	public void setRowHeight(int row, int height) {
		super.setRowHeight(row, height);
		if (mRowNumberTable != null)
			mRowNumberTable.setRowHeight(row, height);
		}

	@Override
	public void addNotify() {
		super.addNotify();

		if (getParent().getParent() instanceof JScrollPane) {
			mScrollPane = (JScrollPane) getParent().getParent();
			mScrollPane.setRowHeaderView(mRowNumberTable);

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
			int maxWidth = mRowNumberTable.getFontMetrics(mRowNumberTable.getFont()).stringWidth(Integer.toString(getRowCount()));
			size.width = Math.max(mRowHeaderMinWidth, maxWidth);
			viewport.setPreferredSize(size);
			}
		}

	@Override
	public void tableChanged(TableModelEvent e) {
		super.tableChanged(e);

		if (mRowNumberTable != null) {
			validateRowNumberColumnWidth();
			mRowNumberTable.invalidate();
			mRowNumberTable.repaint();
			}
		}

	public void rowNumberClicked(int row) {	// overwrite this if you want specific effects
		}

	public void setSelectionModel(ListSelectionModel selectionModel) {
		super.setSelectionModel(selectionModel);
		if (mRowNumberTable != null)
			mRowNumberTable.setSelectionModel(selectionModel);
		}

	// Inner class used to display the row numbers on the left side of the table. Order
	// is considered important enough on this screen to re-enforce it with a visual.

	private class RowNumberTable extends JTable implements ListSelectionListener,MouseListener,MouseMotionListener {
        private static final long serialVersionUID = 0x20060906;

        private static final int cResizeTolerance = 2;

		private boolean mFullRowSelection;
		private int mDragStartY,mDragStartRow,mDragStartRowHeight;

		public RowNumberTable() {
			super();

			setAutoCreateColumnsFromModel(false);
			setModel(new RowNumberTableModel());
			JTableWithRowNumbers.this.setRowHeight(JTableWithRowNumbers.this.getRowHeight());

			setColumnSelectionAllowed(false);
			setRowSelectionAllowed(true);

			TableColumn column = new TableColumn();
			column.setResizable(false);
			column.setCellRenderer(new RowNumberRenderer());
			addColumn(column);

			addMouseListener(this);
			addMouseMotionListener(this);
			}

//		public boolean isFocusTraversable() {
//			return false;
//			}

		@Override
		public void setFont(Font font) {
			super.setFont(font.deriveFont(Font.BOLD));
			}

		public void mouseClicked(MouseEvent e) {
			if (getCursor() != mResizeCursor)
				rowNumberClicked(rowAtPoint(e.getPoint()));
			}

		public void mousePressed(MouseEvent e) {
			mIsResizing = (getCursor() == mResizeCursor);
			mFullRowSelection = !mIsResizing;
			if (mIsResizing) {
				Point p = e.getPoint();
				p.y -= cResizeTolerance;
				int row = rowAtPoint(p);
				mResizingRowY = getCellRect(row, 0, false).y-1;
				mDragStartY = e.getY();
				mDragStartRow = row;
				mDragStartRowHeight = mResizingRowHeight = getRowHeight(row);
				}
			}

		public void mouseReleased(MouseEvent e) {
			if (mIsResizing) {
				mIsResizing = false;
				if (mDragStartRowHeight != mResizingRowHeight)
					JTableWithRowNumbers.this.setRowHeight(mResizingRowHeight);
				}
			}

		public void mouseEntered(MouseEvent e) {
			}

		public void mouseExited(MouseEvent e) {
			mFullRowSelection = false;
			if (!mIsResizing)
				setCursor(mDefaultCursor);
			}

		public void mouseMoved(MouseEvent e) {
			int row = rowAtPoint(e.getPoint());
			Rectangle cellRect = getCellRect(row, 0, false);
			int y = e.getY() - cellRect.y;
			if ((y < cResizeTolerance && row != 0) || (y > cellRect.height-cResizeTolerance) && row != this.getRowCount()-1) {
				setCursor(mResizeCursor);
				this.setSelectionModel(new DefaultListSelectionModel());
				}
			else {
				setCursor(mDefaultCursor);
				this.setSelectionModel(JTableWithRowNumbers.this.getSelectionModel());
				}
			}

		public void mouseDragged(MouseEvent e) {
			if (mIsResizing) {
				int rowHeight = Math.max(16, mDragStartRowHeight + e.getY() - mDragStartY);
				if (mResizingRowHeight != rowHeight) {
					mResizingRowHeight = rowHeight;
					JTableWithRowNumbers.this.setRowHeight(mDragStartRow, mResizingRowHeight);
					}
				}
			}

		public void valueChanged(ListSelectionEvent e) {
			if (!e.getValueIsAdjusting()) {
				int[] selectedRow = getSelectedRows();

				if (selectedRow.length > 0) {
					JTableWithRowNumbers theTable = JTableWithRowNumbers.this;

					if (mFullRowSelection && theTable.getColumnCount() != 0)
						theTable.setColumnSelectionInterval(0, theTable.getColumnCount()-1);
					}
				}

            repaint();
			}
		}

 	private class RowNumberTableModel extends AbstractTableModel {
        private static final long serialVersionUID = 0x20060906;

        public int getColumnCount() {
			return 1;
			}

		public int getRowCount() {
			if (getModel() != null)
				return getModel().getRowCount();

			return 0;
			}

		public Object getValueAt(int r, int c) {
			return ""+(r+1);
			}
		}

	private class RowNumberRenderer extends JPanel implements TableCellRenderer {
        private static final long serialVersionUID = 0x20060906;

        private String mRowHeader;
		private Font mFont;

		public RowNumberRenderer() {
			super();
			}

		public void paintComponent(Graphics g) {
			super.paintComponent(g);

			Color bg = ColorHelper.darker(UIManager.getColor("Table.background"), 0.75f);
			Color fg = LookAndFeelHelper.isDarkLookAndFeel() ?
					  ColorHelper.darker(UIManager.getColor("Table.foreground"), 0.8f)
					: ColorHelper.brighter(UIManager.getColor("Table.foreground"), 0.6f);
			Dimension size = getSize();
			g.setColor(bg);
			g.fillRect(0, 0, size.width, size.height);
			g.setColor(fg);
			g.setFont(mFont);
			int x = (getWidth() - g.getFontMetrics().stringWidth(mRowHeader)) / 2;
			g.drawString(mRowHeader, x, (size.height+mFont.getSize())/2);
			}

		public Component getTableCellRendererComponent(JTable table, Object value,
								boolean isSelected, boolean hasFocus, int row, int col) {
			mRowHeader = Integer.toString(row+1);
			mFont = table.getFont();
			return this;
			}
		}
	}
