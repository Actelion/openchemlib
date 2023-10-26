package com.actelion.research.util.datamodel.table;

import com.actelion.research.calc.Matrix;
import com.actelion.research.util.StringFunctions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * TableModelText
 * Data model for an only text containing table
 * @author Modest von Korff
 * Apr 14, 2015 MvK Start implementation
 */
public class TableModelString {


	private static final String SEP_FIELD = "\t";
	private static final String SEP_LINE = "\n";

	private static final String SEP_FIELD_LATEX = " &\t";
	private static final String SEP_LINE_LATEX = " \\\\\n";

	private List<String> liColName;

	private List<String> liRowName;

	private List<List<String>> liliData;

	private HashMap<String, Integer> hmRowName_Index;
	private HashMap<String, Integer> hmColName_Index;


	/**
	 * 
	 */
	public TableModelString() {
		init();
	}

	public TableModelString(int rows, int cols) {
		init();
		init(rows, cols);
	}

	public TableModelString(List<String> liRowName, List<String> liColName) {
		init();
		init(liRowName.size(), liColName.size());
		for (int i = 0; i < liRowName.size(); i++) {
			setRowName(i, liRowName.get(i));
		}
		for (int i = 0; i < liColName.size(); i++) {
			setColName(i, liColName.get(i));
		}
	}

	public boolean containsColumn(String name){
		return (hmColName_Index.get(name)!=null)?true:false;
	}

	private void init(){
		liColName = new ArrayList<>();

		liRowName = new ArrayList<>();

		liliData = new ArrayList<>();

		hmRowName_Index = new HashMap<>();
		hmColName_Index = new HashMap<>();
	}

	private void init(int rows, int cols) {

		liliData.clear();

		liColName.clear();
		liRowName.clear();

		hmRowName_Index.clear();
		hmColName_Index.clear();

		for (int i = 0; i < rows; i++) {
			List<String> li = new ArrayList<>();

			for (int j = 0; j < cols; j++) {
				li.add(null);
			}
			liliData.add(li);

			liRowName.add("");
		}

		for (int i = 0; i < cols; i++) {
			liColName.add("");
		}
	}

	public String getColName(int col) {
		return liColName.get(col);
	}

	public void setColName(int col, String s) {
		liColName.set(col, s);
		hmColName_Index.put(s, col);

	}

	public String getRowName(int row) {
		return liRowName.get(row);
	}

	public void setRowName(int row, String s) {
		liRowName.set(row, s);
		hmRowName_Index.put(s, row);
	}

	public void set(int row, int col, String s) {
		liliData.get(row).set(col, s);
	}

	public void set(Matrix m, int digits) {

		int r = liRowName.size();
		int c = liColName.size();

		String sFormat="";
		if(digits > 0)
			sFormat += ".";

		for (int i = 0; i < digits; i++) {
			sFormat += "0";
		}

		DecimalFormat nf = new DecimalFormat(sFormat);

		for (int i = 0; i < r; i++) {

			for (int j = 0; j < c; j++) {
				set(i,j,nf.format(m.get(i,j)));
			}
		}
	}

	public void set(String rowName, String colName, String s) {
		int r = hmRowName_Index.get(rowName);
		int c = hmColName_Index.get(colName);
		liliData.get(r).set(c, s);
	}

	public String get(int row, int col) {
		return liliData.get(row).get(col);
	}
	
	
	public int getRows(){
		return liliData.size();
	}
	
	public int getCols(){
		
		if(getRows()==0) return 0;
		
		return liliData.get(0).size();
	}


	@Override
	public String toString() {

		final StringBuilder sb = new StringBuilder();

		int [] arrLengthMax = new int [liColName.size()];

		int lengthMaxRowName = 0;

		for (int i = 0; i < liRowName.size(); i++) {
			lengthMaxRowName = Math.max(lengthMaxRowName, liRowName.get(i).length());
		}

		for (int i = 0; i < liColName.size(); i++) {
			arrLengthMax[i]=Math.max(arrLengthMax[i],liColName.get(i).length());
			for (int j = 0; j < liRowName.size(); j++) {
				String s = get(j,i);
				if(s!=null){
					arrLengthMax[i]=Math.max(arrLengthMax[i],s.length());
				} else {
					System.out.println("TableModelString toString() no value for field: row ["+j+"] " + liRowName.get(j) + ", and col ["+i+"] " + liColName.get(i) + ".");
				}
			}
		}

		// Upper left corner is empty
		for (int i = 0; i < lengthMaxRowName; i++) {
			sb.append(" ");
		}

		sb.append(SEP_FIELD);

		// Column names
		for (int i = 0; i < liColName.size(); i++) {

			String name = liColName.get(i);

			int n = arrLengthMax[i] - name.length();

			for (int j = 0; j < n; j++) {
				sb.append(" ");
			}

			sb.append(name);

			if(i<liColName.size()-1) {
				sb.append(SEP_FIELD);
			}
		}

		sb.append(SEP_LINE);

		for (int i = 0; i < liRowName.size(); i++) {

			String name = liRowName.get(i);

			int n = lengthMaxRowName-name.length();

			for (int j = 0; j < n; j++) {
				sb.append(" ");
			}

			sb.append(name);
			sb.append(SEP_FIELD);
			for (int j = 0; j < liColName.size(); j++) {

				String sData = get(i,j);

				if(sData==null)
					sData="";

				n = arrLengthMax[j]-sData.length();

				for (int k = 0; k < n; k++) {
					sb.append(" ");
				}

				sb.append(sData);

				if(j<liColName.size()-1) {
					sb.append(SEP_FIELD);
				}
			}
			if(i<liRowName.size()-1) {
				sb.append(SEP_LINE);
			}
		}

		return sb.toString();
	}

	/**
	 * \begin{table}[]
	 * \begin{tabular}{lllllllll}
	 *                             & 0.20   & 0.30   & 0.40   & 0.50   & 0.60   & 0.70   & 0.80   & 0.90   \\
	 * Gaussian process regression & 1.3550 & 1.3059 & 1.3818 & 1.4046 & 1.5150 & 1.4804 & 1.4579 & 1.6223 \\
	 * KNN regression              & 1.6776 & 1.7294 & 1.8899 & 2.1573 & 2.3484 & 2.6003 & 2.6931 & 2.8232 \\
	 * Median                      & 1.6730 & 1.7421 & 1.9217 & 2.3378 & 2.7606 & 3.2241 & 3.6768 & 4.5243 \\
	 * PLS                         & 1.4264 & 1.2870 & 1.2131 & 1.2989 & 1.3866 & 1.3830 & 1.4414 & 1.1730 \\
	 * PLS Power                   & 1.4886 & 1.3616 & 1.4954 & 1.7442 & 2.2339 & 2.3904 & 2.1669 & 2.6130 \\
	 * Random Forest regression    & 1.6746 & 1.6848 & 1.8421 & 2.0106 & 2.2479 & 2.3896 & 2.5420 & 2.6523 \\
	 * SVM regression              & 1.5580 & 1.4634 & 1.5879 & 1.6687 & 1.8423 & 1.9769 & 1.8803 & 1.9618
	 * \end{tabular}
	 * \end{table}
	 */


	public String toStringLaTex() {

		final StringBuilder sb = new StringBuilder();

		sb.append("\\begin{table}[]");
		sb.append(SEP_LINE);

		String align = "l";
		for (int i = 0; i < liColName.size(); i++) {
			align += "l";
		}

		sb.append("\\begin{tabular}{"+align+"}");
		sb.append(SEP_LINE);

		int [] arrLengthMax = new int [liColName.size()];

		int lengthMaxRowName = 0;

		for (int i = 0; i < liRowName.size(); i++) {
			lengthMaxRowName = Math.max(lengthMaxRowName, liRowName.get(i).length());
		}

		for (int i = 0; i < liColName.size(); i++) {

			arrLengthMax[i]=Math.max(arrLengthMax[i],liColName.get(i).length());

			for (int j = 0; j < liRowName.size(); j++) {
				String s = get(j,i);
				arrLengthMax[i]=Math.max(arrLengthMax[i],s.length());
			}
		}

		// Upper left corner is empty
		for (int i = 0; i < lengthMaxRowName; i++) {
			sb.append(" ");
		}

		sb.append(SEP_FIELD_LATEX);

		// Column names
		for (int i = 0; i < liColName.size(); i++) {

			String name = liColName.get(i);

			int n = arrLengthMax[i] - name.length();

			for (int j = 0; j < n; j++) {
				sb.append(" ");
			}

			sb.append(name);

			if(i<liColName.size()-1) {
				sb.append(SEP_FIELD_LATEX);
			}
		}

		sb.append(SEP_LINE_LATEX);

		for (int i = 0; i < liRowName.size(); i++) {

			String name = liRowName.get(i);

			int n = lengthMaxRowName-name.length();

			for (int j = 0; j < n; j++) {
				sb.append(" ");
			}

			sb.append(name);
			sb.append(SEP_FIELD_LATEX);
			for (int j = 0; j < liColName.size(); j++) {

				String sData = get(i,j);
				n = arrLengthMax[j]-sData.length();

				for (int k = 0; k < n; k++) {
					sb.append(" ");
				}

				sb.append(sData);

				if(j<liColName.size()-1) {
					sb.append(SEP_FIELD_LATEX);
				}
			}
			if(i<liRowName.size()-1) {
				sb.append(SEP_LINE_LATEX);
			}
		}

		sb.append(SEP_LINE);
		sb.append("\\end{tabular}");
		sb.append(SEP_LINE);
		sb.append("\\end{table}");

		return sb.toString();
	}

	public void write(File fiTxt, String rowName) throws IOException {

		BufferedWriter bw = new BufferedWriter(new FileWriter(fiTxt));

		bw.write(rowName);
		bw.write(SEP_FIELD);
		for (int i = 0; i < liColName.size(); i++) {
			bw.write(liColName.get(i));
			if(i<liColName.size()-1){
				bw.write(SEP_FIELD);
			}
		}
		bw.write(SEP_LINE);

		for (int i = 0; i < liRowName.size(); i++) {

			bw.write(liRowName.get(i));
			bw.write(SEP_FIELD);

			List<String> liData = liliData.get(i);

			if(liData!=null) {
				for (int j = 0; j < liData.size(); j++) {

					String str = liData.get(j);
					if(str!=null)
						bw.write(str);
					if (i < liData.size() - 1) {
						bw.write(SEP_FIELD);
					}
				}
			}

			if (i < liRowName.size() - 1) {
				bw.write(SEP_LINE);
			}
		}

		bw.close();

	}

}
