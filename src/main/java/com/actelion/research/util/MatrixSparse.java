package com.actelion.research.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.StringTokenizer;

import com.actelion.research.calc.Matrix;

/**
 * 
 * 
 * MatrixSparse
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * 2008 MvK: Start implementation
 */
public class MatrixSparse {
	private static final int BUFFER_SIZE = 10000;
	private static final int READ_BUFFER_SIZE = 10^7;
	
	private LinkedHashMap<Integer, LinkedHashMap<Integer,Integer>> hmIndex_hmIndex;
	
	private double [] data;
	
	private int position;
	
	private int rows;
	
	private int cols;
	
	public MatrixSparse(File fi) throws Exception{
		init();
		read(fi);
	}
	
	public MatrixSparse(File fi, int rows) throws Exception{
		init();
		read(fi, rows);
	}

	public MatrixSparse(int rows, int cols){
		init();
		
		this.rows=rows;
		
		this.cols=cols;
	}
	
	
	public void init(){
		hmIndex_hmIndex = new LinkedHashMap<Integer, LinkedHashMap<Integer,Integer>>();
		data = new double [BUFFER_SIZE];
		position = -1;
	}
	
	public void set(int row, int col, double v){
		
		if(row >= rows){
			String m = "Row index " + row + " out of bounds " + rows + ".";
			throw new IndexOutOfBoundsException(m);
		} else if(col >= cols){
			String m = "Col index " + col + " out of bounds " + cols + ".";
			throw new IndexOutOfBoundsException(m);
		} 
		
		position++;
		
		if(position == data.length){
			resizeArray(position+1);
		}
		
		data[position]=v;
		LinkedHashMap<Integer,Integer> lhmRow = hmIndex_hmIndex.get(row);
		if(lhmRow == null){
			lhmRow = new LinkedHashMap<Integer,Integer>();
			hmIndex_hmIndex.put(row, lhmRow);
		}
		
		lhmRow.put(col, position);
		
		
	}
	
	/**
	 * 
	 * @param row
	 * @param col
	 * @return NaN if the value was not set.
	 */
	public double get(int row, int col){
		
		if(hmIndex_hmIndex.containsKey(row)){
			if(hmIndex_hmIndex.get(row).containsKey(col)){
				int pos = hmIndex_hmIndex.get(row).get(col);
				return data[pos];
			}
		}
		
		return Double.NaN;
	}
	
	public double [] getArray(){
		return data;
	}
	public Matrix getMatrixK(){
		Matrix ma = new Matrix(rows(), cols());
		
		for (int i = 0; i < rows(); i++) {
			for (int j = 0; j < cols(); j++) {
				ma.set(i,j, get(i,j));
			}
		}
		return ma;
	}
	
	
	public int cols(){
		return cols;
	}
	public int rows(){
		return rows;
	}
	
	private void read(File file) throws IOException {
		
		FileInputStream iStr = new FileInputStream(file);
		FileChannel fc = iStr.getChannel();
		
		
		ByteBuffer buf = ByteBuffer.allocate(READ_BUFFER_SIZE);
		
		StringBuilder sb = new StringBuilder();
		
		boolean bEOF=false;
		
		while (!bEOF){
			buf.position(0);

			int size = fc.read(buf);
			if(size==-1){
				break;
			}
			
			for (int i = 0; i < size; i++) {
			
				int c = buf.get(i);
				if(c==-1) {
					bEOF=true;
				} else {
					if(c=='\n') {
						read(sb.toString());
						sb.delete(0, sb.length());
					} else {
						sb.append((char)c);
					}
				}
			}
		}
		
		fc.close();

	}
	
	private void read(File file, int rows) throws IOException {
		
		FileInputStream iStr = new FileInputStream(file);
		FileChannel fc = iStr.getChannel();
		
		
		ByteBuffer buf = ByteBuffer.allocate(READ_BUFFER_SIZE);
		
		StringBuilder sb = new StringBuilder();
		
		boolean bEOF=false;
		
		int cc=0;
		while (!bEOF){
			buf.position(0);

			int size = fc.read(buf);
			if(size==-1){
				break;
			}
			
			for (int i = 0; i < size; i++) {
			
				int c = buf.get(i);
				if(c==-1) {
					bEOF=true;
				} else {
					if(c=='\n') {
						read(sb.toString());
						sb.delete(0, sb.length());
						cc++;
						if(cc==rows){
							bEOF = true;
							break;
						}
					} else {
						sb.append((char)c);
					}
				}
			}
		}
		
		fc.close();

	}

	private void read(String line){
		StringTokenizer st = new StringTokenizer(line);
		int row = Integer.parseInt(st.nextToken());
		int col = Integer.parseInt(st.nextToken());
		double v = Double.parseDouble(st.nextToken());
		
		if(row>=rows)
			rows = row + 1;
		
		if(col>=cols)
			cols = col + 1;
		
		set(row, col, v);
		
	}
	
	private void resizeArray(int newsize){

		if(data.length < newsize){
			int si = data.length;
			while(si < newsize) {
				si += BUFFER_SIZE;  
			}
			
			double [] dataTmp = new double [si];
			
			int length = Math.min(position, newsize);
			
			System.arraycopy(data, 0, dataTmp, 0, length);
			
			data = dataTmp;
		}
	}
	
	public String toStringRow(int row, NumberFormat nf) {
		StringBuilder sb = new StringBuilder();
		
		
		
		if(hmIndex_hmIndex.get(row)!=null){
			List<Integer> liCol = new ArrayList<Integer>(hmIndex_hmIndex.get(row).keySet());
			for (Integer col : liCol) {
				
				double v = data[hmIndex_hmIndex.get(row).get(col)];
				String s = row + " " + col + " " + nf.format(v) + "\n"; 
				sb.append(s);
			}
		}
		
		return sb.toString();
	}

	
	public String toString(NumberFormat nf) {
		StringBuilder sb = new StringBuilder();
		
		List<Integer> liRow = new ArrayList<Integer>(hmIndex_hmIndex.keySet()); 

		for (Integer row : liRow) {
			List<Integer> liCol = new ArrayList<Integer>(hmIndex_hmIndex.get(row).keySet());
			for (Integer col : liCol) {
				
				double v = data[hmIndex_hmIndex.get(row).get(col)];
				String s = row + " " + col + " " + nf.format(v) + "\n"; 
				sb.append(s);
			}
		}
		return sb.toString();
	}
	
	public String toString() {
		NumberFormat nf = new DecimalFormat("0.############");
		return toString(nf);
	}

}
