

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

package com.actelion.research.util.datamodel;

import java.awt.Point;
import java.io.*;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;


import com.actelion.research.util.Formatter;

public class ScorePoint extends Point {

	private static final long serialVersionUID = 24052013;
	
	private double score;
	
	public ScorePoint() {
		super(-1,-1);
	}

	public ScorePoint(Point p) {
		super(p);
	}
	
	public ScorePoint(Point p, double value) {
		super(p);
		score = value;
	}
	
	public ScorePoint(int x, int y, double value) {
		super(x,y);
		score = value;
	}

	public ScorePoint(int x, int y) {
		super(x,y);
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}
	
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		
		sb.append("x " + x);
		sb.append(" y " + y);
		sb.append(" " + Formatter.format3(score));
		
		return sb.toString();
	}


	public static List<ScorePoint> read(File fiTxt) throws IOException {

		List<ScorePoint> li = new ArrayList<ScorePoint>();

		BufferedReader br = new BufferedReader(new FileReader(fiTxt));

		String line = null;

		while((line = br.readLine())!= null) {

			String [] arr = line.split("\\t");

			int x = Integer.parseInt(arr[0].trim());
			int y = Integer.parseInt(arr[1].trim());
			double v = java.lang.Double.parseDouble(arr[2].trim());

			li.add(new ScorePoint(x,y,v));

		}

		br.close();

		return li;
	}

	public static int [] extractX(List<ScorePoint> li){

		int [] a = new int[li.size()];

		for (int i = 0; i < li.size(); i++) {
			a[i]=li.get(i).x;
		}

		return a;
	}

	public static Comparator<ScorePoint> getComparatorScore(){

		return new Comparator<ScorePoint>() {

			public int compare(ScorePoint sp1, ScorePoint sp2) {

				if(sp1.score>sp2.score){
					return 1;
				}else if(sp1.score<sp2.score){
					return -1;
				}

				return 0;
			}
		};
	}

	public static Comparator<ScorePoint> getComparatorX(){

		return new Comparator<ScorePoint>() {

			public int compare(ScorePoint sp1, ScorePoint sp2) {

				if(sp1.x>sp2.x){
					return 1;
				}else if(sp1.x<sp2.x){
					return -1;
				}

				return 0;
			}
		};
	}

	public static Comparator<ScorePoint> getComparatorY(){

		return new Comparator<ScorePoint>() {

			public int compare(ScorePoint sp1, ScorePoint sp2) {

				if(sp1.y>sp2.y){
					return 1;
				}else if(sp1.y<sp2.y){
					return -1;
				}

				return 0;
			}
		};
	}


}
