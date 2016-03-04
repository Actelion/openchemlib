/*
 * Created on May 28, 2004
 *
 */
package com.actelion.research.util;

import java.math.RoundingMode;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.StringTokenizer;

/**
 * 
 * @author freyssj
 */
public class Formatter {

	public static final DateFormat dateFormat = new SimpleDateFormat("dd.MM.yy");
	private static final DateFormat timeFormat = new SimpleDateFormat("HH:mm");

	private static final DateFormat[] dateTimeParsers = new SimpleDateFormat[] {		
		new SimpleDateFormat("dd.MM.yy HH:mm:ss"),
		new SimpleDateFormat("dd.MM.yy HH:mm"),
		new SimpleDateFormat("dd.MM.yy"),
		new SimpleDateFormat("MM.yy"),
		new SimpleDateFormat("yy")};
	private static final DateFormat[] dateTimeFormatters = new SimpleDateFormat[] {		
		new SimpleDateFormat("dd.MM.yyyy HH:mm:ss"),
		new SimpleDateFormat("dd.MM.yyyy HH:mm"),
		new SimpleDateFormat("dd.MM.yyyy"),
		new SimpleDateFormat("MM.yyyy"),
		new SimpleDateFormat("yyyy")};
//	private static final DateFormat dateTimeFormat2 = new SimpleDateFormat("dd.MM.yy HH:mm");
//	private static final DateFormat dateTimeFormat3 = new SimpleDateFormat("dd.MM.yy");

	private static final DecimalFormat df0 = new DecimalFormat("0");
	private static final DecimalFormat df1 = new DecimalFormat("0.0");
	private static final DecimalFormat df2 = new DecimalFormat("0.00");
	private static final DecimalFormat dfmax2 = new DecimalFormat("0.##");
	private static final DecimalFormat df3 = new DecimalFormat("0.000");
	private static final DecimalFormat dfmax3 = new DecimalFormat("0.###");
	private static final DecimalFormat df4 = new DecimalFormat("0.0000");
	private static final DecimalFormat df8 = new DecimalFormat("0.00000000");
	private static final DecimalFormat dfE = new DecimalFormat("0.00E0");
	
	public static final DecimalFormat dfI2 = new DecimalFormat("00");
	public static final DecimalFormat dfI3 = new DecimalFormat("000");

	static {
		//Decimal Format should round like Excel, i.e half up
		df0.setRoundingMode(RoundingMode.HALF_UP);
		df1.setRoundingMode(RoundingMode.HALF_UP);
		df2.setRoundingMode(RoundingMode.HALF_UP);
		dfmax2.setRoundingMode(RoundingMode.HALF_UP);
		df3.setRoundingMode(RoundingMode.HALF_UP);
		dfmax3.setRoundingMode(RoundingMode.HALF_UP);
		df4.setRoundingMode(RoundingMode.HALF_UP);
		df8.setRoundingMode(RoundingMode.HALF_UP);
		dfE.setRoundingMode(RoundingMode.HALF_UP);
	}
	
	public static final String format0(Double d) {
		if(d==null) return "";
		return df0.format(d);
	}

	public static final String format1(Double d) {
		if(d==null) return "";
		return df1.format(d);
	}

	public static final String format2(Double d) {
		if(d==null) return "";
		return df2.format(d);
	}

	public static final String formatMax2(Double d) {
		if(d==null) return "";
		return dfmax2.format(d);
	}

	public static final String format3(Double d) {
		if(d==null) return "";
		return df3.format(d);
	}

	public static final String formatMax3(Double d) {
		if(d==null) return "";
		return dfmax3.format(d);
	}
	
	public static final String format4(Double d) {
		if(d==null) return "";
		return df4.format(d);
	}

	public static final String format8(Double d) {
		if(d==null) return "";
		return df8.format(d);
	}
	
	public static final String formatE(Double d) {
		if(d==null) return "";
		return dfE.format(d);
	}

	public static final String format(Object value) {
		if(value==null) {
			return "";
		} else if(value instanceof Double) {
			return format3((Double) value);
		} else if(value instanceof Date) {
			return formatDateTime((Date) value);
		} else {
			return "" + value;
		}
	}


	

	public static final String formatDate(Date d) {
		return d == null ? "" : dateFormat.format(d);
	}

	public static final String formatTime(Date d) {
		return d == null ? "" : timeFormat.format(d);
	}

	public static final String formatDateTime(Date d) {
		if(d==null) return "";
		Calendar cal = Calendar.getInstance();
		cal.setTime(d);
		if(cal.get(Calendar.SECOND)!=0) return dateTimeFormatters[0].format(d);
		if(cal.get(Calendar.HOUR_OF_DAY)!=0 || cal.get(Calendar.MINUTE)!=0) return dateTimeFormatters[1].format(d);
		if(cal.get(Calendar.DAY_OF_MONTH)!=1) return dateTimeFormatters[2].format(d);
		if(cal.get(Calendar.MONTH)!=0) return dateTimeFormatters[3].format(d);
		return dateTimeFormatters[4].format(d);
	}
	
	public static final String formatDateTimeShort(Date d) {
		if(d==null) return "";
		return dateTimeParsers[2].format(d);
	}

	public static final Date parseDate(String s) {
		if(s==null || s.length()==0) return null;
		try {			
			return dateFormat.parse(s);
		} catch (Exception e) {
			return null;
		}
	}
	
	public static final Date parseDateTime(String s) {
		if(s==null || s.length()==0) return null;
		
		
//		if(true) {
//			s = s.replace('/', '.');
//			for (DateFormat dateFormat : dateTimeParsers) {
//				try {			
//					return dateFormat.parse(s);
//				} catch (Exception e) {
//				}			
//			}
//			return null;
//		}

		
		Calendar cal = Calendar.getInstance();
		
		
		int ndate = 0;
		int ntime = 0;
		int date[] = new int[3];
		int time[] = new int[3];
		boolean inDate = true;
		StringTokenizer st = new StringTokenizer(s.trim(), ":/. ", true);
		while(st.hasMoreTokens()) {
			String tok = st.nextToken();
			if(tok.equals("/") || tok.equals(".")) {
				if(inDate) {
					ndate++;
					if(ndate>2) return null; 
				} else {
					return null;
				}
			} else if(tok.equals(":")) {
				if(inDate) {
					return null;
				} else {
					ntime++;
					if(ntime>2) return null; 
				}				
			} else if(tok.equals(" ")) {
				if(inDate) {
					ndate++;
					inDate = false;
				} else {
					return null;
				}				
			} else {
				for (int i = 0; i < tok.length(); i++) {
					if(!Character.isDigit(tok.charAt(i))) return null;					
				}
				if(inDate) {
					date[ndate] = Integer.parseInt(tok);
				} else {
					time[ntime] = Integer.parseInt(tok);					
				}			
			}
		}
		if(inDate) {
			ndate++;
		} else {
			ntime++;
		}
//		System.out.println("Formatter.parseDateTime() "+ndate+">"+Arrays.toString(date)+" "+Arrays.toString(time));

		if(ndate==0) {
			cal.set(Calendar.DAY_OF_MONTH, 1);
			cal.set(Calendar.MONTH, 0);
			cal.set(Calendar.YEAR, 0);
		} else if(ndate==1) {
			cal.set(Calendar.DAY_OF_MONTH, 1);
			cal.set(Calendar.MONTH, 0);
			cal.set(Calendar.YEAR, date[0]<100? date[0]+2000: date[0]);
		} else if(ndate==2) {
//			if(date[0]<1 || date[0]>12) return null;
			cal.set(Calendar.DAY_OF_MONTH, 1);
			cal.set(Calendar.MONTH, date[0]-1);
			cal.set(Calendar.YEAR, date[1]<100? date[1]+2000: date[1]);
		} else if(ndate==3) {
//			if(date[0]<1 || date[0]>31) return null;
//			if(date[1]<1 || date[1]>12) return null;
			cal.set(Calendar.DAY_OF_MONTH, date[0]);
			cal.set(Calendar.MONTH, date[1]-1);
			cal.set(Calendar.YEAR, date[2]<100? date[2]+2000: date[2]);
		} else {
			return null; //impossible case
		}
		
		if(ntime==0) {
			cal.set(Calendar.HOUR_OF_DAY, 0);
			cal.set(Calendar.MINUTE, 0);
			cal.set(Calendar.SECOND, 0);
		} else if(ntime==1) {
			return null; //only hour or minutes?
		} else if(ntime==2) {
			cal.set(Calendar.HOUR_OF_DAY, time[0]);
			cal.set(Calendar.MINUTE, time[1]);
			cal.set(Calendar.SECOND, 0);
		} else if(ntime==3) {
			cal.set(Calendar.HOUR_OF_DAY, time[0]);
			cal.set(Calendar.MINUTE, time[1]);
			cal.set(Calendar.SECOND, time[2]);
		} else {
			return null; //impossible case
		}
		

		
		
		return cal.getTime();
	}

	public static final String cleanDate(String s) {
		if(s==null || s.length()==0) return s;
		try {
			return formatDate(parseDate(s));
		} catch (Exception e) {
			return null;
		}
			
	}
	
	public static final String cleanDateTime(String s) {
		if(s==null || s.length()==0) return s;
		Date d = parseDateTime(s);
		if(d==null) return "";
		return formatDateTime(d);
	}
	
	public static void main(String[] args) {
		System.out.println(parseDateTime("10.10.2013 12:23:20"));
		System.out.println(parseDateTime("10.10.2013 12:23"));
		System.out.println(parseDateTime("10.10.2013"));
		System.out.println(parseDateTime("10.2013"));
		System.out.println(parseDateTime("2013"));
		System.out.println();
		
		System.out.println(formatDateTime(parseDateTime("10.10.2013 12:23:20")));
		System.out.println(formatDateTime(parseDateTime("10.10.2013 12:23")));
		System.out.println(formatDateTime(parseDateTime("10.10.2013")));
		System.out.println();
		
		System.out.println("10.10.2013 12:23:20 -> "+cleanDateTime("10.10.2013 12:23:20"));
		System.out.println("10/1/2013 12:60 -> "+cleanDateTime("10.10.2013 12:60"));
		System.out.println("31.12.2013 23:60 -> "+cleanDateTime("31.12.2013 23:60"));
		System.out.println("10.12.2013 -> "+cleanDateTime("10.12.2013"));
		System.out.println("12.2013 -> "+cleanDateTime("12.2013"));
		System.out.println("2013 -> "+cleanDateTime("2013"));
		System.out.println("13 -> "+cleanDateTime("13"));
		System.out.println();
		System.out.println(cleanDateTime("10/10/2013"));
		System.out.println(cleanDateTime("toto"));

	}
}
