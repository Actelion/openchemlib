package com.actelion.research.util;


import java.util.ArrayList;
import java.util.Date;
import java.util.List;

/**
 * <p>TimeDelta: </p>
 * <p>Description:  </p>
 * <p>Company: Actelion Ltd. </p>
 * @author Modest von Korff
 * @version 1.0
 * 22.04.2005 Start implementation
 */
public class TimeDelta {



	 public static final int PRECISION_MS = 0;
	
	 public static final int PRECISION_SECONDS = 1;
	 
	 public static final int PRECISION_MINUTES = 2;
	 
	 public static final int PRECISION_HOURS = 3;

    // So many nano seconds has a milli second
    public static final long NANO_MS = 1000 * 1000;

	 // So many nano seconds has a second
     public static final long NANO_SECOND = NANO_MS * 1000;

     public static final long MS_SECOND = 1000;

     public static final long MS_MINUTE = MS_SECOND * 60;

     public static final long MS_HOUR = MS_MINUTE * 60;
     
     public static final long MS_DAY = MS_HOUR * 24;
     
     /**
      * A month has 30 days here.
      */
     public static final long MS_MONTH = MS_DAY * 30;
     
     public static final long MS_YEAR = MS_DAY * 365;

     // For logging medium long running programs.
     public static final long [] ARR = {MS_MINUTE, MS_MINUTE*5, MS_MINUTE*30, MS_HOUR};



    // Time in milliseconds
    private long time;
    
    private int precision;

    public TimeDelta(long l) {
        time = l;
        precision = PRECISION_SECONDS;
    }
    public TimeDelta(long l, int precision) {
        time = l;
        this.precision = precision;
    }

    public TimeDelta(int precision) {
    	this.precision = precision;
    }
    
    
    public static void main(String [] args) {
        
        List<Long> li = new ArrayList<Long>();
        li.add(0l);
        li.add(1l);
        li.add(1120l);
        li.add(60000l);
        li.add(101010101l);
        li.add(100000000001l);
    	
        for (long d : li) {
        	System.out.println(toString(d,0));	
		}
    }

    public String format(long milliseconds){
    	return toString(milliseconds, precision);
    }
    
    public String format(Date date){
    	return toString(date.getTime(), precision);
    }
    
    public void setPrecision(int prec) {
    	precision = prec;
    }
    public String toString() {
        return toString(time,precision);
    }
    
    /**
     * Precision of the output depends on the input value.
     * @param milliseconds
     * @return
     */
    public static String toString(long milliseconds) {
    	
    	int precision = PRECISION_MS;
    	
    	if(milliseconds > MS_YEAR) {
    		precision = PRECISION_HOURS;
    	} else if(milliseconds > MS_DAY) {
    		precision = PRECISION_MINUTES;
    	} else if(milliseconds > MS_HOUR) {
    		precision = PRECISION_SECONDS;
    	}
    	    	
    	return toString(milliseconds, precision);
    }
    
    public static String toString(long milliseconds, int precision) {
        String str = "";

        if(milliseconds==0) {
            return milliseconds + " Millisec";
        }

        long millisecRemaining = milliseconds;

        if(millisecRemaining / MS_YEAR >= 1) {
        	String s = "";
        	if(millisecRemaining / MS_YEAR == 1){
        		s = "1 Year ";
        	} else {
        		s = (millisecRemaining / MS_YEAR) + " Years ";
        	}
            str += s;
            millisecRemaining = millisecRemaining - (MS_YEAR * (millisecRemaining / MS_YEAR));
        }
        
        if(millisecRemaining / MS_DAY >= 1) {
        	String s = "";
        	if(millisecRemaining / MS_DAY == 1){
        		s = "1 Day ";
        	} else {
        		s = (millisecRemaining / MS_DAY) + " Days ";
        	}
        	
            str += s;
            millisecRemaining = millisecRemaining - (MS_DAY * (millisecRemaining / MS_DAY));
        }
        
        if(millisecRemaining / MS_HOUR >= 1) {
        	String s = "";
        	if(millisecRemaining / MS_HOUR == 1){
        		s = "1 Hour ";
        	} else {
        		s = (millisecRemaining / MS_HOUR) + " Hours ";
        	}
        	
            str += s;
            millisecRemaining = millisecRemaining - (MS_HOUR * (millisecRemaining / MS_HOUR));
        }
        
        if(precision == PRECISION_HOURS)
        	return str.trim();
        
        if(millisecRemaining / MS_MINUTE >= 1) {
        	String s = "";
        	if(millisecRemaining / MS_MINUTE == 1){
        		s = "1 Minute ";
        	} else {
        		s = (millisecRemaining / MS_MINUTE) + " Minutes ";
        	}
        	
            str += s;
            millisecRemaining = millisecRemaining - (MS_MINUTE * (millisecRemaining / MS_MINUTE));
        }
        
        if(precision == PRECISION_MINUTES)
        	return str.trim();
        if(millisecRemaining / MS_SECOND >= 1) {
        	String s = "";
        	if(millisecRemaining / MS_SECOND == 1){
        		s = "1 Second";
        	} else {
        		s = (millisecRemaining / MS_SECOND) + " Seconds";
        	}
            
            str += s;
            millisecRemaining = millisecRemaining - (MS_SECOND * (millisecRemaining / MS_SECOND));
        }
        
        if(precision == PRECISION_SECONDS)
        	return str.trim();
        
        if(millisecRemaining > 0) {
            String s = " " + millisecRemaining + " Millisec";
            str += s;
        }

        return str.trim();
    }

    public static boolean isOlderThanOneWeek(long t){
    	
    	return isOlderThanDays(t, 7);
    }
    
    /**
     * 
     * @param t
     * @param days
     * @return true if the date t is more than <code>days</code> days before the current date.
     */
    public static boolean isOlderThanDays(long t, int days){
    	boolean older = false;
    	long diff = new Date().getTime() - t;
    	if(diff > MS_DAY * days){
    		older = true;
    	}
    	return older;
    }

    public static boolean isOlderThanMilliseconds(long t, long ms){
    	boolean older = false;
    	long diff = new Date().getTime() - t;
    	if(diff > ms){
    		older = true;
    	}
    	return older;
    }

    public static boolean isOlderThanHours(long t, int hours){

    	boolean older = false;

    	long diff = new Date().getTime() - t;

    	if(diff > MS_HOUR * hours){
    		older = true;
    	}

    	return older;
    }


};