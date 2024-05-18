package com.actelion.research.util;

import java.util.Calendar;
import java.util.GregorianCalendar;

public class DateAnalysis {
	private static final String TIME_FORMAT = "\\s\\d\\d:\\d\\d:\\d\\d";
	private static final int TIME_FORMAT_LENGTH = 9;

	private final int[] max,min,value;
	private int day,month,year,currentYear;
	private Calendar calendar;
	private boolean isTwoDigitYear;

	public DateAnalysis() {
		value = new int[3];
		max = new int[3];
		min = new int[3];
		min[0] = Integer.MAX_VALUE;
		min[1] = Integer.MAX_VALUE;
		min[2] = Integer.MAX_VALUE;
		day = -1;
		month = -1;
		year = -1;
		isTwoDigitYear = false;
		currentYear = new GregorianCalendar().get(Calendar.YEAR);
	}

	/**
	 * Checks, whether entry consists of three elements separated by one or more non-word characters
	 * (neither digit nor letter) that might be day, month and year.
	 * @param entry
	 * @return false if entry cannot be a date
	 */
	public boolean analyse(String entry) {
		if (endsWithTime(entry))
			entry = entry.substring(0, entry.length() - TIME_FORMAT_LENGTH);

		if (!interpreteValues(entry))
			return false;

		for (int i=0; i<3; i++) {
			if (min[i] > value[i])
				min[i] = value[i];
			if (max[i] < value[i])
				max[i] = value[i];
		}

		return true;
	}

	private boolean interpreteValues(String entry) {
		String[] item = entry.split("[^\\w]+");
		if (item.length != 3)
			return false;

		for (int i=0; i<3; i++) {
			value[i] = interprete(item[i], i);
			if (value[i] < 0 || value[i] > 3000)
				return false;
		}

		return true;
	}

	/**
	 * Call this only once after calling analyse() once or multiple times.
	 * If this method returns true, then you can use this object to interpret
	 * date strings of same format with getDateMillis(), e.g. the one passed
	 * with analyse().
	 * @return
	 */
	public boolean isConclusive() {
		for (int i=0; i<3; i++) {
			if (max[i] > 31 || min[i] == 0) {
				if (year != -1 && year != i)
					return false;
				year = i;
			}
		}

		// if there is no obvious year, we may have a 2-digit year as third of three entries
		if (year == -1) {
			if (max[2] < 100
					&& ((max[0] <= 12 && max[1] <= 31)
					|| (max[0] <= 31 && max[1] <= 12))) {
				year = 2;
				isTwoDigitYear = true;
			}
		}

		for (int i=0; i<3; i++) {
			if (i != year && (max[i] > 12 || (month != -1 && month != i))) {
				if (day != -1 && day != i)
					return false;
				day = i;
			}
		}

		for (int i=0; i<3; i++) {
			if (i != year && i != day) {
				if (month != -1 && month != i)
					return false;
				month = i;
			}
		}

		if (day == -1 || month == -1 || year == -1)
			return false;

		calendar = Calendar.getInstance();
		return true;
	}

	public long getDateMillis(String entry) {
		int hour = 12;
		int minute = 0;
		int second = 0;

		if (endsWithTime(entry)) {
			int index = entry.length()-TIME_FORMAT_LENGTH;
			hour   = Integer.parseInt(entry.substring(index + 1, index + 3));
			minute = Integer.parseInt(entry.substring(index + 4, index + 6));
			second = Integer.parseInt(entry.substring(index + 7, index + 9));
			entry = entry.substring(0, entry.length() - TIME_FORMAT_LENGTH);
		}

		if (!interpreteValues(entry))
			return -1;

		if (isTwoDigitYear)
			value[year] += (value[year] <= currentYear) ? 2000 : 1900;

		calendar.set(value[year], value[month]-1, value[day], hour, minute, second);
		return calendar.getTimeInMillis();
	}

	public float getDateFloat(String entry) {
		long millis = getDateMillis(entry);
		return millis == -1L ? Float.NaN : (float)((millis+43200000L)/86400000L);
	}

	private int interprete(String item, int index) {
		try {
			return Integer.parseInt(item);
		}
		catch (Exception e) {
			int value = interpreteMonth(item);
			if (value != -1) {
				if (month != -1 && month != index)
					return -1;	// month found in different positions
				month = index;
			}
			return value;
		}
	}

	public static int interpreteMonth(String item) {
		item = item.toLowerCase();
		if (item.startsWith("jan"))
			return 1;
		if (item.startsWith("feb"))
			return 2;
		if (item.startsWith("mar") || item.startsWith("mär") || item.startsWith("mae"))
			return 3;
		if (item.startsWith("apr"))
			return 4;
		if (item.startsWith("may") || item.startsWith("mai"))
			return 5;
		if (item.startsWith("jun"))
			return 6;
		if (item.startsWith("jul"))
			return 7;
		if (item.startsWith("aug"))
			return 8;
		if (item.startsWith("sep"))
			return 9;
		if (item.startsWith("oct") || item.startsWith("okt"))
			return 10;
		if (item.startsWith("nov"))
			return 11;
		if (item.startsWith("dec") || item.startsWith("dez"))
			return 12;
		return -1;
	}

	private boolean endsWithTime(String entry) {
		return entry != null
			&& entry.length() > TIME_FORMAT_LENGTH
			&& entry.substring(entry.length()-TIME_FORMAT_LENGTH).matches(TIME_FORMAT);
	}
}