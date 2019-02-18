/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of ActelionMMFF94.
 * 
 * ActelionMMFF94 is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * ActelionMMFF94 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with ActelionMMFF94.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Paolo Tosco,Daniel Bergmann
 */

package mmff;

/**
 * Interface that must be implemented for array objects that can be
 * searched by the binary search function.
 */
public interface Searchable {
    /**
     * This should get an integer value given a column and row. The binary
     * search function only searches columns containing integers for a
     * value.
     *  @param row The row in the table.
     *  @param col The column in the row to return.
     *  @return The value at 'col' in 'row'.
     */
    public int get(int row, int col);

    /**
     * This function should return the total number of rows in a
     * searchable table. This is normally just the length of the array.
     *  @return The number of elements that can be searched.
     */
    public int length();
}
