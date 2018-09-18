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

package com.actelion.research.util;

import java.awt.*;
import java.util.*;
import java.util.List;

public class PointUtils {


    public static List<Point> getPointsBetween(Point p1, Point p2){

        List<Point> li = new ArrayList<Point>();

        int xStart = Math.min(p1.x, p2.x);
        int yStart = Math.min(p1.y, p2.y);

        int xEnd = Math.max(p1.x, p2.x) + 1;
        int yEnd = Math.max(p1.y, p2.y) + 1;

        for (int i = xStart; i < xEnd; i++) {

            for (int j = yStart; j < yEnd; j++) {

                Point p = new Point(i,j);

                li.add(p);
            }

        }

        return li;

    }

    public static Comparator<Point> getComparatorX(){

        return new Comparator<Point>() {
            @Override
            public int compare(Point p1, Point p2) {

                int cmp = 0;

                if(p1.x > p2.x){
                    cmp=1;
                }else if(p1.x < p2.x){
                    cmp=-1;
                }

                return cmp;
            }
        };
    }


}
