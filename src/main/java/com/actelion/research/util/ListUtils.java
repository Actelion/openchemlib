package com.actelion.research.util;

import java.util.ArrayList;
import java.util.List;

/**
 * ListUtils
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 25.10.18.
 */
public class ListUtils {

    /**
     * Order dependent.
     * @param li1
     * @param li2
     * @return
     */
    public static boolean equals(List<String> li1, List<String> li2) {

        if(li1.size()!=li2.size()){
            return false;
        }

        boolean eq = true;
        for (int i = 0; i < li1.size(); i++) {

            String s1 = li1.get(i);
            String s2 = li2.get(i);

            if(!s1.equals(s2)){
                eq=false;
                break;
            }
        }

        return eq;
    }


    public static String toString(List<String> li){
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < li.size(); i++) {
            sb.append(li.get(i));
            if(i<li.size()-1){
                sb.append(ConstantsDWAR.SEP_VALUE);
            }
        }

        return sb.toString();
    }

    public static String toStringInteger(List<Integer> li){
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < li.size(); i++) {
            sb.append(li.get(i));
            if(i<li.size()-1){
                sb.append(ConstantsDWAR.SEP_VALUE);
            }
        }

        return sb.toString();
    }


    public static  List<Integer> createIndexList(int n){
        List<Integer> liIndex = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            liIndex.add(i);
        }
        return liIndex;
    }


}
