package com.actelion.research.chem.io.pdb.parser;

import java.util.List;

/**
 * ListInteger
 * <p>Copyright: Idorsia Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Idorsia Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * Created by korffmo1 on 22.03.18.
 */
public class ListInteger<T> {

    private int id;

    private List<T> li;

    public ListInteger(List<T> li, int id) {
        this.id = id;
        this.li = li;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public List<T> getLi() {
        return li;
    }

    public void setLi(List<T> li) {
        this.li = li;
    }
}
