package com.actelion.research.util.datamodel;

/**
 * IModelClonable
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 03.01.19.
 */
public interface IModelCloneable<T> {

    /**
     *
     * @return a deep copy of this is needed.
     */
    IModelCloneable<T> getDeepClone();

}
