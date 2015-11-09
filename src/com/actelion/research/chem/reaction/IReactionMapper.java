/*
 * Project: DD_core
 * @(#)IReactionMapper.java
 *
 * Copyright (c) 1997- 2015
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All Rights Reserved.
 *
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 *
 * Author: Christian Rufener
 */

package com.actelion.research.chem.reaction;

import com.actelion.research.chem.SSSearcher;

public interface IReactionMapper
{
    Reaction matchReaction(Reaction r, SSSearcher sss);

}
