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

package com.actelion.research.share.gui.editor.geom;

import java.awt.*;

/**
 * Project:
 * User: rufenec
 * Date: 11/24/2014
 * Time: 3:21 PM
 */
public interface ICursor
{

    int DEFAULT = java.awt.Cursor.DEFAULT_CURSOR;
    int SE_RESIZE = java.awt.Cursor.SW_RESIZE_CURSOR;
    int CROSSHAIR = java.awt.Cursor.CROSSHAIR_CURSOR;

    int MAX_DEFAULT_CURSOR =  Cursor.MOVE_CURSOR;
    int TOOL_CURSOR_BASE = Cursor.MOVE_CURSOR + 32;

    int TOOL_CHAINCURSOR = TOOL_CURSOR_BASE + 0;
    int TOOL_DELETECURSOR = TOOL_CURSOR_BASE + 1;
    int TOOL_HANDCURSOR = TOOL_CURSOR_BASE + 2;
    int TOOL_HANDPLUSCURSOR = TOOL_CURSOR_BASE + 3;
    int TOOL_FISTCURSOR = TOOL_CURSOR_BASE + 4;
    int TOOL_LASSOCURSOR = TOOL_CURSOR_BASE + 5;
    int TOOL_LASSOPLUSCURSOR = TOOL_CURSOR_BASE + 6;
    int TOOL_SELECTRECTCURSOR = TOOL_CURSOR_BASE + 7;
    int TOOL_SELECTRECTPLUSCURSOR = TOOL_CURSOR_BASE + 8;
    int TOOL_ZOOMCURSOR = TOOL_CURSOR_BASE + 9;

    int TOOL_POINTERCURSOR = Cursor.DEFAULT_CURSOR;
    int TOOL_CTEXTCURSOR = Cursor.TEXT_CURSOR;
}
