package com.darkfuturestudios.ndktest

/**
 * Probably not necessary
 * Leads to garbage collector problems when used for bitmaps
 */
class ARGBColor(var a: Int = 0, var r: Int = 0, var g: Int = 0, var b: Int = 0) {

    /**
     * Stacks an ARBGColor onto this one by adding ARGB values
     */
    fun stack(color: ARGBColor) {
        a += color.a
        r += color.r
        g += color.g
        b += color.b
    }

    /**
     * Stacks an ARBGColor onto this one by adding ARGB values
     */
    fun stack(pa: Int, pr: Int, pg: Int, pb: Int) {
        a += pa
        r += pr
        g += pg
        b += pb
    }

}