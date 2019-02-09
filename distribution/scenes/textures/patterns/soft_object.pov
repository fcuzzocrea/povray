// This work is licensed under the Creative Commons Attribution 3.0 Unported License.
// To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
// or send a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
// California, 94041, USA.
//
// Persistence Of Vision Ray Tracer ('POV-Ray') sample file.
//
// soft_object pattern example soft_object.pov.
//
// +w450 +h300 -a0.3     // +a0.3 takes approx 5x longer

#version 3.8;
global_settings { assumed_gamma 1 }
#declare Grey50 = srgbft <0.5,0.5,0.5,0,0>;
background { color Grey50 }
#declare Camera00 = camera {
    perspective
    location <2.5,2.5,-2.501>
    sky <0,1,0>
    angle 35
    right x*(image_width/image_height)
    look_at <0,0.25,0>
}
#declare White = srgbft <1,1,1,0,0>;
#declare Light00 = light_source {
    <50,150,-250>, White
}
#declare Red = srgbft <1,0,0,0,0>;
#declare CylinderX = cylinder {
    <-0.6,0,0>, <0.6,0,0>, 0.01
    pigment { color Red }
}
#declare Green = srgbft <0,1,0,0,0>;
#declare CylinderY = cylinder {
    <0,-0.6,0>, <0,0.6,0>, 0.01
    pigment { color Green }
}
#declare Blue = srgbft <0,0,1,0,0>;
#declare CylinderZ = cylinder {
    <0,0,-0.6>, <0,0,0.6>, 0.01
    pigment { color Blue }
}
#declare Text00 = text {
    ttf "timrom.ttf" "Water" 0.05, 0.001 translate <0.02,0.02,0>
}
#declare Grey_Green = srgbft <0.2863,0.3176,0.06667,0,0>;
#declare TransformTextPos = transform { translate <-1.3,0,0.3> };
#declare TransformIsoOffsetZ = transform { translate <0,0,-0.6> };
#declare ObjectText = object {
    Text00
    pigment { color Grey_Green }
    transform { TransformTextPos }
}
#include "functions.inc"
#declare VarSpacing = 0.009;
#declare VarTurb = 0.03;
#declare VarContainFuzz = (5*VarSpacing)+(2*VarTurb);
#declare FnSoftObj = function {
    pattern { soft_object { ObjectText }
              spacing VarSpacing
              strength 1.00
              warp { turbulence VarTurb lambda 5 }
    }
}
#declare FnSoftObj01 = function (x,y,z) {
    0.015-FnSoftObj(x,y,z)
}
#declare Grenadier = srgbft <0.8353,0.2745,0,0,0>;
#declare IsoText = isosurface {
    function { FnSoftObj01(x,y,z) }
    contained_by {
        box { min_extent(ObjectText)-VarContainFuzz,
              max_extent(ObjectText)+VarContainFuzz
        }
    }
    threshold 0
    accuracy 0.001
    max_gradient 3.3 // Less than max found for performance & OK look
    max_trace 1
    pigment { color Grenadier }
    transform { TransformIsoOffsetZ }
}

//--- Version of text by character
#declare VarCharWidth = 0.4; // Guessing average
#declare Text_W = text { ttf "timrom.ttf" "W" 0.05, 0 }
#declare Text_a = text { ttf "timrom.ttf" "a" 0.05, 0 }
#declare Text_t = text { ttf "timrom.ttf" "t" 0.05, 0 }
#declare Text_e = text { ttf "timrom.ttf" "e" 0.05, 0 }
#declare Text_r = text { ttf "timrom.ttf" "r" 0.05, 0 }
#declare ObjectText_W = object {
    Text_W
    transform { TransformTextPos }
}
#declare ObjectText_a = object {
    Text_a
    translate <VarCharWidth*2.2,0,0>
    transform { TransformTextPos }
}
#declare ObjectText_t = object {
    Text_t
    translate <VarCharWidth*3.4,0,0>
    transform { TransformTextPos }
}
#declare ObjectText_e = object {
    Text_e
    translate <VarCharWidth*4.3,0,0>
    transform { TransformTextPos }
}
#declare ObjectText_r = object {
    Text_r
    translate <VarCharWidth*5.4,0,0>
    transform { TransformTextPos }
}
#declare UnionTextChars = union {
    object { ObjectText_W }
    object { ObjectText_a }
    object { ObjectText_t }
    object { ObjectText_e }
    object { ObjectText_r }
    pigment { color Grey_Green }
}
//--- W   (Yes, could be a macro)
#declare FnSoftObj_W = function {
    pattern { soft_object { ObjectText_W }
              spacing VarSpacing
              strength 1.00
              warp { turbulence VarTurb lambda 5 }
    }
}
#declare FnSoftObj01_W = function (x,y,z) {
    0.015-FnSoftObj_W(x,y,z)
}
#declare IsoText_W = isosurface {
    function { FnSoftObj01_W(x,y,z) }
    contained_by {
        box { min_extent(ObjectText_W)-VarContainFuzz,
              max_extent(ObjectText_W)+VarContainFuzz
        }
    }
    threshold 0
    accuracy 0.001
    max_gradient 3.3
    max_trace 1
    pigment { color Grenadier }
    transform { TransformIsoOffsetZ }
}
//--- a
#declare FnSoftObj_a = function {
    pattern { soft_object { ObjectText_a }
              spacing VarSpacing
              strength 1.00
              warp { turbulence VarTurb lambda 5 }
    }
}
#declare FnSoftObj01_a = function (x,y,z) {
    0.015-FnSoftObj_a(x,y,z)
}
#declare IsoText_a = isosurface {
    function { FnSoftObj01_a(x,y,z) }
    contained_by {
        box { min_extent(ObjectText_a)-VarContainFuzz,
              max_extent(ObjectText_a)+VarContainFuzz
        }
    }
    threshold 0
    accuracy 0.001
    max_gradient 3.3
    max_trace 1
    pigment { color Grenadier }
    transform { TransformIsoOffsetZ }
}
//--- t
#declare FnSoftObj_t = function {
    pattern { soft_object { ObjectText_t }
              spacing VarSpacing
              strength 1.00
              warp { turbulence VarTurb lambda 5 }
    }
}
#declare FnSoftObj01_t = function (x,y,z) {
    0.015-FnSoftObj_t(x,y,z)
}
#declare IsoText_t = isosurface {
    function { FnSoftObj01_t(x,y,z) }
    contained_by {
        box { min_extent(ObjectText_t)-VarContainFuzz,
              max_extent(ObjectText_t)+VarContainFuzz
        }
    }
    threshold 0
    accuracy 0.001
    max_gradient 3.3
    max_trace 1
    pigment { color Grenadier }
    transform { TransformIsoOffsetZ }
}
//--- e
#declare FnSoftObj_e = function {
    pattern { soft_object { ObjectText_e }
              spacing VarSpacing
              strength 1.00
              warp { turbulence VarTurb lambda 5 }
    }
}
#declare FnSoftObj01_e = function (x,y,z) {
    0.015-FnSoftObj_e(x,y,z)
}
#declare IsoText_e = isosurface {
    function { FnSoftObj01_e(x,y,z) }
    contained_by {
        box { min_extent(ObjectText_e)-VarContainFuzz,
              max_extent(ObjectText_e)+VarContainFuzz
        }
    }
    threshold 0
    accuracy 0.001
    max_gradient 3.3
    max_trace 1
    pigment { color Grenadier }
    transform { TransformIsoOffsetZ }
}
//--- r
#declare FnSoftObj_r = function {
    pattern { soft_object { ObjectText_r }
              spacing VarSpacing
              strength 1.00
              warp { turbulence VarTurb lambda 5 }
    }
}
#declare FnSoftObj01_r = function (x,y,z) {
    0.015-FnSoftObj_r(x,y,z)
}
#declare IsoText_r = isosurface {
    function { FnSoftObj01_r(x,y,z) }
    contained_by {
        box { min_extent(ObjectText_r)-VarContainFuzz,
              max_extent(ObjectText_r)+VarContainFuzz
        }
    }
    threshold 0
    accuracy 0.001
    max_gradient 3.3
    max_trace 1
    pigment { color Grenadier }
    transform { TransformIsoOffsetZ }
}

//---
camera { Camera00 }
light_source { Light00 }
object { CylinderX }
object { CylinderY }
object { CylinderZ }
// --- Simpler set 4x slower than by char method below
// object { ObjectText }
// object { IsoText }
// --- Messier by char set up is ~4x faster than single text object.
object { UnionTextChars }
object { IsoText_W }
object { IsoText_a }
object { IsoText_t }
object { IsoText_e }
object { IsoText_r }

