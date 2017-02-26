// This work is licensed under the Creative Commons Attribution 3.0 Unported License.
// To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
// or send a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
// California, 94041, USA.
//
// Persistence Of Vision Ray Tracer ('POV-Ray') sample file.
//
// soft_object pattern example soft_object.pov.
//
// +w450 +h300 +a0.3

#version 3.72;
global_settings { assumed_gamma 1 }
#default { finish {ambient 0.005 diffuse 0.45} }
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
    ttf "timrom.ttf" "Water"
    0.05, 0.001
    translate <0.02,0.02,0>
}
#declare Grey_Green = srgbft <0.2863,0.3176,0.06667,0,0>;
#declare ObjectText = object {
    Text00
    pigment { color Grey_Green }
    translate <-1.3,0,0.3>
}
#declare VarTurb = 0.09;
#declare VarContainFuzz = (2*VarTurb);
#include "functions.inc"
#declare FnHardObj = function {
    pattern { hard_object { ObjectText }
        radius 0.04
        recursion_limit 10
        samples 22
        turbulence VarTurb octaves 3 lambda 3
    }
}
#declare Fn01 = function (x,y,z) {
    0.015-FnHardObj(x,y,z)
}
#declare Grenadier = srgbft <0.8353,0.2745,0,0,0>;
#declare Iso99 = isosurface {
    function { Fn01(x,y,z) }
    contained_by {
        box { min_extent(ObjectText)-VarContainFuzz,
              max_extent(ObjectText)+VarContainFuzz
        }
    }
    threshold 0
    accuracy 0.001
    max_gradient 83.3  // Less than max found for performance & OK look
    max_trace 1
    pigment { color Grenadier }
    translate <0,0,-0.6>
}

//---
camera { Camera00 }
light_source { Light00 }
object { CylinderX }
object { CylinderY }
object { CylinderZ }
object { Iso99 }
object { ObjectText }
