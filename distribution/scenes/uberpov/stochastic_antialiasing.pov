// This work is licensed under the Creative Commons Attribution 3.0 Unported License.
// To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
// or send a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
// California, 94041, USA.
//
// UberPOV Raytracer sample file.
// Created by Christoph Lipka - 2013-08-23
// Updated by Christoph Lipka - 2016-08-24
// This file demonstrates the stochastic anti-aliasing mode.
//
// +w800 +h600 -a
// +w800 +h600 +am2 +a0.1
// +w800 +h600 +am3 +a0.1  +ac0.9 +r3
// +w800 +h600 +am3 +a0.03 +ac0.9 +r4

#version 3.7;

// also toy around with this:
#declare Focal_Blur = no;

global_settings{ assumed_gamma 1.0 }

camera {
  perspective angle 50
  location  <0.0, 3,-5.0>
  right     x*image_width/image_height
  look_at   <0.0, 1.5, 0.0>
  #if (Focal_Blur)
    focal_point <0.0, 1.0, -0.5>
    blur_samples 1, 64
    aperture 0.1
    confidence 0.9
    variance   0.1
  #end
}

light_source{ <3000,3000,-3000> color rgb 1 }

sky_sphere {
  pigment {
    gradient <0,1,0>
    color_map {
      [0.00 rgb <0.6,0.7,1.0>]
      [0.35 rgb <0.1,0.0,0.8>]
      [0.65 rgb <0.1,0.0,0.8>]
      [1.00 rgb <0.6,0.7,1.0>] 
    } 
    scale 2
  }
}

plane{ <0,1,0>, 0 
  texture{
    pigment{ checker color rgb 1 color rgb 0 }
    finish {
      ambient 0.1
      diffuse 0.9
      phong 0.1
    }
  }
}

sphere { <0,1,0>, 1
  texture {
    pigment { rgb <0.2, 0.2, 0.2> }
    finish {
      ambient 0.1
      diffuse 0.7
      brilliance 6.0
      reflection 0.6
      phong 0.8
      phong_size 120
    }
  }
}

  
