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
// +w800 +h600 +am3 +a0.1  +ac0.9  +r3
// +w800 +h600 +am3 +a0.01 +ac0.99 +r6

#version 3.7;

// also toy around with this:
#declare Focal_Blur = no;

global_settings {
  assumed_gamma 1.0
  mm_per_unit 10
  subsurface { samples 100,10 }
}

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

light_source {
  <30,30,-30> color rgb 1
  area_light 5*x, 5*y, 9,9 adaptive 1 circular orient jitter
}

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
      diffuse 0.7
      specular albedo 0.1
      roughness 0.001
    }
  }
}

sphere { <0,0,0>, 1
  texture {
    checker
    texture {
      pigment { rgbt 1 }
      finish { ambient 0 diffuse 0 specular 0 }
    }
    texture {
      pigment { rgb <0.2, 0.2, 0.2> }
      finish {
        ambient 0.1
        diffuse 0.7
        brilliance 6.0
        reflection 0.6
        specular albedo 0.8
        roughness 0.001
      }
    }
    scale <0.25,0.05,0.1>
    warp { spherical }
  }
  interior_texture {
    checker
    texture {
      pigment { rgbt 1 }
      finish { ambient 0 diffuse 0 specular 0 }
    }
    texture {
      pigment { rgb 1 }
    }
    scale <0.25,0.05,0.1>
    warp { spherical }
  }
  translate <0.5,1,0>
}

sphere { <0,0,0>, 1
  texture {
    pigment { rgb <1.0, 0.7, 0.2> }
    finish {
      ambient 0.1
      diffuse 0.7
      specular albedo 0.1
      roughness 0.02
      subsurface {
        translucency 0.2
      }
    }
  }
  interior { ior 1.5 }
  translate z*2
  rotate -y*30
  translate <0.5,1,0>
}
