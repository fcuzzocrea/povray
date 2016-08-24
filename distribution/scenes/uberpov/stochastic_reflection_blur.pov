// This work is licensed under the Creative Commons Attribution 3.0 Unported License.
// To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
// or send a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
// California, 94041, USA.
//
// UberPOV Raytracer sample file.
// Created by Christoph Lipka - 2013-08-23
// Updated by Christoph Lipka - 2016-08-24
// This file demonstrates reflection blur.
//
// +w800 +h600 -a
// +w800 +h600 +am2 +a0.1
// +w800 +h600 +am3 +a0.1  +ac0.9  +r3
// +w800 +h600 +am3 +a0.01 +ac0.99 +r6

#version unofficial patch 3.7;
#patch "upov-reflection-roughness" 0.9;

// also toy around with this:
#declare Focal_Blur = no;

global_settings {
  max_trace_level 25
  assumed_gamma 1.0

  radiosity {
    pretrace_start 0.08           // start pretrace at this size
    pretrace_end   0.01           // end pretrace at this size
    count 35                      // higher -> higher quality (1..1600) [35]

    nearest_count 20              // higher -> higher quality (1..10) [5]
    error_bound 0.8               // higher -> smoother, less accurate [1.8]
    recursion_limit 1             // how much interreflections are calculated (1..5+) [3]

    low_error_factor .5           // reduce error_bound during last pretrace step
    gray_threshold 0.0            // increase for weakening colors (0..1) [0]
    minimum_reuse 0.015           // reuse of old radiosity samples [0.015]
  }
}

camera {
  perspective angle 50
  location  <0.0, 3,-5.0>
  right     x*image_width/image_height
  look_at   <0.0, 1.0, 0.0>
  #if (Focal_Blur)
    focal_point <0.0, 1.0, -0.5>
    blur_samples 1, 64
    aperture 0.1
    confidence 0.9
    variance   0.1
  #end
}

light_source {
  <3000,3000,0> color rgb 1
  area_light x*500,y*500, 9,9 adaptive 1 circular orient
}

sky_sphere {
  pigment {
    gradient <0,1,0>
    color_map {
      [0.00 srgb <0.6,0.7,1.0>]
      [0.35 srgb <0.1,0.0,0.8>]
      [0.65 srgb <0.1,0.0,0.8>]
      [1.00 srgb <0.6,0.7,1.0>] 
    } 
    scale 2
  }
}

plane{ <0,1,0>, 0 
  texture{
    pigment{ checker color rgb 1 color rgb 0.1 }
    finish {
      ambient 0.1
      diffuse 0.7
      reflection { 1.0 fresnel roughness 0.005 }
      specular albedo 0.1
      roughness 0.005
    }
  }
  interior { ior 1.5 }
}

sphere { <-1,1,2>, 1
  texture {
    pigment { rgb <0.4, 0.4, 0.5> }
    finish {
      ambient 0.0
      diffuse 0.0
      reflection { 1.0 metallic roughness 0.001 }
      specular albedo 0.1
      roughness 0.001
      metallic
    }
  }
  interior { ior 1.5 }
}

sphere { <0.5,1,0>, 1
  texture {
    pigment { rgb <1, 0.8, 0.2> }
    finish {
      ambient 0.1
      diffuse 0.7
      reflection { 1.0 fresnel roughness 0.005 }
      specular albedo 0.1
      roughness 0.005
    }
  }
  interior { ior 1.5 }
}
