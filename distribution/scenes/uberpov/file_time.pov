// This work is licensed under the Creative Commons Attribution 3.0 Unported License.
// To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
// or send a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
// California, 94041, USA.
//
// UberPOV Raytracer sample file.
// Created by Christoph Lipka - 2013-08-20
// Updated by Christoph Lipka - 2016-08-24
// This file demonstrates the file_time feature
//
// +w800 +h600 +a0.3

#version unofficial patch 3.7;
#patch "upov-file_time" 0.9;

global_settings {
  assumed_gamma 1.0
  max_trace_level 25
  radiosity {
    pretrace_start 0.1
    pretrace_end   0.01
    count 160
    nearest_count 10
    error_bound 0.8
    recursion_limit 1
    low_error_factor .5
    gray_threshold 0.0
    minimum_reuse 0.005
    maximum_reuse 0.05
    brightness 1
    adc_bailout 0.005
  }
}

default {
  finish { diffuse 0.7 specular 0.3 }
}

camera {
  right x*image_width/image_height
  location  <0,3,-4>
  look_at   <0,0,0>
  angle 20
  rotate -y*30
}

light_source {
  <100,200,100>
  color rgb 1
  area_light 300*x,300*y,17,17 adaptive 1 circular orient jitter
}

sky_sphere {
  pigment {
    gradient y
    color_map {
      [0.0 rgb <0.6,0.7,1.0>]
      [0.7 rgb <0.0,0.1,0.8>]
    }
  }
}

plane {
  y, 0
  texture {
    pigment { color srgb 0.7 }
  }
}

// ----------------------------------------

union {

  text {
    ttf "crystal.ttf", "file_time.pov was last modified on:", 0.1, 0.0
    translate <0,2>
  }

  text {
    ttf "crystal.ttf", datetime(file_time("file_time.pov"), "%Y-%m-%d %H:%M"), 0.1, 0.0
    translate <10,1>
  }

  text {
    ttf "crystal.ttf", "file_time.png was last rendered on:", 0.1, 0.0
    translate <0,-1>
  }

  text {
    ttf "crystal.ttf", /* #if (file_exists("file_time.png")) */ datetime(file_time("file_time.png"), "%Y-%m-%d %H:%M") /* #else concat("-never-") #end */ , 0.1, 0.0
    translate <10,-2>
  }
  
  texture {
    pigment { color srgb 0.7 }
  }
  
  rotate x*90
  translate y*0.2
  translate x*-10
  scale 0.1
}

// ----------------------------------------



#debug concat(datetime(file_time("file_time.pov"), "%Y-%m-%d"),"\n")
#debug concat(datetime(file_time("file_time.png"), "%Y-%m-%d"),"\n")
