// This work is licensed under the Creative Commons Attribution 3.0 Unported License.
// To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
// or send a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
// California, 94041, USA.
//
// UberPOV Raytracer sample file.
// Created by Christoph Lipka - 2013-08-20
// Updated by Christoph Lipka - 2016-08-24
// This file demonstrates the plaintext read feature
//
// +w800 +h600 +a0.3

#version unofficial patch 3.7;
#patch "upov-read-text" 1.0;

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

#macro SUBSTR(Str,P,L)
  #local Pos = max(1,P);
  #local Len = min(L,strlen(Str)-(Pos-1));
  #if (Len>0)
    #local Return = substr(Str,Pos,Len);
  #else
    #local Return = "";
  #end
  Return
#end

#macro NEXT_TOKEN(Str)
  #if (SUBSTR(Str,1,2) = "//")
    #local Return = Str;
  #elseif(SUBSTR(Str,1,1) = " ")
    #local L = 1;
    #while(SUBSTR(Str,1+L,1) = " ") #local L=L+1; #end
    #local Return = SUBSTR(Str,1,L)
  #elseif(SUBSTR(Str,1,1) = "\"")
    #local L = 1;
    #while((L < strlen(Str)) & (SUBSTR(Str,1+L,1) != "\"")) #if(SUBSTR(Str,1+L,1) = "\\") #local L=L+1; #end #local L=L+1; #end
    #local L=L+1;
    #local Return = SUBSTR(Str,1,L)
  #elseif(((SUBSTR(Str,1,1) >= "a") & (SUBSTR(Str,1,1) <= "z")) | ((SUBSTR(Str,1,1) >= "A") & (SUBSTR(Str,1,1) <= "Z")) | (SUBSTR(Str,1,1) = "_"))
    #local L = 1;
    #while(((SUBSTR(Str,1+L,1) >= "a") & (SUBSTR(Str,1+L,1) <= "z")) | ((SUBSTR(Str,1+L,1) >= "A") & (SUBSTR(Str,1+L,1) <= "Z")) | ((SUBSTR(Str,1+L,1) >= "0") & (SUBSTR(Str,1+L,1) <= "9")) | (SUBSTR(Str,1+L,1) = "_")) #local L=L+1; #end
    #local Return = SUBSTR(Str,1,L)
  #elseif((SUBSTR(Str,1,1) >= "0") & (SUBSTR(Str,1,1) <= "9"))
    #local L = 1;
    #while(((SUBSTR(Str,1+L,1) >= "0") & (SUBSTR(Str,1+L,1) <= "9")) | (SUBSTR(Str,1+L,1) = ".")) #local L=L+1; #end
    #local Return = SUBSTR(Str,1,L)
  #else
    #local Return = SUBSTR(Str,1,1);
  #end
  Return
#end

#macro TOKEN_COLOR(Str)
  #if((Str="adc_bailout")|(Str="adaptive")|(Str="angle")|(Str="area_light")|(Str="assumed_gamma")|(Str="brightness")
     |(Str="camera")|(Str="circular")|(Str="color")|(Str="color_map")|(Str="count")|(Str="declare")|(Str="default")|(Str="defined")|(Str="diffuse")
     |(Str="else")|(Str="elseif")|(Str="end")|(Str="error_bound")|(Str="finish")|(Str="fclose")|(Str="fopen")
     |(Str="global_settings")|(Str="gradient")|(Str="gray_threshold")|(Str="if")
     |(Str="light_source")|(Str="local")|(Str="location")|(Str="look_at")|(Str="low_error_factor")
     |(Str="macro")|(Str="max")|(Str="max_extent")|(Str="max_trace_level")|(Str="maximum_reuse")|(Str="min")|(Str="minimum_reuse")
     |(Str="nearest_count")|(Str="object")|(Str="orient")|(Str="patch")|(Str="pigment")|(Str="plane")|(Str="pretrace_start")|(Str="pretrace_end")
     |(Str="radiosity")|(Str="read")|(Str="recursion_limit")|(Str="rgb")|(Str="right")|(Str="rotate")
     |(Str="scale")|(Str="sky_sphere")|(Str="specular")|(Str="srgb")|(Str="strlen")|(Str="substr")
     |(Str="texture")|(Str="text")|(Str="translate")|(Str="ttf")|(Str="union")|(Str="version")|(Str="while")|(Str="x")|(Str="y")|(Str="z"))
    <0.6,0.2,0.6>
  #elseif(SUBSTR(Str,1,1) = "#")
    <0.6,0.2,0.6>
  #elseif((SUBSTR(Str,1,1) >= "0") & (SUBSTR(Str,1,1) <= "9"))
    <0.2,0.5,0.5>
  #elseif(SUBSTR(Str,1,1) = "\"")
    <0.8,0,0>
  #elseif((Str="+")|(Str="-")|(Str="*")|(Str="/")|(Str="=")|(Str="<")|(Str=">")|(Str="&")|(Str="|"))
    <0.8,0,0>
  #elseif((Str="{")|(Str="}"))
    <0,0,0.9>
  #elseif(SUBSTR(Str,1,2) = "//")
    <0,0.6,0>
  #else
    0.3
  #end
#end

#declare TT  = text { ttf "crystal.ttf", "TT",  0.1, 0.0 }
#declare T_T = text { ttf "crystal.ttf", "T T", 0.1, 0.0 }
#declare Width = max_extent(T_T).x-max_extent(TT).x;

union {

  #declare Y  = 0;
  #declare Y0 = 0;
  #fopen F "read_text.pov" read
  #while (defined(F))

    #read (F,text,S) //!

    #if (SUBSTR(S,strlen(S)-2,3) = "//!") #declare Y0 = Y; #end

    #declare X = 0;
    #if(Y > 120) #while (strlen(S) > 0)
      #declare Token = NEXT_TOKEN(S);
      text {
        ttf "crystal.ttf", Token, 0.1, 0.0
        #if ((Token="text") & (SUBSTR(S,strlen(S)-2,3) = "//!")) pigment { color srgb <0.6,0.2,0.6> } finish { emission <5,0,0> }
        #else pigment { color srgb TOKEN_COLOR(Token) }
        #end
        translate <X,-Y>
      }
      #declare S = substr(S,strlen(Token)+1,strlen(S)-strlen(Token));
      #declare X = X + Width * strlen(Token);
    #end #end
    #declare Y = Y + 1;

  #end
  #fclose F
  
  texture {
    pigment { color srgb 0.7 }
  }
  
  rotate x*90
  translate y*0.2
  translate x*-8
  translate z*Y0
  scale 0.1
}

// ----------------------------------------
