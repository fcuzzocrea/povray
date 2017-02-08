// This work is licensed under the Creative Commons Attribution-ShareAlike 3.0 Unported License.
// To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/3.0/ or send a
// letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

// Persistence of Vision Ray Tracer version 3.7 Include File
// File: arraycoupleddf3s.inc
// Last updated: Decemeber 12 2016  WFP
// Description: Macros & functions for working with arrays and df3 files.
//
//   ARRAYS_ReadDF3             - Complement to arrays.inc's ARRAYS_WriteDF3.
//   ConvertImagesToDF3         - 2D images to df3 file.
//   ConvertArrayOfStringsToDF3 - Array of strings to df3 file.
//   ConvertFunctionToDF3       - Given a function creates a df3 file.
//   NxNx9DBLtoArrayEntry       - Convert POV-Ray DBL to 9 deep NxN array
//                                usually ahead of ARRAYS_WriteDF3 write and
//                                access to the values as +-64 bit depth via
//                                FnctNxNx9ArrayEntryToDBL.
//   FnctNxNx9ArrayEntryToDBL   - Function to work with +-64 bit, df3 stored,
//                                integer values.
//------------------------------------------------------------------------------

#ifndef(ARRAYCOUPLEDDF3S_INC_TEMP)
#declare ARRAYCOUPLEDDF3S_INC_TEMP = version;
#version 3.71;

#ifdef(View_POV_Include_Stack)
   #debug "including arraycoupleddf3s.inc\n"
#end

#include "arrays.inc"

//--- Generic df3 read. The df3 write was added in 3.7.0 to arrays.inc
//    as ARRAYS_WriteDF3(). The incoming array must have been pre-defined
//    at a size matching the df3 file to accurately read and store all
//    voxel values.
#macro ARRAYS_ReadDF3(_Ary,_FileName)
   #local Fn_df3   = function{pattern{density_file df3 _FileName interpolate 0}}
   #local _X       = dimension_size(_Ary,1);
   #local _Y       = dimension_size(_Ary,2);
   #local _Z       = dimension_size(_Ary,3);
   #local Dx       = 1/_X;
   #local DxVoxCtr = Dx/2;
   #local Dy       = 1/_Y;
   #local DyVoxCtr = Dy/2;
   #local Dz       = 1/_Z;
   #local DzVoxCtr = Dz/2;
   #for (Z,0,_Z-1)
     #for (Y,0,_Y-1)
       #for (X,0,_X-1)
         #local Xi=(X*Dx)+DxVoxCtr;
         #local Yi=(Y*Dy)+DyVoxCtr;
         #local Zi=(Z*Dz)+DzVoxCtr;
         #local _Ary[X][Y][Z]=Fn_df3(Xi,Yi,Zi).red;
       #end
     #end
   #end
#end

//--- Utility function with with ConvertArrayOfStringsToDF3.
#macro arraycoupleddf3sAsciiToVal(_AsciiVal)
   #switch (_AsciiVal)
      #range (asc("a"), asc("z")) #local ReturnVal = (_AsciiVal-asc("a"))/52.0; #break
      #range (asc("A"), asc("Z")) #local ReturnVal = (_AsciiVal-asc("A")+26)/52.0; #break
      #range (asc("0"), asc("9")) #local ReturnVal = (_AsciiVal-asc("0"))/10.0; #break
      #case (asc ("*")) #local ReturnVal = 1.0; #break
      #else #local ReturnVal = 0.0; #break
   #end
   (ReturnVal)
#end

//--- Translate a 'y,z' 2D array of text strings into df3 file. String length
//    provides the 'x' range for the df3.
//
//    0-9    characters => 0.0 to 0.9 df3 values.
//    *      character  => 1.0.
//    a-zA-Z characters => 52 levels from a=1/52 to Z=52/52.
//    Other  characters => 0.0.
//
//    _Ary                - The array containing text strings to convert.
//    _DF3FileName        - Name of the output df3 file.
//    _Depth              - 8 or 16. Best resolution from text only 1/52 so usually 8.
//    _Pad                - 0|<val>. If >0, we will pad all sides of the df3 with that
//                          count of 0.0 values.
//    _PadOnLeft          - 0 or 1. If 1 defined array df3 is padded on the left with
//                          0s to a width equal to the incoming right side. This is
//                          done for convenience when creating df3 for lathe work, but
//                          also due the 256 character string limit in POV-Ray SDL.
//
// Example creating a df3 profile for bowl. Isosurface later scales flatter.
//
//    #include "arrays.inc"           // For ARRAYS_WriteDF3()
//    #include "arraycoupleddf3s.inc" // Array coupled DF3 macros & functions.
//
//    #declare AryStrs = array[25][1]
//    #declare AryStrs[0][0]  = "-------------------------------------2------------";
//    #declare AryStrs[1][0]  = "------------------------------------3-------------";
//    #declare AryStrs[2][0]  = "-----------------------------------2--------------";
//    #declare AryStrs[3][0]  = "----------------------------------2---------------";
//    #declare AryStrs[4][0]  = "---------------------------------2----------------";
//    #declare AryStrs[5][0]  = "--------------------------------3-----------------";
//    #declare AryStrs[6][0]  = "-------------------------------4------------------";
//    #declare AryStrs[7][0]  = "-------------------------------5------------------";
//    #declare AryStrs[8][0]  = "------------------------------4-------------------";
//    #declare AryStrs[9][0]  = "-----------------------------5--------------------";
//    #declare AryStrs[10][0] = "----------------------------6---------------------";
//    #declare AryStrs[11][0] = "----------------------------6---------------------";
//    #declare AryStrs[12][0] = "----------------------------6---------------------";
//    #declare AryStrs[13][0] = "----------------------------6---------------------";
//    #declare AryStrs[14][0] = "---------------------------7----------------------";
//    #declare AryStrs[15][0] = "--------------------------7-----------------------";
//    #declare AryStrs[16][0] = "-------------------------7------------------------";
//    #declare AryStrs[17][0] = "------------------------7-------------------------";
//    #declare AryStrs[18][0] = "-----------------------7--------------------------";
//    #declare AryStrs[19][0] = "-------------------6666---------------------------";
//    #declare AryStrs[20][0] = "------------------7-------------------------------";
//    #declare AryStrs[21][0] = "-----------------7--------------------------------";
//    #declare AryStrs[22][0] = "----------------7---------------------------------";
//    #declare AryStrs[23][0] = "---------------7----------------------------------";
//    #declare AryStrs[24][0] = "666666666666666-----------------------------------";
//    //
//    //                         <---------------------- 50 ---------------------->
//
//    ConvertArrayOfStringsToDF3(AryStrs,"StringsToDF3.df3",8,5,1)
//
#macro ConvertArrayOfStringsToDF3(_Ary,_DF3FileName,_Depth,_Pad,_PadOnLeft)
   #ifndef (_Ary)
       #error "ConvertArrayOfStringsToDF3: Array passed is not defined.\n"
   #end
   #if (dimensions(_Ary)!=2)
       #error "ConvertArrayOfStringsToDF3: Array/DF3 must have 2 dimensions.\n"
   #end
   #local _Y = dimension_size(_Ary,1);
   #local _Z = dimension_size(_Ary,2);
   #if (_Y<1)
       #error "ConvertArrayOfStringsToDF3: Array/DF3 dimension 1 must be >=1.\n"
   #end
   #if (_Z<1)
       #error "ConvertArrayOfStringsToDF3: Array/DF3 dimension 2 must be >=1.\n"
   #end
   #local _X = strlen(_Ary[0][0]);
   #if (_X<1)
       #error "ConvertArrayOfStringsToDF3: String lengths must be >=1.\n"
   #end
   #if (_Depth!=8 & _Depth!=16)
       #error "ConvertArrayOfStringsToDF3: Depth parameter must be 8 or 16.\n"
   #end
   #if (_Pad<0 | _Pad>=pow(2,16))
       #error "ConvertArrayOfStringsToDF: Pad value must be >=0 and less than max df3 range.\n"
   #end
   #if (_PadOnLeft!=0 & _PadOnLeft!=1)
       #error "ConvertImagesToDF3: PadOnLeft value must be 0 or 1.\n"
   #end

   #local _Xn  = _X+(2*_Pad);
   #local _Yn  = _Y+(2*_Pad);
   #local _Zn  = _Z+(2*_Pad);
   #if (_PadOnLeft)
     #local _Xn   = _Xn+_X;
     #local padBy = _X;
   #else
     #local padBy = 0;
   #end
   #local DF3Ary = array[_Xn][_Yn][_Zn];

   #local LookAtStringPlaneCnt=0;
   #for (Z,0,_Zn-1)
     #if (Z>=_Pad & Z<_Zn-_Pad)
       #debug concat("zPlane ",str(Z,-1,0)," StringPlane ",str(LookAtStringPlaneCnt,-1,0),"\n")
     #else
       #debug concat("zPlane ",str(Z,-1,0)," Padding with zeroes.\n")
     #end

     #for (Y,0,_Yn-1)
       #if (Y>=_Pad & Y<_Yn-_Pad)
         #local _Yr=Y-_Pad;
       #else
         #local _Yr=0;
       #end
       #if (strlen(_Ary[_Yr][min(LookAtStringPlaneCnt,_Z-1)])!=_X)
         #debug concat("\n","ConvertArrayOfStringsToDF3: String length of array entry at:","\n")
         #debug concat("[",str(_Yr,-1,0),"][",str(min(LookAtStringPlaneCnt,_Z-1),-1,0),"] must match [0][0] length of ",str(_X,-1,0),".","\n\n")
         #error "Stop."
       #end
       #for (X,0,_Xn-1)
         #if (X>=_Pad+padBy & Y>=_Pad & Z>=_Pad & X<_Xn-_Pad & Y<_Yn-_Pad & Z<_Zn-_Pad)
           #local _Xr=X-_Pad;
           // #debug concat(str(_Xr,-1,0)," ",str(_Yr,-1,0)," ",str(LookAtStringPlaneCnt,-1,0),"\n")
           #local Character=asc(substr(_Ary[_Yr][LookAtStringPlaneCnt],_Xr-padBy+1,1));
           #local DF3Ary[X][(_Yn-1)-Y][Z]=arraycoupleddf3sAsciiToVal(Character);
         #else
           #local DF3Ary[X][(_Yn-1)-Y][Z]=0.0;
         #end
       #end
     #end
     #if (Z>=_Pad & Z<_Zn-_Pad)
       #local LookAtStringPlaneCnt=LookAtStringPlaneCnt+1;
     #end
   #end
   ARRAYS_WriteDF3(DF3Ary,_DF3FileName,_Depth)
#end

//--- Translate one or more images into df3 file filtering by linear distance
//    from the passed color as a 0-1 ratio in the color unit cube.
//
//    _X                  - The x voxel grid count.
//    _Y                  - The y voxel grid count.
//    _BaseVal            - The base value to use when color is not matched.
//                          This is normally 0.0.
//    _UseVal             - To get specific value in the df3 on matching color
//                          set explicitly to >=0.0, <=1.0 value.
//                          Otherwise -1 will set voxel to the image's '.grey'.
//    _Color              - The color to match as 3 value vector.
//    _ColorNearness      - 0-1 value specifies how near image color must be to
//                          be considered a match.
//    _Weight             - 0|1 If 1, samples 12 times in circular pattern to
//                          a half voxel radius and multiplies voxel center
//                          sampled value by ratio matched by color-wise near.
//    _ImageFileNamesAry  - An array of 1 to N image file names. Image files
//                          can repeat for depth of representation in the df3.
//    _DF3FileName        - Name of the output df3 file.
//    _Depth              - 8 or 16.
//    _Pad                - 0|<val>. If >0, we will pad all sides of the df3
//                          with that count of 0.0 values.
//
// Example code to create a df3 where clouds defined white/grey on black
// background in 3 image files:
//
//    #include "arrays.inc"           // For ARRAYS_WriteDF3()
//    #include "arraycoupleddf3s.inc" // Array coupled DF3 macros & functions.
//
//    #declare ImageFileNamesAry=array[3] {
//       "CloudFront.png",
//       "CloudMiddle.png",
//       "CloutBack.png"
//    }
//    ConvertImagesToDF3(500,500,0.0,-1,<1,1,1>,1.5,0,ImageFileNamesAry,"Cloud.df3",8,5)
//
//    #error "Stopping after parse stage for DF3 creation from image."
//
#macro ConvertImagesToDF3(_X,_Y,_BaseVal,_UseVal,_Color,_ColorNearness,_Weight,
                          _ImageFileNamesAry,_DF3FileName,_Depth,_Pad)
   #if (_X<1 | _Y<1)
       #error "ConvertImagesToDF3: Passed X and Y dimensions must be >=1.\n"
   #end
   #if (_BaseVal<0.0 | _BaseVal>1.0)
       #error "ConvertImagesToDF3: Base value must be >=0.0 and <=1.0 (normally 0.0).\n"
   #end
   #if ((_UseVal!=-1) & (_UseVal<0.0 | _UseVal>1.0))
       #error "ConvertImagesToDF3: Use value must -1 or be >=0.0 and <=1.0.\n"
   #end
   #local _tmpColor = color _Color;
   #if (_tmpColor.red<0.0)
       #error "ConvertImagesToDF3: red component of color must be >=0.0.\n"
   #end
   #if (_tmpColor.green<0.0)
       #error "ConvertImagesToDF3: green component of color must be >=0.0.\n"
   #end
   #if (_tmpColor.blue<0.0)
       #error "ConvertImagesToDF3: blue component of color must be >=0.0.\n"
   #end
   #if (_ColorNearness<0.0)
       #error "ConvertImagesToDF3: ColorNearness value must be >=0.0.\n"
   #end
   #if (_Weight!=0 & _Weight!=1)
       #error "ConvertImagesToDF3: Weight value must be 0 or 1.\n"
   #end
   #ifndef (_ImageFileNamesAry)
       #error "ConvertImagesToDF3: Image array passed is not defined.\n"
   #end
   #if (dimensions(_ImageFileNamesAry)!=1)
       #error "ConvertImagesToDF3: Image array must have 1 dimension - a list.\n"
   #end
   #if (dimension_size(_ImageFileNamesAry,1)<1)
       #error "ConvertImagesToDF3: Image array dimension 1 must be >=1.\n"
   #end
   #if (_Depth!=8 & _Depth!=16)
       #error "ConvertImagesToDF3: Depth parameter must be 8 or 16.\n"
   #end
   #if (_Pad<0 | _Pad>=pow(2,16))
       #error "ConvertImagesToDF3: Pad value must be >=0 and less than max df3 range.\n"
   #end



   #local _Xn  = _X+(2*_Pad);
   #local _Yn  = _Y+(2*_Pad);
   #local _Zi  = dimension_size(_ImageFileNamesAry,1);
   #local _Zn  = _Zi+(2*_Pad);
   #local Ary  = array[_Xn][_Yn][_Zn];

   #if(file_exists(_DF3FileName) & (version>3.7))
       #local tmpDensity = density { density_file df3 _DF3FileName interpolate 0 }
       #local tmpVarDf3MapRange = max_extent(tmpDensity);                          // Requires 3.71
       #if ((_Xn!=tmpVarDf3MapRange.x) | (_Yn!=tmpVarDf3MapRange.y) |(_Zn!=tmpVarDf3MapRange.z))
           #error concat("ConvertImagesToDF3: Existing df3 ",_DF3FileName," ranges do not match calling ranges.\n"
       #end
       ARRAYS_ReadDF3(Ary,_DF3FileName)
       #local tmpOverwrite = false;
   #else
       #local tmpOverwrite = true;
   #end


   #local minXY = min(_Xn,_Yn);
   #local AryWeightVecs=array[12];
   #local AryWeightVecs[0]  = <+0.500,+0.000,0.0>*(1/minXY);
   #local AryWeightVecs[1]  = <+0.433,+0.250,0.0>*(1/minXY);
   #local AryWeightVecs[2]  = <+0.250,+0.433,0.0>*(1/minXY);
   #local AryWeightVecs[3]  = <+0.000,+0.500,0.0>*(1/minXY);
   #local AryWeightVecs[4]  = <-0.250,+0.433,0.0>*(1/minXY);
   #local AryWeightVecs[5]  = <-0.433,+0.250,0.0>*(1/minXY);
   #local AryWeightVecs[6]  = <-0.500,+0.000,0.0>*(1/minXY);
   #local AryWeightVecs[7]  = <-0.433,-0.250,0.0>*(1/minXY);
   #local AryWeightVecs[8]  = <-0.250,-0.433,0.0>*(1/minXY);
   #local AryWeightVecs[9]  = <-0.000,-0.500,0.0>*(1/minXY);
   #local AryWeightVecs[10] = <+0.250,-0.433,0.0>*(1/minXY);
   #local AryWeightVecs[11] = <+0.433,-0.250,0.0>*(1/minXY);

   #local _LookAtImageCnt=0;
   #for (Z,0,_Zn-1)
     #if (Z>=_Pad & Z<_Zn-_Pad)
       #debug concat("zPlane ",str(Z,-1,0)," ",_ImageFileNamesAry[_LookAtImageCnt],"\n")
       #ifdef (FnImagePigment)
         #undef FnImagePigment
       #end
       #local FnImagePigment=function{ pigment {image_map {concat(_ImageFileNamesAry[_LookAtImageCnt]) gamma 1.0 map_type 0 interpolate 2} } }
       #local _LookAtImageCnt=_LookAtImageCnt+1;
     #else
       #debug concat("zPlane ",str(Z,-1,0)," Padding with zeroes.\n")
     #end

     #for (Y,0,_Yn-1)
       #local _Yr=(Y-_Pad)/(_Yn-(2*_Pad));
       #for (X,0,_Xn-1)
         #if (X>=_Pad & Y>=_Pad & Z>=_Pad & X<_Xn-_Pad & Y<_Yn-_Pad & Z<_Zn-_Pad)
           #local _Xr=(X-_Pad)/(_Xn-(2*_Pad));
           #local rDiff=_Color.red  -FnImagePigment(_Xr,_Yr,0).red;
           #local gDiff=_Color.green-FnImagePigment(_Xr,_Yr,0).green;
           #local bDiff=_Color.blue -FnImagePigment(_Xr,_Yr,0).blue;
           #if (sqrt(rDiff*rDiff+gDiff*gDiff+bDiff*bDiff)<=_ColorNearness)
             #if (_Weight=1)
                #local weightMult = 1.0; // Count center as 1 as this must be true to be at this line
                #for (I,0,11)
                  #local rDiff=_Color.red  -FnImagePigment(_Xr+AryWeightVecs[I].x,_Yr+AryWeightVecs[I].y,0).red;
                  #local gDiff=_Color.green-FnImagePigment(_Xr+AryWeightVecs[I].x,_Yr+AryWeightVecs[I].y,0).green;
                  #local bDiff=_Color.blue -FnImagePigment(_Xr+AryWeightVecs[I].x,_Yr+AryWeightVecs[I].y,0).blue;
                  #if (sqrt(rDiff*rDiff+gDiff*gDiff+bDiff*bDiff)<=_ColorNearness)
                    #local weightMult = weightMult + 1;
                  #end
                #end
                #local weightMult = weightMult/13;
             #else
                #local weightMult = 1.0;
             #end
             #if (_UseVal>=0 & _UseVal<=1)
               #local Ary[X][Y][Z]=_UseVal*weightMult;
             #else
               #local Ary[X][Y][Z]=FnImagePigment(_Xr,_Yr,0).gray*weightMult;
             #end
           #else
             #if (tmpOverwrite=true)
               #local Ary[X][Y][Z]=_BaseVal;
             #end
           #end
         #else
           #if (tmpOverwrite=true)
             #local Ary[X][Y][Z]=0.0;
           #end
         #end
       #end
     #end
   #end
   ARRAYS_WriteDF3(Ary,_DF3FileName,_Depth)
#end

//--- Translate a function into a df3 file.
//
//    _X                  - The x voxel grid count to be created.
//    _Y                  - The y voxel grid count.
//    _Z                  - The z voxel grid count.
//    _BaseVal            - The base value to use when function value is not
//                          matched. This is normally 0.0.
//    _UseVal             - For specific value set explicitly as >=0.0, <=1.0.
//                          Otherwise -1 to use absolute of matched function
//                          value clamped to 0.0 - 1.0 range.
//    _FnctMinValue       - The minimum value returned by function to be
//                          considered an match.
//    _FnctMaxValue       - The maximum value for a match.
//    _Function           - The function to run over the unit cube to a
//                          resolution of _X,_Y,_Z.
//    _DF3FileName        - Name of the output df3 file.
//    _Depth              - 8 or 16.
//    _Pad                - 0|<val>. If >0, pads all sides of the df3 with
//                          that count of 0.0 values.
//
// Example code:
//
//    #include "arrays.inc"           // For ARRAYS_WriteDF3()
//    #include "arraycoupleddf3s.inc" // Array coupled DF3 macros & functions.
//
//    ConvertFunctionToDF3(50,50,50,0.0,0.8,-9e10,0.0,Fnct00,"Cloud.df3",8,5)
//
//    #error "Stopping after parse stage for DF3 creation from image."
//
#macro ConvertFunctionToDF3(_X,_Y,_Z,_BaseVal,_UseVal,_FnctMinVal,_FnctMaxVal,
                            _Function,_DF3FileName,_Depth,_Pad)
   #if (_X<1 | _Y<1 | _Z<1)
       #error "ConvertFunctionToDF3: Passed X, Y and Z dimensions must be >=1.\n"
   #end
   #if (_BaseVal<0.0 | _BaseVal>1.0)
       #error "ConvertFunctionToDF3: Base value must be >=0.0 and <=1.0 (normally 0.0).\n"
   #end
   #if ((_UseVal!=-1) & (_UseVal<0.0 | _UseVal>1.0))
       #error "ConvertFunctionToDF3: Use value must -1 or be >=0.0 and <=1.0.\n"
   #end
   #if (_FnctMinVal > _FnctMaxVal)
       #error "ConvertFunctionToDF3: Function min value must be <=0.0 max value.\n"
   #end
   #ifndef (_Function)
       #error "ConvertFunctionToDF3: Function passed is not defined.\n"
   #end
   #if (_Depth!=8 & _Depth!=16)
       #error "ConvertFunctionToDF3: Depth parameter must be 8 or 16.\n"
   #end

   #local _Xn = _X+(2*_Pad);
   #local _Yn = _Y+(2*_Pad);
   #local _Zn = _Z+(2*_Pad);
   #local Ary = array[_Xn][_Yn][_Zn];
   //---
   #for (Z,0,_Zn-1)
     #if (Z>=_Pad & Z<_Zn-_Pad)
       #debug concat("zPlane ",str(Z,-1,0)," Evaluating\n")
     #else
       #debug concat("zPlane ",str(Z,-1,0)," Padding with zeroes.\n")
     #end
     #local _Zr=(Z-_Pad)/_Z;

     #for (Y,0,_Yn-1)
       #local _Yr=(Y-_Pad)/_Y;
       #for (X,0,_Xn-1)
         #if (X>=_Pad & Y>=_Pad & Z>=_Pad & X<_Xn-_Pad & Y<_Yn-_Pad & Z<_Zn-_Pad)
           #local _Xr=(X-_Pad)/_X;
           #local tmpVal=_Function(_Xr+((1/_X)/2),_Yr+((1/_Y)/2),_Zr+((1/_Z)/2));
           #if (tmpVal>=_FnctMinVal & tmpVal<=_FnctMaxVal)
             #if (_UseVal>=0 & _UseVal<=1)
               #local Ary[X][Y][Z]=_UseVal;
             #else
               #local Ary[X][Y][Z]=min(1,max(0,abs(tmpVal)));
             #end
           #else
             #local Ary[X][Y][Z]=_BaseVal;
           #end
         #else
           #local Ary[X][Y][Z]=0.0;
         #end
       #end
     #end
   #end
   ARRAYS_WriteDF3(Ary,_DF3FileName,_Depth)
#end

//--- NxNx9 DF3 (coupled to Array) entry of converted POV-Ray DBL.
//    Expected depth is Nine 0.0 to 1.0 values being least to
//    most significant bytes of 64bits in depths 0-7 and
//    a 0 or >zero value in depth 8 with 0 meaning positive.
//    Result is DBL like in accuracy, but of course not range.
//    Macro defines a function which returns the DBL normalized
//    to +- the passed limited range. Function set up has
//    the parameters i,j,DeltaZ,range. Intended for use with one
//    byte per voxel DF3 file storage and retrieval of DBL, 64bit
//    accurate values.
//
// Example set up of function FnctDepth follows. Macro call passed
// couple and already defined array and function->density_file pattern.
//
//    #declare FnctDepth=FnctNxNx9ArrayEntryToDBL(Ary,PairedDF3Fnct);
//
//    _i     - What entry in i (x in df3, perhaps not geometrical info stored...)
//    _j     - what entry in j (y in df3)
//    _Db    - Depth of byte   (z in df3)
//    _range - Range over which stored values previously normalized.
//
// An example scene finding distance, or depth map, values is
// included in the documentation.
//
#macro FnctNxNx9ArrayEntryToDBL(_Ary,_FnctDF3)
    function (_i,_j,_Db,_range) {
      (_FnctDF3(_i,_j,(0*_Db)+(_Db/2))*255          +_FnctDF3(_i,_j,(1*_Db)+(_Db/2))*255*pow(2, 8)+
       _FnctDF3(_i,_j,(2*_Db)+(_Db/2))*255*pow(2,16)+_FnctDF3(_i,_j,(3*_Db)+(_Db/2))*255*pow(2,24)+
       _FnctDF3(_i,_j,(4*_Db)+(_Db/2))*255*pow(2,32)+_FnctDF3(_i,_j,(5*_Db)+(_Db/2))*255*pow(2,40)+
       _FnctDF3(_i,_j,(6*_Db)+(_Db/2))*255*pow(2,48)+_FnctDF3(_i,_j,(7*_Db)+(_Db/2))*255*pow(2,56)
      )/pow(2,56)*_range*select((_FnctDF3(_i,_j,(8*_Db)+(_Db/2))-0.5),1,-1)
    }
    #ifndef (_Ary)
        #error "FnctNxNx9ArrayEntryToDBL: Array passed is not defined.\n"
    #end
    #ifndef (_FnctDF3)
        #error "FnctNxNx9ArrayEntryToDBL: Function passed is not defined.\n"
    #end
    #if (dimensions(_Ary)!=3)
        #error "FnctNxNx9ArrayEntryToDBL: Array/DF3 must have 3 dimensions.\n"
    #end
    #local _II = dimension_size(_Ary,1);
    #if (_II<1)
        #error "FnctNxNx9ArrayEntryToDBL: Array/DF3 dimension 1 must be >=1.\n"
    #end
    #local _JJ = dimension_size(_Ary,2);
    #if (_JJ<1)
        #error "FnctNxNx9ArrayEntryToDBL: Array/DF3 dimension 2 must be >=1.\n"
    #end
    #local _B  = dimension_size(_Ary,3);
    #if (_B!=9)
        #error "FnctNxNx9ArrayEntryToDBL: Array/DF3 dimension 3 must be 9.\n"
    #end
#end

//--- NxNx9  POV-Ray DBL to 9 array entries (64bit+sign) conversion.
//    Expected depth is Nine 0.0 to 1.0 values being Least to
//    most significant bytes of 64bits in depths 0-7.
//    A 0 or >zero value in byte depth 8. 0 meaning the value is positive.
//    Result is > DBL like in accuracy, but of course not range.
//    Macro converts the DBL normalized to +- the passed limited range to a
//    set of 9 array entries.  The parameters are array,i,j,range,dbl.
//    Intended for use with one byte per voxel DF3 file storage
//    and retrieval of DBL, 64bit accurate values. Values larger than
//    the passed +-range are set to 0.
//
// Example macro call:
//
//    NxNx9DBLtoArrayEntry(DepthAry,I,J,Range,Value)
//
// An example scene finding distance, or depth map, values is
// included in the documentation.
//
#macro NxNx9DBLtoArrayEntry(_Ary,_i,_j,_range,_value)
    #ifndef (_Ary)
        #error "FnctNxNx9ArrayEntryToDBL: Array passed is not defined.\n"
    #end
    #if (dimensions(_Ary)!=3)
        #error "FnctNxNx9ArrayEntryToDBL: Array must have 3 dimensions.\n"
    #end
    #local _II = dimension_size(_Ary,1);
    #if (_II<1)
        #error "FnctNxNx9ArrayEntryToDBL: Array dimension 1 must be >=1.\n"
    #end
    #local _JJ = dimension_size(_Ary,2);
    #if (_JJ<1)
        #error "FnctNxNx9ArrayEntryToDBL: Array dimension 2 must be >=1.\n"
    #end
    #local _B  = dimension_size(_Ary,3);
    #if (_B!=9)
        #error "FnctNxNx9ArrayEntryToDBL: Array dimension 3 must be 9.\n"
    #end
    #if (_value<0)
        #local _Ary[_i][_j][8] = 1;
    #else
        #local _Ary[_i][_j][8] = 0;
    #end
    #if (abs(_value)>_range)
        #local _AbsFullBitsVal=0;
    #else
        #local _AbsFullBitsVal=((abs(_value)/_range)*pow(2,56));
    #end
    #declare _Ary[_i][_j][0] = bitwise_and( _AbsFullBitsVal           ,255)/255;
    #declare _Ary[_i][_j][1] = bitwise_and((_AbsFullBitsVal/pow(2, 8)),255)/255;
    #declare _Ary[_i][_j][2] = bitwise_and((_AbsFullBitsVal/pow(2,16)),255)/255;
    #declare _Ary[_i][_j][3] = bitwise_and((_AbsFullBitsVal/pow(2,24)),255)/255;
    #declare _Ary[_i][_j][4] = bitwise_and((_AbsFullBitsVal/pow(2,32)),255)/255;
    #declare _Ary[_i][_j][5] = bitwise_and((_AbsFullBitsVal/pow(2,40)),255)/255;
    #declare _Ary[_i][_j][6] = bitwise_and((_AbsFullBitsVal/pow(2,48)),255)/255;
    #declare _Ary[_i][_j][7] = bitwise_and((_AbsFullBitsVal/pow(2,56)),255)/255;
#end

//---------------------------------------------------------------------------
#version ARRAYCOUPLEDDF3S_INC_TEMP;
#end
