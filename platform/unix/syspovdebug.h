//******************************************************************************
///
/// @file platform/unix/syspovdebug.h
///
/// This header file is included by povdebug.h, which is in turn included by
/// all C++ files in POV-Ray. povdebug.h is the last header file included
/// (with the possible exception of files that do not declare anything that
/// could clash with declarations in povdebug.h or this file).
///
/// As a rule, system header files are not safe to include after this file.
///
/// @author Christopher J. Cason
///
/// @copyright
/// @parblock
///
/// Persistence of Vision Ray Tracer ('POV-Ray') version 3.8.
/// Copyright 1991-2019 Persistence of Vision Raytracer Pty. Ltd.
///
/// POV-Ray is free software: you can redistribute it and/or modify
/// it under the terms of the GNU Affero General Public License as
/// published by the Free Software Foundation, either version 3 of the
/// License, or (at your option) any later version.
///
/// POV-Ray is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU Affero General Public License for more details.
///
/// You should have received a copy of the GNU Affero General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.
///
/// ----------------------------------------------------------------------------
///
/// POV-Ray is based on the popular DKB raytracer version 2.12.
/// DKBTrace was originally written by David K. Buck.
/// DKBTrace Ver 2.0-2.12 were written by David K. Buck & Aaron A. Collins.
///
/// @endparblock
///
//******************************************************************************

#ifndef POVRAY_UNIX_SYSPOVDEBUG_H
#define POVRAY_UNIX_SYSPOVDEBUG_H

#include <cstring>
#include <fstream>
#include <cstdio>
#include <chrono>
#include <ctime>

#include <boost/format.hpp>

inline void InitStringToFile( std::string logMsg )
{
    std::string filePath = "depth.pgm";

    std::ofstream ofs(filePath.c_str(), std::ios_base::out);
    ofs << logMsg << '\n';
    ofs.close();
}

inline void DebugStringToFile( std::string logMsg )
{
    std::string filePath = "depth.pgm";

    std::ofstream ofs(filePath.c_str(), std::ios_base::out | std::ios_base::app );
    ofs << logMsg;  // No << '\n' with .pgm format as we need to control newlines.
    ofs.close();
}

#endif // POVRAY_UNIX_SYSPOVDEBUG_H
