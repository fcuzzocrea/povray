//******************************************************************************
///
/// @file parser/patch.h
///
/// @todo   What's in here?
///
/// @copyright
/// @parblock
///
/// UberPOV Raytracer version 1.37.
/// Copyright 2013-2016 Christoph Lipka.
///
/// UberPOV is free software: you can redistribute it and/or modify
/// it under the terms of the GNU Affero General Public License as
/// published by the Free Software Foundation, either version 3 of the
/// License, or (at your option) any later version.
///
/// UberPOV is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU Affero General Public License for more details.
///
/// You should have received a copy of the GNU Affero General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.
///
/// ---------------------------------------------------------------------------
///
/// UberPOV is based on the popular Persistence of Vision Ray Tracer
/// ('POV-Ray') version 3.7, Copyright 1991-2013 Persistence of Vision
/// Raytracer Pty. Ltd.
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

#ifndef POVRAY_PARSER_PATCH_H
#define POVRAY_PARSER_PATCH_H

namespace pov
{

/// Checks whether the given patch is supported by this branch, and at which version.
unsigned int GetPatchVersion(char* patchName);

} // end of pov namespace

#endif // POVRAY_PARSER_PATCH_H
