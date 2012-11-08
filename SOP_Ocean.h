//    Hey emacs, this is -*- c++ -*- code.
//
//    SOP_Ocean.h - Implements a modifier SOP and some Vex functions
//                    for building Ocean waves (see Ocean.h for more details).
//
//     March 2005.
//     Drew.Whitehouse@anu.edu.au
//
//     $Id: SOP_Ocean.C 88 2005-06-17 03:04:46Z drw900 $
//

//     Houdini Ocean Toolkit
//     Copyright (C) 2005  Drew Whitehouse, ANU Supercomputer Facility

//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.

//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.

//     You should have received a copy of the GNU General Public License
//     along with this program; if not, write to the Free Software
//     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef __drw_sop_h
#define __drw_sop_h

#include "Ocean.h"
#include <UT/UT_Math.h>
#include <SOP/SOP_Node.h>

class SOP_Ocean : public SOP_Node
{
 public:

  static OP_Node *myConstructor(OP_Network*,const char *,OP_Operator *);

  static PRM_Template myTemplateList[];

  static	int oceanChanged(void *,int,float,const PRM_Template *);

  // oceanChanged can call this to flag a rebuild of _ocean
  void    oceanNeedsRebuild() { _ocean_needs_rebuild = true; }

 protected:

  void getHelpText (UT_String &help, int level, bool *is_html_ptr);

  // These guys are called for us, hence the protected status
  SOP_Ocean(OP_Network *net, const char *name, OP_Operator *op);
  virtual ~SOP_Ocean();

  virtual unsigned disableParms();

  virtual OP_ERROR cookMySop(OP_Context &context);

  // This is where all the wave action takes place
  drw::Ocean        *_ocean;
  drw::OceanContext *_ocean_context;
  float              _ocean_scale;

  // If this is true cook will create a new instance of drw::Ocean
  // next time it runs.
  bool _ocean_needs_rebuild;

 private:

  int     GRID_RES(float t)       { return evalInt   (0,0,t); }

  float   GRID_SIZE(float t)    { return evalFloat (1,0,t); }

  float   V(float t)        { return evalFloat (2,0,t); }

  float   W(float t)        { return UTdegToRad(evalFloat (3,0,t)); }

  float   L(float t)        { return evalFloat (4,0,t); }

  float   SCALE(float t)    { return evalFloat (5,0,t); }

  int     SEED(float t)     { return evalInt (6,0,t); }

  int     CHOP(float t)     { return evalInt(7,0,t); }

  float   CHOPAMOUNT(float t){ return evalFloat (8,0,t); }

  float   DAMP(float t)     { return evalFloat (9,0,t); }

  float   FOOALIGN(float t)    { return evalFloat (10,0,t); }

  float   DEPTH(float t)    { return evalFloat (11,0,t); }

  float   TIME(float t)    { return evalFloat (12,0,t); }

  int     INTERP(float t)     { return evalInt(13,0,t); }

  int     NORMALS(float t)      { return evalInt(14,0,t); }

  int     JACOBIAN(float t)     { return evalInt(15,0,t); }

};

#endif // __drw_sop_h
