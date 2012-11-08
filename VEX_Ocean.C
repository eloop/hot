//    Yo emacs, this is -*- c++ -*- code.
//
//    VEX_Ocean.cpp - Vex functions for building Ocean waves (see
//    Ocean.h for more details).
//
//     March 2005.
//     Drew.Whitehouse@anu.edu.au
//
//     $Id: SOP_Ocean.C 132 2005-08-23 04:56:06Z drw900 $
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
//
//		 Modified by Hans Joergen Kjaernet <hansj@stormstudios.no>, 2012.03.12, to work in H12.
//		 The modification is rather small and based upon pieces of code from Christian Schnellhammers multithreaded hot port to H12.
//		 This modification is not multi-threaded.

#include <limits.h>
#include <string>

#include "Ocean.h"

#include <UT/UT_DSOVersion.h>
#include <CMD/CMD_Manager.h>
#include <CMD/CMD_Args.h>
#include <VEX/VEX_VexOp.h>
#include <OP/OP_Director.h>
#include <UT/UT_Math.h>
#include <UT/UT_Interrupt.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <CH/CH_LocalVariable.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>

//
// vex functions implementing Tessendorf's ocean model
//
struct OceanHolder
{
    drw::Ocean        *ocean;
    drw::OceanContext *context;
    float              normalize_factor;
		float							now;

    OceanHolder() : ocean(0),context(0)
    {
        // nothing
    }

    ~OceanHolder()
    {
        if (ocean) delete ocean;
        if (context) delete context;
    }
};

static void*
ocean_init()
{
    OceanHolder *data = new OceanHolder;
    //std::cout << "Ocean VEX Init() " << std::endl << std::flush;
    return reinterpret_cast<void*>(data);
}

static void
ocean_cleanup(void* data)
{
  //std::cout << "Ocean VEX Cleanup()\n" << std::flush; // never called ?
  if (data)
    {
      delete reinterpret_cast<OceanHolder*>(data);
    }
}

static drw::Ocean*
ocean_from_argv(void *argv[])
{

    // consumes "IFFFFFFFI" ...
    int   res  = 1 << *(int *)argv[0]; // I
    float size = *(float *)argv[1];    // F
    float V    = *(float *)argv[2];    // F
    float l    = *(float *)argv[3];    // F
    float w    = *(float *)argv[4];    // F
    float damp = *(float *)argv[5];    // F
    float align= *(float *)argv[6];    // F
    float depth= *(float *)argv[7];    // F
    int   seed = *(int *)argv[8];      // I

    return new drw::Ocean(res,res,size/float(res),size/float(res),
                          V,l,0.000001,UTdegToRad(w),
                          1-damp,align,depth,seed);


}
static void
ocean_eval_ij(int, void *argv[], void *data)
{
    OceanHolder *oh = reinterpret_cast<OceanHolder*>(data);

    int    i              = *(int *)argv[0];  // F
    int    j              = *(int *)argv[1];  // F
    float  now            = *(int *)argv[2];  // F
    float  height_scale   = *(int *)argv[3];  // F

    int    do_chop        = *(int *)argv[4];    // I
    float  chop_amount    = *(float *)argv[5];  // F
    float *displacement   =  (float *)argv[6];  // &V

    int    do_normal      = *(int*)argv[7];     // I
    float *normal         =  (float*)argv[8];   // &V

    int    do_jacobian    = *(int*)argv[9];    // I
    float *Jminus         =  (float*)argv[10]; // &F
    float *Jplus          =  (float*)argv[11]; // &F
    float *Eminus         =  (float*)argv[12]; // &V
    float *Eplus          =  (float*)argv[13]; // &V

    if (!oh->ocean)
    {
        oh->ocean = ocean_from_argv(argv+14);
        oh->normalize_factor = oh->ocean->get_height_normalize_factor();
        oh->context = oh->ocean->new_context(true,do_chop,do_normal,do_jacobian);
        oh->ocean->update (now,*oh->context, true,do_chop,do_normal,do_jacobian,
                           height_scale * oh->normalize_factor,chop_amount);
				oh->now = now;

    }

		if ( oh->now != now ) {
        oh->ocean->update (now,*oh->context, true,do_chop,do_normal,do_jacobian,
                           height_scale * oh->normalize_factor,chop_amount);
				oh->now = now;
		}

    oh->context->eval_ij(i,j);

    displacement[0] = oh->context->disp[0];
    displacement[1] = oh->context->disp[1];
    displacement[2] = oh->context->disp[2];

    if (do_normal)
    {
        normal[0] = oh->context->normal[0];
        normal[1] = oh->context->normal[1];
        normal[2] = oh->context->normal[2];
    }
    if (do_jacobian)
    {
        *Jminus  = oh->context->Jminus;
        *Jplus   = oh->context->Jplus;

        Eminus[0] = oh->context->Eminus[0];
        Eminus[1] = oh->context->Eminus[1];
        Eminus[2] = oh->context->Eminus[2];

        Eplus[0] = oh->context->Eplus[0];
        Eplus[1] = oh->context->Eplus[1];
        Eplus[2] = oh->context->Eplus[2];

    }

}


static void
ocean_eval(int, void *argv[], void *data)
{
    OceanHolder *oh = reinterpret_cast<OceanHolder*>(data);

    float  x              = *(float *)argv[0];  // F
    float  z              = *(float *)argv[1];  // F
    float  now            = *(float *)argv[2];  // F
    float  height_scale   = *(float *)argv[3];  // F

    int    do_chop        = *(int *)argv[4];    // I
    float  chop_amount    = *(float *)argv[5];  // F
    float *displacement   =  (float *)argv[6];  // &V

    int    do_normal      = *(int*)argv[7];     // I
    float *normal         =  (float*)argv[8];   // &V

    int    do_jacobian    = *(int*)argv[9];    // I
    float *Jminus         =  (float*)argv[10]; // &F
    float *Jplus          =  (float*)argv[11]; // &F
    float *Eminus         =  (float*)argv[12]; // &V
    float *Eplus          =  (float*)argv[13]; // &V

    if (!oh->ocean)
    {
        oh->ocean = ocean_from_argv(argv+14);
        oh->normalize_factor = oh->ocean->get_height_normalize_factor();
        oh->context = oh->ocean->new_context(true,do_chop,do_normal,do_jacobian);
        oh->ocean->update (now,*oh->context, true,do_chop,do_normal,do_jacobian,
                           height_scale * oh->normalize_factor,chop_amount);

				oh->now = now;

    }

		if ( oh->now != now ) {
        oh->ocean->update (now,*oh->context, true,do_chop,do_normal,do_jacobian,
                           height_scale * oh->normalize_factor,chop_amount);
				oh->now = now;
		}

    // We always choose the catmull rom version here, the linear
    // option is primarily for the SOP where realtime feedback is more
    // important.
    oh->context->eval2_xz(x,z);

    displacement[0] = oh->context->disp[0];
    displacement[1] = oh->context->disp[1];
    displacement[2] = oh->context->disp[2];

    if (do_normal)
    {
        normal[0] = oh->context->normal[0];
        normal[1] = oh->context->normal[1];
        normal[2] = oh->context->normal[2];
    }
    if (do_jacobian)
    {
        *Jminus  = oh->context->Jminus;
        *Jplus   = oh->context->Jplus;

        Eminus[0] = oh->context->Eminus[0];
        Eminus[1] = oh->context->Eminus[1];
        Eminus[2] = oh->context->Eminus[2];

        Eplus[0] = oh->context->Eplus[0];
        Eplus[1] = oh->context->Eplus[1];
        Eplus[2] = oh->context->Eplus[2];

    }
}


void
newVEXOp(void *)
{
    new VEX_VexOp("ocean_eval_ij@IIFFIF&VI&VI&F&F&V&VIFFFFFFFI",
                  ocean_eval_ij,
                  VEX_ALL_CONTEXT,
                  ocean_init,ocean_cleanup);
    new VEX_VexOp("ocean_eval@FFFFIF&VI&VI&F&F&V&VIFFFFFFFI",
                  ocean_eval,
                  VEX_ALL_CONTEXT,
                  ocean_init,ocean_cleanup);
}
