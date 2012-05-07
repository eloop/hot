//    Hey emacs, this is -*- c++ -*- code.
//
//    SOP_Ocean.cpp - Implements a modifier SOP for building Ocean
//    waves (see Ocean.h for more details).
//
//     March 2005.
//     Drew.Whitehouse@anu.edu.au
//
//     $Id: SOP_Ocean.C 203 2007-11-20 01:49:42Z drw900 $ 
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


#include <limits.h>
#include <vector>

#include "SOP_Ocean.h"

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
#include <GU/GU_PrimPart.h>
#include <GA/GA_AttributeRef.h>

// --------------------------------- SOP paramaters --------------------------------------------

// simulation grid resolution, we want it to be power of two only 16 - 2048 allowed
static PRM_Range       resolutionRange(PRM_RANGE_RESTRICTED,4,PRM_RANGE_RESTRICTED,11);  
static PRM_Name        resolutionName("res","Ocean Resolution");
static PRM_Default     resolutionDefaults[] = {PRM_Default(8)};

// grid size
static PRM_Range       gridRange(PRM_RANGE_RESTRICTED,1.0,PRM_RANGE_FREE);  
static PRM_Name        gridName("s","Ocean Size (m)");
static PRM_Default     gridDefaults[] = {PRM_Default(200.0)};

// wind speed
static PRM_Range       windspeedRange(PRM_RANGE_RESTRICTED,0.0,PRM_RANGE_FREE);  
static PRM_Name        windspeedName("V","Windspeed (m/s)");
static PRM_Default     windspeedDefault = PRM_Default(30.0);

// wind direction
static PRM_Range       winddirRange(PRM_RANGE_RESTRICTED,0.0,PRM_RANGE_RESTRICTED,360.0);  
static PRM_Name        winddirName("w","Wind direction");
static PRM_Default     winddirDefault = PRM_Default(0.0);

// smallest wave size
static PRM_Range       smallestRange(PRM_RANGE_RESTRICTED,0.02,PRM_RANGE_FREE);  
static PRM_Name        smallestName("l","Shortest Wavelength(m)");
static PRM_Default     smallestDefault = PRM_Default(1.0);

// fiddle 
static PRM_Range       fiddleRange(PRM_RANGE_RESTRICTED,0.0001,PRM_RANGE_UI,100.0);  
static PRM_Name        fiddleName("A","Approx. Waveheight(m)");
static PRM_Default     fiddleDefault = PRM_Default(1.0);

// seed for the random number generator
static PRM_Name        seedName("seed","Seed");

// chop toggle
static PRM_Name        chopToggleName("chop_toggle","Chop");
static PRM_Default     chopToggleDefault = PRM_Default(0); 

// chop amount
static PRM_Range       chopRange(PRM_RANGE_RESTRICTED,0.0,PRM_RANGE_FREE,100.0);  
static PRM_Name        chopName("chop","Choppyness");
static PRM_Default     chopDefault = PRM_Default(1.0);

// damp amount
static PRM_Range       dampRange(PRM_RANGE_RESTRICTED,0.0,PRM_RANGE_RESTRICTED,1.0);  
static PRM_Name        dampName("damp","Damp Reflections");
static PRM_Default     dampDefault = PRM_Default(1.0);

// wind_align amount
static PRM_Range       wind_alignRange(PRM_RANGE_RESTRICTED,1.0,PRM_RANGE_RESTRICTED,10.0);  
static PRM_Name        wind_alignName("wind_align","Wind Alignment");
static PRM_Default     wind_alignDefault = PRM_Default(2.0);

// depth
static PRM_Range       depthRange(PRM_RANGE_RESTRICTED,0.01,PRM_RANGE_RESTRICTED,3000.0);  
static PRM_Name        depthName("ocean_depth","Ocean Depth");
static PRM_Default     depthDefault = PRM_Default(200.0);

// time
static PRM_Range       timeRange(PRM_RANGE_FREE,0.0,PRM_RANGE_FREE,3000.0);  
static PRM_Name        timeName("time","Time");
static PRM_Default     timeDefault = PRM_Default(0,"$T");

// interp toggle
static PRM_Name        interpToggleName("interp_toggle","Catmull-Rom Interpolation");
static PRM_Default     interpToggleDefault = PRM_Default(1); 

// normals toggle
static PRM_Name       normalsToggleName("normals_toggle","Normals");
static PRM_Default     normalsToggleDefault = PRM_Default(0); 

// jacobian toggle
static PRM_Name       jacobianToggleName("jacobian_toggle","Jacobian");
static PRM_Default     jacobianToggleDefault = PRM_Default(0); 

PRM_Template
SOP_Ocean::myTemplateList[] = 
{
  PRM_Template(PRM_INT,1,&resolutionName, resolutionDefaults,0,&resolutionRange,oceanChanged),

  PRM_Template(PRM_XYZ,1,&gridName, gridDefaults,0,&gridRange,oceanChanged),

  PRM_Template(PRM_FLT,1,&windspeedName, &windspeedDefault,0,&windspeedRange,oceanChanged),

  PRM_Template(PRM_ANGLE,1,&winddirName,&winddirDefault,0,&winddirRange,oceanChanged),

  PRM_Template(PRM_FLT,1,&smallestName, &smallestDefault,0,&smallestRange,oceanChanged),

  PRM_Template(PRM_FLT,1,&fiddleName, &fiddleDefault,0,&fiddleRange,oceanChanged),

  PRM_Template(PRM_INT,1,&seedName, 0,0,0,oceanChanged),

  PRM_Template(PRM_TOGGLE,1,&chopToggleName, &chopToggleDefault,0,0,oceanChanged), 

  PRM_Template(PRM_FLT,1,&chopName, &chopDefault,0,&chopRange), 

  PRM_Template(PRM_FLT,1,&dampName, &dampDefault,0,&dampRange,oceanChanged), 

  PRM_Template(PRM_FLT,1,&wind_alignName, &wind_alignDefault,0,&wind_alignRange,oceanChanged), 

  PRM_Template(PRM_FLT,1,&depthName, &depthDefault,0,&depthRange,oceanChanged), 

  PRM_Template(PRM_FLT_J,1,&timeName, &timeDefault,0,&timeRange,0), 

  PRM_Template(PRM_TOGGLE,1,&interpToggleName, &interpToggleDefault,0,0), 

  PRM_Template(PRM_TOGGLE,1,&normalsToggleName, &normalsToggleDefault,0,0,oceanChanged), 

  PRM_Template(PRM_TOGGLE,1,&jacobianToggleName, &jacobianToggleDefault,0,0,oceanChanged), 

  PRM_Template()
};

OP_Node *
SOP_Ocean::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
  return new SOP_Ocean(net, name, op);
}

SOP_Ocean::SOP_Ocean(OP_Network *net, const char *name, OP_Operator *op) : SOP_Node(net, name, op)
{
  _ocean = 0;
  _ocean_context = 0;
  _ocean_scale = 1.0f;

  _ocean_needs_rebuild = true;

}

SOP_Ocean::~SOP_Ocean()
{
  if (_ocean)
  {
    delete _ocean;
  }
  if (_ocean_context)
  {
    delete _ocean_context;
  }
}

unsigned
SOP_Ocean::disableParms()
{
  unsigned changes = SOP_Node::disableParms();

  changes += enableParm("normals_toggle",CHOP(0) ? 0 : 1);

  return changes;
}

int
SOP_Ocean::oceanChanged(void * data, int,  float, const PRM_Template *)
{
  SOP_Ocean* o = reinterpret_cast<SOP_Ocean*>(data);
  o->oceanNeedsRebuild();
  return 1;
}

OP_ERROR
SOP_Ocean::cookMySop(OP_Context &context)
{
  float now = context.getTime();

  //std::cout << "cook ocean, t = " << now << std::endl;

  // lock inputs
  if (lockInputs(context) >= UT_ERROR_ABORT )
  {
    return error();
  }


  GEO_Point        *ppt;
  UT_Interrupt    *boss;

  // Check to see that there hasn't been a critical error in cooking the SOP.
  if (error() < UT_ERROR_ABORT)
  {
    boss = UTgetInterrupt();

    // Start the interrupt server
    boss->opStart("Updating Ocean");

    duplicatePointSource(0,context);

    int   gridres  = 1 << int(GRID_RES(now));
    float stepsize = GRID_SIZE(now) / (float)gridres;

    bool do_chop     = CHOP(now);
    bool do_jacobian = JACOBIAN(now);
    bool do_normals  = NORMALS(now) && !do_chop;

    if (!_ocean || _ocean_needs_rebuild)
    {
      if (_ocean)
      {
        delete _ocean;
      }

      if (_ocean_context)
      {
        delete _ocean_context;
      }

      _ocean = new drw::Ocean(gridres,gridres,stepsize,stepsize,
                              V(0),L(0),1.0,W(0),1-DAMP(0),FOOALIGN(0),
                              DEPTH(0),SEED(0));
      _ocean_scale   = _ocean->get_height_normalize_factor();

      _ocean_context = _ocean->new_context(true,do_chop,do_normals,do_jacobian);

      _ocean_needs_rebuild = false;
      //             std::cout << "######### SOP, rebuilt ocean, norm_factor = " << _ocean_scale 
      //                       << " chop = " << do_chop 
      //                       << " norm = " << do_normals
      //                       << " jacobian = " << do_jacobian
      //                       << std::endl;
    }

    float chop_amount = CHOPAMOUNT(now);

    // sum up the waves at this timestep
    _ocean->update(TIME(now),*_ocean_context,true,do_chop,do_normals,do_jacobian,
                   _ocean_scale * SCALE(now),chop_amount);

    bool linterp = ! INTERP(now);


    // get our attribute indices
    GA_RWAttributeRef normal_index;
    GA_RWAttributeRef jminus_index;
    GA_RWAttributeRef eminus_index;

    if (do_normals)
    {
      normal_index = gdp->addNormalAttribute(GEO_POINT_DICT);
    }
    if (do_jacobian)
    {
      // jminus_index = gdp->addPointAttrib("mineigval",sizeof(float),GB_ATTRIB_FLOAT,0);
      // eminus_index = gdp->addPointAttrib("mineigvec",sizeof(UT_Vector3),GB_ATTRIB_VECTOR,0);
      jminus_index = gdp->addTuple(GA_STORE_REAL32,GA_ATTRIB_POINT,"mineigval",1,GA_Defaults(0));
      eminus_index = gdp->addFloatTuple(GA_ATTRIB_POINT,"mineigvec",1,GA_Defaults(0));
    }

    // this is not that fast, can it be done quicker ???
    GA_FOR_ALL_GPOINTS(gdp, ppt)
    {
      UT_Vector4 p = ppt->getPos();

                
      if (linterp)
      {
        _ocean_context->eval_xz(p(0),p(2));
      }
      else
      {
        _ocean_context->eval2_xz(p(0),p(2));
      }

      if (do_chop) 
      {
        p.assign( p(0) + _ocean_context->disp[0],
                  p(1) + _ocean_context->disp[1],
                  p(2) + _ocean_context->disp[2] );
      }
      else
      {
        p.assign(p(0), p(1)+_ocean_context->disp[1], p(2));
      }

      if (do_normals)
      {
        /*
          
          UT_Vector3* normal = (UT_Vector3*) ppt->castAttribData<UT_Vector3>(normal_index);
          normal->assign(_ocean_context->normal[0],
          _ocean_context->normal[1],
          _ocean_context->normal[2]);
          normal->normalize();
        */
        ppt->getValue<UT_Vector3>(normal_index).assign(_ocean_context->normal[0],
                                                       _ocean_context->normal[1],
                                                       _ocean_context->normal[2]);
        ppt->getValue<UT_Vector3>(normal_index).normalize();
      }

      if (do_jacobian)
      {/*
         float *js = (float*)ppt->castAttribData<float>(jminus_index);
         *js = _ocean_context->Jminus;
         UT_Vector3* eminus = (UT_Vector3*)ppt->castAttribData<UT_Vector3>(eminus_index);
         eminus->assign(_ocean_context->Eminus[0],0,_ocean_context->Eminus[1]);
       */
        ppt->setValue<float>(jminus_index,_ocean_context->Jminus);
        ppt->getValue<UT_Vector3>(eminus_index).assign(_ocean_context->Eminus[0],0,_ocean_context->Eminus[1]);
      }
      ppt->setPos(p);
    }


    gdp->notifyCache(GU_CACHE_ALL);

    // Tell the interrupt server that we've completed
    boss->opEnd();
  }

  unlockInputs();

  OP_Node::flags().setTimeDep(true);

  return error();
}

//
// Here, we sub-class off of OP_Operator so that we can embed our help in the C
// file.
//
class OP_OceanOperator : public OP_Operator 
{

 public:

  OP_OceanOperator();
  virtual ~OP_OceanOperator();

  /*     virtual void    getOperatorHelp(const char *pathname, UT_String &str, */
  /*                                     int level, bool *is_html_ptr, */
  /*                                     bool try_other_levels) const; */

};

OP_OceanOperator::OP_OceanOperator(): OP_Operator("ocean",
                                                  "Ocean",
                                                  SOP_Ocean::myConstructor,
                                                  SOP_Ocean::myTemplateList,
                                                  1,
                                                  1,
                                                  0)        
{
  // nothing
}

OP_OceanOperator::~OP_OceanOperator() 
{
  // nothing
}

// todo: make this more sylish, possibly copy the houdini help style ?
void SOP_Ocean::getHelpText (UT_String &help,int, bool *)
{

  help  = "<h3>Ocean SOP\n\n</h3>";  

  help  += "<h4>Description:</h4>";

  help  += "" 
           "<p>"
           "This SOP is an implementation of the Ocean model described in the Siggraph course "
           "<a href=\"http://www.finelightvisualtechnology.com/docs/coursenotes2004.pdf\">Simulating Ocean Surfaces</a>,"
           " by <a href=\"http://www.finelightvisualtechnology.com/index.php\">Jerry Tessendorf</a>."
           "It generates a wave surface by randomly "
           "generating a set of waves with a spectrum approximated by the heuristic equation described by "
           "Phillips. The waves are then rapidly summed together using a fast fourier transform(FFT)."
           "</p>";
  help  += "" 
           "<p>"
           "The SOP is a modifier that additively adjusts the incoming points, with wave height acting in the y direction and "
           "chop displacement in the xz plane."
           " The modifier simply wraps the computation grid onto the input grid. The periodic nature of FFT's guarantees"
           " the result tiles correctly."
           "</p>";

  help  += "<h4>Parameters:</h4>";
  
  help  += "<p>";

  help  += "Ocean Resolution - defines the resolution of the compututional grid that we calculate the waves on which "
           "is limited to be power of two,so e.g. 10 => 2^10 => 1024. Think of it as like a the number of pixels in an image texture.<br>";

  help  += "Ocean Size (m) - the physical size of the computation grid. Think of this as defining the texture coordinates on the plane, the ocean will tile with blocks of this size.<br>";
  help  += "Wind Speed (m/s)- the speed of the wind generating the wave field.<br>";
  help  += "Wind Direction - the angle to the x axis on the xz plane.<br>";
  help  += "Smallest Wavelength (m)- ignore waves smaller than this length, helps to anti-alias the waves being superimposed on a course mesh.<br>";
  help  += "Wave Height (m) - The waves are generated to range in the [-1,1] range, this scales it to the desired height.<br>";
  help  += "Seed - seeds the random number generator.<br>";
  help  += "Chop - toggles the generation of chop displacement that gives the waves a peaked appearance.<br>";
  help  += "Choppyness - scales the chop displacement, you will need to play with thil to ensure the waves are't self intersecting.<br>";
  help  += "Time - time that the surface is evaluated for. You probably want to put the expression $T in here<br>";

  help  += "<p>See <a href=\"http://www.odforce.net/wiki/index.php/HoudiniOceanToolkit\">the homepage</a> for the latest details and releases.";
}
void
newSopOperator(OP_OperatorTable *table)
{
#if 0
  // If we want to use the "standard" operator class, we can do this here.
  // However, the OP_OceanOperator has the help built in (and doesn't require
  // an external file).
  table->addOperator(
      new OP_Operator("Ocean",            // Internal name
                      "Ocean",            // UI name
                      SOP_Ocean::myConstructor,    // How to build the SOP
                      SOP_Ocean::myTemplateList,    // My parameters
                      0,                // Min # of sources
                      0,                // Max # of sources
                      SOP_Ocean::myVariables,    // Local variables
                      OP_FLAG_GENERATOR)        // Flag it as generator
                     );
#else
  table->addOperator(new OP_OceanOperator());
#endif
}


