/*
 *
 * NAME:	The Cleave SOP version 1.0
 *
 * AUTHOR:      Stuart Ramsden, ANUSF Vizlab, Stuart.Ramsden@anu.edu.au
 * 
 * COMMENTS:	This SOP Cleaves the geometry recursively into multiple pieces
 *
 *
 * Things to do in future:
 *               Oriented cleaving for Isotropy (compute OBB)

 *               Profile cleaving (cookie cutter) eg. jagged edge
 *               Tesselation cleaving (generalized bricker)
 *               Polygon amalgamation, either topology or color 
 *               Eventually: 3D Solid cleaving

 Copyright (C) 2003  Stuart Ramsden
 
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details,
 at http://www.gnu.org/copyleft/gpl.html
*/

// #include <UT/UT_DSOVersion.h>
#include <UT/UT_Math.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <UT/UT_Math.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_XformOrder.h>
#include <GU/GU_Detail.h>
#include <GU/GU_Align.h>
#include <GQ/GQ_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <GEO/GEO_PrimType.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>

#include <GA/GA_AttributeRef.h> 

//#include <SOP/SOP_Guide.h>
#include <tools/tools.h>
#include "SOP_Cleave.h"
#include <iostream.h>

class OP_CleaveOperator : public OP_Operator 
{
    public:
        OP_CleaveOperator();
        virtual ~OP_CleaveOperator();

        virtual void getOperatorHelp(const char *pathname, UT_String &str,int level) const;

};

OP_CleaveOperator::OP_CleaveOperator()
  : OP_Operator("cleave",
                "Cleave",
                SOP_Cleave::myConstructor,
                SOP_Cleave::myTemplateList,
                1,
                1,
                0)
{
}

OP_CleaveOperator::~OP_CleaveOperator()
{
}

void
OP_CleaveOperator::getOperatorHelp(const char *, UT_String &help, int /*level*/) const
{
    help  = "This SOP is useful for generating both regular\n";
    help  += "and randomized subdivisions of any polyconst char *, UT_String &help, int level) constgonal input\n";
    help  += "It recursively cleaves polygons\n";
    help  += "into sub-polygons.\n";
    help += "\n\n";
    help += "Parameters:\n";
    help += "   Frequency:             Number of recursive subdivisions\n";
    help += "\n";
    help += "   Initial Rotation:      Angle of first cleave\n";
    help += "   Cummulative Rotation:  Angle added to each recursive subdivision\n";
    help += "   Random Rotation:       Maximum angle perturbation added to each\n";
    help += "                          recursive subdivision\n";
    help += "\n";
    help += "   Seed:                  Seed for Randomness\n";
    help += "\n";
    help += "   Splits:                Number of polygons to cleave into at\n";
    help += "                          each subdivision step - 2 gives a Binary\n";
    help += "                          Space Partition (BSP), higher values cleave\n";
    help += "                          into N pie slices\n";
    help += "   Use Random Splits:     choose a split value randomly between\n";
    help += "                          2 and 'Splits' for each subdivision\n";
    help += "\n";
    help += "   Channel Width:         Thickness of channel to remove in subdivision\n";
    help += "   Channel Scale:         Amount to fractally decrease the channel width\n";
    help += "                          for successive subdivisions\n";
    help += "\n";
    help += "   Subdiv Threshold:      Set the area size limit of polygons, below which \n";
    help += "                          no subdivision will occur.\n";

}

void
newSopOperator(OP_OperatorTable *table)
{
     table->addOperator(new OP_CleaveOperator());
}


static PRM_Name
names[] = 
{

    PRM_Name("freq",  "Frequency"),
    PRM_Name("init",  "Initial Rotation"),
    PRM_Name("delta", "Cummulative Rotation"),
    PRM_Name("rand", "Random Rotation"),
    PRM_Name("seed", "Seed"),

    PRM_Name("splits", "Splits"),
    PRM_Name("randsplits",  "Use Random Splits"),

    PRM_Name("dist",	"Channel Width"),
    PRM_Name("dist_scale",	"Channel Scale"),

    PRM_Name("subdiv_thresh",	"Subdiv threshold"),


};

PRM_Template
SOP_Cleave::myTemplateList[] = 
{
    PRM_Template(PRM_STRING,    1, &PRMgroupName, 0, &SOP_Node::primGroupMenu),
    PRM_Template(PRM_INT_J,     1, &names[0], PRMtwoDefaults, 0, &PRMfrequency10Range),
    PRM_Template(PRM_ANGLE_J,   1, &names[1], 0, 0, &PRMangleRange),
    PRM_Template(PRM_ANGLE_J,   1, &names[2], 0, 0, &PRMangleRange),
    PRM_Template(PRM_ANGLE_J,   1, &names[3], 0, 0, &PRMangleRange),
    PRM_Template(PRM_INT_J,     1, &names[4], 0, 0, &PRMdivisionRange),
    PRM_Template(PRM_INT_J,	1, &names[5], PRMthreeDefaults, 0, &PRMorderRange),
    PRM_Template(PRM_TOGGLE,    1, &names[6]),
    PRM_Template(PRM_FLT_J,	1, &names[7], PRMzeroDefaults, 0, &PRMtoleranceRange),
    PRM_Template(PRM_FLT_J,	1, &names[8], PRMoneDefaults,  0, &PRMtoleranceRange),
    PRM_Template(PRM_FLT_J,	1, &names[9], PRMzeroDefaults, 0, &PRMtoleranceRange),
    PRM_Template(),
};


OP_Node *
SOP_Cleave::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_Cleave(net, name, op);
}

SOP_Cleave::SOP_Cleave(OP_Network *net, const char *name, OP_Operator *op)
  : SOP_Node(net, name, op) 
{
}

SOP_Cleave::~SOP_Cleave() 
{
}

OP_ERROR
SOP_Cleave::cookMySop(OP_Context &context)
{

    const GA_PrimitiveGroup  *polyGroup; 

    GEO_Primitive     	*prim;
    GQ_Detail           *gqd;
    int                  i,j,k;
    UT_Vector4           np,p;

    // Before we do anything, we must lock our inputs.  Before returning,
    //	we have to make sure that the inputs get unlocked.
    if (lockInputs(context) >= UT_ERROR_ABORT) return error();

    float now = context.getTime();
    duplicateSource(0, context, 0, 1);

    // Here we determine which groups we have to work on.  We only
    //	handle poly groups.

    UT_String groups;
    getGroups(groups);

    if (groups.isstring()) polyGroup = parsePrimitiveGroups(groups);
    else                   polyGroup = 0;
    
    if (error() >= UT_ERROR_ABORT) {
        unlockInputs();
        return error();
    }

    UT_Interrupt* boss = UTgetInterrupt();

    // Start the interrupt server
    boss->opStart("Cleaving Polys");


    // separate out all polys to be cleaved	
    GA_PrimitiveGroup* cleave_group = gdp->newPrimitiveGroup("cleave",1);
    GA_PrimitiveGroup* not_cleave_group = gdp->newPrimitiveGroup("not_cleave",1);

    if (polyGroup) {

        GA_FOR_ALL_PRIMITIVES(gdp,prim)
        {

            if ( (prim->getPrimitiveId()==GEO_PrimTypeCompat::GEOPRIMPOLY) && (polyGroup->contains(prim)!=0))
                cleave_group->add(prim);
            else
                not_cleave_group->add(prim);
        }

    } else {

        GA_FOR_ALL_PRIMITIVES(gdp,prim)
        {
            if (prim->getPrimitiveId()==GEO_PrimTypeCompat::GEOPRIMPOLY)
                cleave_group->add(prim);
            else
                not_cleave_group->add(prim);
        }
    }

    GU_Detail* scratch_gdp = new GU_Detail();
    scratch_gdp->merge(*gdp,cleave_group);

    GU_Detail* untouched_gdp = new GU_Detail();
    untouched_gdp->merge(*gdp,not_cleave_group);

    gdp->clearAndDestroy();
    gdp->merge(*scratch_gdp);
    delete scratch_gdp;

    // create rest pos attribute
    // float zeros[3] = {1,1,1};
    //GA_RWAttributeRef cleave_rest_pos = gdp->addPointAttrib("CleaveRestPos", sizeof(float)*3,
	//														GB_ATTRIB_FLOAT, (void *)zeros);
    GA_RWAttributeRef cleave_rest_pos = gdp->addFloatTuple(GA_ATTRIB_POINT, "CleaveRestPos",
														   1, GA_Defaults(1.0));

    // Pointer to the color data
    float* rpos = 0;
    GEO_Point* ppt = 0;

    UT_Vector4 pos;
    UT_Vector3 Anorm;
    UT_Vector3 Bnorm;

    // disassociate the polygons
    gdp->uniquePoints();

    // assign rest pos attribute
    for (i = 0; i < gdp->points().entries(); i++) {

        ppt = gdp->points()(i);

        // Use the index that we were returned in the example above
        pos = ppt->getPos();
        // pre 9.1 rpos = (float *)ppt->getAttribData(cleave_rest_pos);
		/*
        rpos = (float *)ppt->castAttribData<float>(cleave_rest_pos);

        rpos[0] = pos.x();
        rpos[1] = pos.y();
        rpos[2] = pos.z();
		*/
        // rpos = ppt->getValue<float*>(cleave_rest_pos);
        ppt->setValue<float>(cleave_rest_pos,pos.x(),0);
        ppt->setValue<float>(cleave_rest_pos,pos.y(),1);
        ppt->setValue<float>(cleave_rest_pos,pos.z(),2);
        // rpos[0] = pos.x();
        // rpos[1] = pos.y();
        // rpos[2] = pos.z();
    }
	

    float dist = DIST(now)/5.;
    float dscale = DSCALE(now);
    float sthresh = STHRESH(now)/10.;
    float deg = ANG(now);
    int splits = SPLITS(now);

    int curr_splits = splits;

    // create align geometry
    UT_Matrix4 mat;
    GU_AlignParms gap;
    gap.bias = 1;
    mat.identity();
    gap.postxform = mat;

    GU_Detail* align_gdp = new GU_Detail();
    gap.auxprim = align_gdp->polyGrid(2,2);

    // and align them to the XY plane
    gdp->align(gap);

    UT_XformOrder order("xyz");

    unsigned int seed = SEED(now);

    UT_Matrix4 mat4;

    float rval;

    float area = 0;
    GEO_Point* apt = 0;
    GEO_Point* bpt = 0;
    GEO_Point* cpt = 0;

    UT_Vector4 base;
    UT_Vector4 rel;

    float delta_deg;
    float arad;
    float brad;

    GU_Detail* current_gdp = 0;
    GU_Detail* orig_gdp = 0;

    for(j=0;j<FREQ(now);j++) {

        if (boss->opInterrupt()) break;

        if (RSPLITS()) {
            curr_splits = int(floor(2.5+float(splits-2)*UTfastRandom(seed)));
        }

        current_gdp = new GU_Detail;

        // loop with area calculation
        GA_FOR_ALL_PRIMITIVES(gdp,prim)
        {
            
            // Area computation
            area = 0.0;
            apt = prim->getVertex(0).getPt();

            for (i = prim->getVertexCount()-1; i >= 1; i--) {

                bpt = prim->getVertex(i).getPt();
                cpt = prim->getVertex(i-1).getPt();
 		    
                base = bpt->getPos() - apt->getPos();
                rel = cpt->getPos() - apt->getPos();
                area += (rel[0]*base[1] - rel[1]*base[0])/2.;

            }

            if(area <= sthresh) {

                current_gdp->merge(*prim);
                continue;

            }

            pos = prim->baryCenter();		
            mat4.identity();		

            rval = 2*UTfastRandomZero(seed);
		
            if (j!=0) mat4.rotate(0,0,UTdegToRad(DANG(now)+rval*RANG(now)),order);
            else      mat4.rotate(0,0,UTdegToRad(rval*RANG(now)),order);
		
             // center and rotate polys 
            for (i = prim->getVertexCount()-1; i >= 0; i--) {

                ppt = prim->getVertex(i).getPt();
                /*
                ppt->getPos() -= pos;				    
                ppt->getPos().multiply3(mat4);
*/
                ppt->setPos(ppt->getPos() - pos);
                UT_Vector4 temp_pos = ppt->getPos();
                temp_pos.multiply3(mat4);
                ppt->setPos(temp_pos);

            }

            scratch_gdp = new GU_Detail();
            scratch_gdp->merge(*prim);
            scratch_gdp->uniquePoints();
	    
            if ((curr_splits%2==0)&&(dist==0)) {

                // fast cleaving via creasing - only even cleaves, without channels
                delta_deg = 360./curr_splits;
                gqd = new GQ_Detail(scratch_gdp);

                for(k=0;k<curr_splits/2;k++) {

                    arad = UTdegToRad(deg+k*delta_deg);
                    Anorm.assign(cos(arad),sin(arad),0);
                    gqd->crease(Anorm, 0, 0, 0);

                }

                delete gqd;

                // scratch_gdp->removeUnused(); NICHOLAS
                scratch_gdp->uniquePoints();

                // accumulate result
                current_gdp->merge(*scratch_gdp);
                current_gdp->uniquePoints();

            } else {

                orig_gdp = new GU_Detail;
                orig_gdp->merge(*scratch_gdp);

                delta_deg = 360./curr_splits;
                for(k=0;k<curr_splits;k++) {

                    arad = UTdegToRad(deg+k*delta_deg);
                    Anorm.assign(cos(arad),sin(arad),0);

                    gqd = new GQ_Detail(scratch_gdp);
                    gqd->clip(Anorm, -dist, 0);
                    delete gqd;

                    // second clip - necessary for any splits > 2
                    if (splits>2) {

                        brad = UTdegToRad(deg+(k+1)*delta_deg);
                        Bnorm.assign(-cos(brad),-sin(brad),0);
                        gqd = new GQ_Detail(scratch_gdp);
                        gqd->clip(Bnorm, -dist, 0);
                        delete gqd;
                    }

                    // accumulate result
                    current_gdp->merge(*scratch_gdp);
                    scratch_gdp->clearAndDestroy();
                    if(k != (curr_splits-1)) scratch_gdp->merge(*orig_gdp);

                }

                delete orig_gdp;

            }

            delete scratch_gdp;

        }

        gdp->clearAndDestroy();
        gdp->merge(*current_gdp);

        dist *= dscale;

    }

    delete align_gdp;
    // restore rest positions
    for (i = 0; i < gdp->points().entries(); i++) {

        ppt = gdp->points()(i);

        // Use the index that we were returned in the example above
        pos = ppt->getPos();
        // pre 9.1 - rpos = (float *)ppt->getAttribData(cleave_rest_pos);
		/*
        rpos = (float *)ppt->castAttribData<float>(cleave_rest_pos);
        ppt->getPos()(0) = rpos[0];
        ppt->getPos()(1) = rpos[1];
        ppt->getPos()(2) = rpos[2];
		*/
        ppt->setValue<float>(cleave_rest_pos,rpos[0],0);
        ppt->setValue<float>(cleave_rest_pos,rpos[1],1);
        ppt->setValue<float>(cleave_rest_pos,rpos[2],2);

    }

    // delete rest pos
    gdp->destroyPointAttrib("CleaveRestPos");

    // restore non-poly primitives and untouched polys
    gdp->merge(*untouched_gdp);
    delete untouched_gdp;

    // Tell the interrupt server that we've completed
    boss->opEnd();
	
    unlockInputs();
    return error();

}


const char *
SOP_Cleave::inputLabel(unsigned) const
{
    return "Geometry to Cleave";
}
