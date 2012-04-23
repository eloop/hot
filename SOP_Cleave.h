/*
 * NAME:	SOP library (C++)
 *
 * COMMENTS:	CLEAVE SOP
 */


#ifndef __SOP_Cleave_h__
#define __SOP_Cleave_h__

#include <SOP/SOP_Node.h>

class SOP_Cleave : public SOP_Node
{
    public:

        SOP_Cleave(OP_Network *net, const char *name, OP_Operator *op);
        virtual ~SOP_Cleave();

        static PRM_Template myTemplateList[];

        static OP_Node* myConstructor(OP_Network*, const char *, OP_Operator *);

    protected:

        virtual const char* inputLabel(unsigned idx) const;

        // Method to cook geometry for the SOP
        virtual OP_ERROR cookMySop(OP_Context &context);

    private:

        void getGroups(UT_String &str) { evalString(str, 0, 0, 0); }
        int  FREQ(float t)             { return evalInt(1, 0, t); }

        float ANG(float t)             { return evalFloat(2, 0, t); }
        float DANG(float t)            { return evalFloat(3, 0, t); }
        float RANG(float t)            { return evalFloat(4, 0, t); }

        int SEED(float t)              { return evalInt(5, 0, t); }
        int SPLITS(float t)            { return evalInt(6, 0, t); }
        int RSPLITS(void)              { return evalInt(7, 0, 0); }

        float DIST(float t)            { return evalFloat(8, 0, t); }
        float DSCALE(float t)          { return evalFloat(9, 0, t); }
        float STHRESH(float t)         { return evalFloat(10, 0, t); }

};

#endif
