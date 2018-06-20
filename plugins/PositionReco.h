#ifndef __POSITION_RECO__
#define __POSITION_RECO__

#include "interface/PluginBase.h"
#include "interface/PositionTree.h"

using namespace std;

class PositionReco: public PluginBase
{
public:
    //---ctors---
    PositionReco() {};

    //---dtor---
    ~PositionReco() {};

    // //---utils---
    bool Begin(CfgManager& opts, uint64* index);
    bool ProcessEvent(const H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(CfgManager& opts) { return true; };
    
private:
    PositionTree     positionTree_;
    int              nPlanes_=2;
    int              nAxis_=2;
    int              nFibers_=64;
    int              fiberorder_[2][32]={{
                        31, 29, 23, 21,  5,  7, 15, 13,
                        1,  3, 11,  9,  6,  8, 16, 14,
                        17, 27, 19, 25, 24, 22, 32, 30,
                        4,  2, 12, 10, 20, 18, 28, 26,
                        },{
                        54, 64, 56, 62, 49, 51, 59, 57,
                        53, 55, 63, 61, 45, 47, 37, 39,
                        34, 42, 36, 44, 50, 52, 58, 60,
                        38, 48, 40, 46, 41, 43, 33, 35,
                        }};
    double           offset_[2][2]={{3.21, 0.00},{1.02,0.01}};

};

DEFINE_PLUGIN(PositionReco);

#endif
