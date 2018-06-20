#include "PositionReco.h"

//**********Utils*************************************************************************
//----------Begin*************************************************************************
bool PositionReco::Begin(CfgManager& opts, uint64* index)
{
    //---create a position tree
    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : true;
    RegisterSharedData(new TTree("position", "position_tree"), "position_tree", storeTree);
    positionTree_ = PositionTree(index, (TTree*)data_.back().obj, nPlanes_);
    positionTree_.Init();
    
    return true;
}

bool PositionReco::ProcessEvent(const H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{

  // Find position beam
    int pattern[2][2][2];
    for(unsigned int ip=0; ip<h4Tree.nPatterns; ip++) {
      int plane=-1, axis=-1, lsb=-1;
      if(h4Tree.patternBoard[ip]==0x08020001)plane=1;
      if(h4Tree.patternBoard[ip]==0x08020002)plane=0;
      if(h4Tree.patternChannel[ip]<2) axis=1;
      if(h4Tree.patternChannel[ip]>=2) axis=0;
      lsb=h4Tree.patternChannel[ip]%2;
      if(plane<0 || axis<0 || lsb<0) {
        printf("Strange pattern info : %d %d %d\n",plane,axis,lsb);
      } else
        pattern[plane][axis][lsb]=h4Tree.pattern[ip];
    }

    double pos[nPlanes_][nAxis_][nFibers_];
    int fiber[nPlanes_][nAxis_][nFibers_];
    int nhit[nPlanes_][nAxis_]={{0,0},{0,0}}; 

    for(int iaxis=0; iaxis<nAxis_; iaxis++){
      for(int iplane=0; iplane<nPlanes_; iplane++) {
        // cout<<"plane "<<iplane<<", axis "<<iaxis<<": "<<pattern[iplane][iaxis][0]<<" "<<pattern[iplane][iaxis][1]<<endl;
        for(int ib=0; ib<nFibers_; ib++) {
          fiber[iplane][iaxis][ib]=0;
          pos[iplane][iaxis][ib]=0;
        }
        for(int ib=0; ib<(nFibers_/2); ib++) {
          int ipos=fiberorder_[0][ib]-1;
          if(((pattern[iplane][iaxis][0]>>ib)&0b1)==1)fiber[iplane][iaxis][ipos]=1;
          ipos=fiberorder_[1][ib]-1;
          if(((pattern[iplane][iaxis][1]>>ib)&0b1)==1)fiber[iplane][iaxis][ipos]=1;
        }
        // for(int ib=31; ib>=0; ib--) cout<<fiber[iplane][iaxis][(nFibers_/2)+ib];
        // cout<<" ";
        // for(int ib=31; ib>=0; ib--) cout<<fiber[iplane][iaxis][ib];
        // cout<<endl;
        
      // Search for clusters and compute positions
        int cluster=0;
        int clust_start[nFibers_], clust_size[nFibers_];
        for(int ib=0; ib<nFibers_; ib++) {
          clust_start[ib]=0;
          clust_size[ib]=0;
          // cout<<"iplane="<<iplane<<", iaxis="<<iaxis<<", ib="<<ib<<", fiber[iplane][iaxis][ib]="<<fiber[iplane][iaxis][ib]<<endl;
          if(fiber[iplane][iaxis][ib]==1) {
            if(clust_size[nhit[iplane][iaxis]]==0) clust_start[nhit[iplane][iaxis]]=ib;
            clust_size[nhit[iplane][iaxis]]++;
            cluster=1;
          } else if(fiber[iplane][iaxis][ib]==0 && cluster==1) {
            nhit[iplane][iaxis]++;
            clust_size[nhit[iplane][iaxis]]=0;
            cluster=0;
          } else
            cluster=0;
        }
        // cout<<"iplane="<<iplane<<", iaxis="<<iaxis<<", clustSize="<<clust_size[0]<<", clustStart="<<clust_start[0]<<endl;

        if (iaxis==0) {
            positionTree_.nFibresOnX[iplane] = clust_size[0]; 
            // cout<<"nFibresOnX["<<iplane<<"]="<<positionTree_.nFibresOnX[iplane]<<endl;
        } else if (iaxis==1) {
            positionTree_.nFibresOnY[iplane] = clust_size[0]; 
            // cout<<"nFibresOnY["<<iplane<<"]="<<positionTree_.nFibresOnY[iplane]<<endl;
        }

        // cout<<"iplane "<<iplane<<", iaxis "<<iaxis<<": ";
        for(int ihit=0; ihit<nhit[iplane][iaxis]; ihit++) {
          pos[iplane][iaxis][ihit]=0.5*((clust_start[ihit]-(nFibers_/2.))+clust_size[ihit]/2.)+offset_[iplane][iaxis];
          // cout<<"hit "<<ihit<<" = "<<pos[iplane][iaxis][ihit] <<" mm, ";
        }
        // cout<<endl;
      }
    }

    for(int iplane=0; iplane<nPlanes_; iplane++) {
        positionTree_.X[iplane] = 0;
        positionTree_.Y[iplane] = 0;
    	if(nhit[iplane][0]==1) positionTree_.X[iplane] = pos[iplane][0][0];
    	if(nhit[iplane][1]==1) positionTree_.Y[iplane] = pos[iplane][1][0];
    }
    // cout<<"X0="<<pos[0][0][0]<<", X1="<<pos[1][0][0]<<", Y0="<<pos[0][1][0]<<", Y1="<<pos[1][1][0]<< endl;
    // cout<<"tree: X0="<<positionTree_.X[0]<<", X1="<<positionTree_.X[1]<<", Y0="<<positionTree_.Y[0]<<", Y1="<<positionTree_.Y[1]<< endl;


  // Ask for hits in plane 0 and plane 1 at the same place
    int i0_best[2]={-1,-1}, i1_best[2]={-1,-1};
    double pos_diff[2]={999.,999.};
    double pos_best[2]={-99.,-99.}; 
    for(int iaxis=0; iaxis<nAxis_; iaxis++) {
      for(int i0=0; i0<nhit[0][iaxis]; i0++) {
        for(int i1=0; i1<nhit[1][iaxis]; i1++) {
          if(fabs(pos[0][iaxis][i0]-pos[1][iaxis][i1])<pos_diff[iaxis]) {
            pos_diff[iaxis]=fabs(pos[0][iaxis][i0]-pos[1][iaxis][i1]);
            i0_best[iaxis]=i0;
            i1_best[iaxis]=i1;
            pos_best[iaxis]=(pos[0][iaxis][i0]+pos[1][iaxis][i1])/2.;
          }
        }
      }
    }
    positionTree_.iX_best[0] = i0_best[0];
    positionTree_.iY_best[0] = i0_best[1];
    positionTree_.iX_best[1] = i1_best[0];
    positionTree_.iY_best[1] = i1_best[1];
    positionTree_.posX_best = pos_best[0];
    positionTree_.posY_best = pos_best[1];
    // cout<<"Best position found X : i0="<<i0_best[0]<<", i1="<<i1_best[0]<<", dx="<<pos_diff[0]<<", X="<<pos_best[0]<<endl;
    // cout<<"Best position found Y : i0="<<i0_best[1]<<", i1="<<i1_best[1]<<", dy="<<pos_diff[1]<<", Y="<<pos_best[1]<<endl;


    //---fill output tree
    positionTree_.Fill();
    
    return true;
}
