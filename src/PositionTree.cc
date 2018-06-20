#include "interface/PositionTree.h"

PositionTree::PositionTree(uint64* idx, TTree* tree, int nPlanes)
{
    tree_ = tree ? tree : new TTree();

    n_planes = nPlanes;
    n_hitsX=0;
    n_hitsY=0;
    posX_best=-99.;
    posY_best=-99.;

    index=idx;
}

void PositionTree::Init()
{
    //---global branches
    tree_->Branch("index", index, "index/l");
    X = new float[n_planes];
    Y = new float[n_planes];
    nFibresOnX = new int[n_planes];
    nFibresOnY = new int[n_planes];
    iX_best = new float[n_planes];
    iY_best = new float[n_planes];

    //---position branches
    tree_->Branch("n_planes", &n_planes, "n_planes/I");
    tree_->Branch("X", X, "X[n_planes]/F");
    tree_->Branch("Y", Y, "Y[n_planes]/F");
    tree_->Branch("nFibresOnX", nFibresOnX, "nFibresOnX[n_planes]/I");
    tree_->Branch("nFibresOnY", nFibresOnY, "nFibresOnY[n_planes]/I");
    tree_->Branch("iX_best", iX_best, "iX_best[n_planes]/F");
    tree_->Branch("iY_best", iY_best, "iY_best[n_planes]/F");
    tree_->Branch("posX_best", &posX_best, "posX_best/F");
    tree_->Branch("posY_best", &posY_best, "posY_best/F");

}
