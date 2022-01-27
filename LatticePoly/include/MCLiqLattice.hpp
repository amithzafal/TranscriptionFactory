//
//  MCLiqLattice.hpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCLiqLattice_hpp
#define MCLiqLattice_hpp

#include "MCLattice.hpp"


class MCLiqLattice: public MCLattice
{
public:
	void Init(int);

	void ToVTK(int);
	void FromVTK(int);

	void TrialMove(double*);
	void AcceptMove();
	
	double GetCouplingEnergy(const int[Ntot]) const;

	int spinTable[Ntot];
	int OriginCheck(std::vector<int>);
	
	int nLiq;
	bool stop_update = false;


private:
	void GenerateRandom();
	void GenerateDroplets();
	
	void DisplaceSpins();
	
	double GetSpinEnergy() const;

	int lookupTable[Ntot];

	int v1;
	int v2;
	std::vector<int> spinConf;
	std::vector<double3> spinDisp;
	std::vector<int> SpinLocked;
};


#endif /* MCLiqLattice_hpp */
