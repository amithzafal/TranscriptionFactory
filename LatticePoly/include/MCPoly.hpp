//
//  MCPoly.hpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCPoly_hpp
#define MCPoly_hpp

#include "MCTadUpdater.hpp"


class MCPoly
{
public:
	MCPoly(MCLattice*);
	virtual ~MCPoly();

	void Init(int);
	void GenerateRandom(int);
	
	void ToVTK(int);
	void FromVTK(int);

	void TrialMove(double*);
	void AcceptMove();
	void OriginMove();
	void ForkMove();
	
	int Ntad;
	int Nbond;
	
	MCLattice* lat;
	MCTadUpdater* tadUpdater;
		
protected:
	MCTad* tadTrial;
	
	std::vector<MCTad> tadConf;
	std::vector<MCBond> tadTopo;
	
	void CreateBond(MCBond&);
	void FixPBCPair(std::vector<double3>&, int, int);
	
	virtual std::vector<double3> GetPBCConf();
	double3 GetPBCCenterMass(std::vector<double3>::iterator, std::vector<double3>::iterator);
	
	double3 centerMass;
};


#endif /* MCPoly_hpp */
