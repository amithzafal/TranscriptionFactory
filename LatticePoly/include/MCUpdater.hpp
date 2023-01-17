//
//  MCUpdater.hpp
//  LatticePoly
//
//  Created by mtortora on 28/02/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#ifndef MCUpdater_hpp
#define MCUpdater_hpp

#include "MCLivingPoly.hpp"
#include "MCReplicPoly.hpp"


template<class lattice>
inline bool MetropolisMove(lattice* lat, double dE)
{
	if ( dE > 0. )
	{
		double rnd = lat->rngDistrib(lat->rngEngine);
		return (rnd < exp(-dE));
	}
	
	return true;
}

template<class lattice>
inline bool RateMove(lattice* lat, double dE)
{
	if ( dE > 0. )
	{
		double rnd = lat->rngDistrib(lat->rngEngine);
		return ( rnd < dE );
		
	}
	
	return false;
}


/* struct containers are required to circumvent function template partial specialization */
// UpdateTAD template specialisations
template<class lattice, class polymer>
struct UpdateTADImpl
{
	static inline void _(lattice* lat, polymer* pol, unsigned long long* acceptCount, unsigned long long* acceptCountTopo)
	{
		double dE;
		acceptCountTopo = 0;

		pol->TrialMove(&dE);

		if ( pol->tadUpdater->legal )
		{
			double dEcpl = pol->GetCouplingEnergy(lat->spinTable);
			bool acceptMove = MetropolisMove(lat, dE+dEcpl);
			
			if ( acceptMove )
			{
				pol->AcceptMove();
				++(*acceptCount);
			}
		}
	}
};

template<class polymer>
struct UpdateTADImpl<MCLattice, polymer>
{
	static inline void _(MCLattice* lat, polymer* pol, unsigned long long* acceptCount, unsigned long long* acceptCountTopo)
	{
		double dE;
		acceptCountTopo = 0;

		pol->TrialMove(&dE);

		if ( pol->tadUpdater->legal )
		{
			double dEeff = pol->GetEffectiveEnergy();
			bool acceptMove = MetropolisMove(lat, dE+dEeff);
		
			if ( acceptMove )
			{
				pol->AcceptMove();
				++(*acceptCount);
			}
		}
	}
};

template<>
struct UpdateTADImpl<MCLattice, MCPoly>
{
	static inline void _(MCLattice* lat, MCPoly* pol, unsigned long long* acceptCount, unsigned long long* acceptCountTopo)
	{
		double dE;
		
		pol->TrialMove(&dE);

		if ( pol->tadUpdater->legal )
		{
			bool acceptMove = MetropolisMove(lat, dE);
		
			if ( acceptMove )
			{
				pol->AcceptMove();
				++(*acceptCount);
			}
		}
		

		pol->TrialMoveTopo();
				
		if ( pol->tadUpdater->legalTopo2 )
		{	
			bool acceptTopoMove = RateMove(lat, TopoRate);
			if ( acceptTopoMove )
			{
				pol->AcceptMoveTopo();
				++(*acceptCountTopo);
			}	
		}

	}
};


// UpdateSpin template specialisations
template<class lattice, class polymer>
struct UpdateSpinImpl
{
	static inline void _(lattice* lat, polymer* pol, unsigned long long* acceptCount)
	{
		double dE;
			
		lat->TrialMove(&dE);
		
		double dEcpl = lat->GetCouplingEnergy(pol->hetTable);
		bool acceptMove = MetropolisMove(lat, dE+dEcpl);

		if ( acceptMove )
		{
			lat->AcceptMove();
			++(*acceptCount);
		}
	}
};

template<class polymer>
struct UpdateSpinImpl<MCLattice, polymer>
{
	static inline void _(MCLattice*, polymer*, unsigned long long*) {}
};


// Wrapper functions
template<class lattice, class polymer>
inline void UpdateTAD(lattice* lat, polymer* pol, unsigned long long* acceptCount,unsigned long long* acceptCountTopo)
{
	UpdateTADImpl<lattice, polymer>::_(lat, pol, acceptCount,acceptCountTopo);
}

template<class lattice, class polymer>
inline void UpdateSpin(lattice* lat, polymer* pol, unsigned long long* acceptCount)
{
	UpdateSpinImpl<lattice, polymer>::_(lat, pol, acceptCount);
}


#endif /* MCUpdater_hpp */
