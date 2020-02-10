#include "initial_state.h"
#include <itensor/all.h>
#include <iostream>

using namespace std;
using namespace itensor;

//----------------------------------------------------------------------
//all spins UP along X

void 
initial_state_all_UP( const SpinHalf sites , MPS* psi , const int N)
	{
	for(int i=1; i<=N; i++)
		{
		auto si = sites(i);
		auto wf = ITensor(si);
		
		wf.set(si(1), 1/sqrt(2));
		wf.set(si(2), 1/sqrt(2));
	
		(*psi).setA(i,wf);
		}
	}
	
//----------------------------------------------------------------------
//all spins DOWN along X

void 
initial_state_all_DOWN( const SpinHalf sites , MPS* psi , const int N)
	{
	for(int i=1; i<=N; i++)
		{
		auto si = sites(i);
		auto wf = ITensor(si);
		
		wf.set(si(1), 1/sqrt(2));
		wf.set(si(2), -1/sqrt(2));
	
		(*psi).setA(i,wf);
		}
	}
	
//----------------------------------------------------------------------
//half chain DOWN along X, half chain UP along X (domain wall with single kink)
	
void
initial_state_DOMAIN_WALL( const SpinHalf sites , MPS* psi , const int N)
	{
	for(int i=1; i<=N; i++)
		{
		auto si = sites(i);
		auto wf = ITensor(si);
		if( i <= N/2 )
			{
			wf.set(si(1), 1/sqrt(2));
			wf.set(si(2), 1/sqrt(2));
			}
		else
			{
			wf.set(si(1), 1/sqrt(2));
			wf.set(si(2), -1/sqrt(2));
			}
		(*psi).setA(i,wf);
		}
	}
	

	

	

	
