#ifndef __SORTJETS_H__
#define __SORTJETS_H__

#include <cmath>
#include <iostream>
#include "TLorentzVector.h"
#include <string>

using namespace std;
#include <iostream>

class sortJets
{	
	private:
		
		void findBestJetComb(void);

		std::vector<TLorentzVector> superJet1;
		std::vector<TLorentzVector> superJet2;
		std::vector<TLorentzVector> miscJets;
		int nMiscJets;
	public:
		sortJets(std::vector<TLorentzVector>,std::vector<TLorentzVector>,std::vector<TLorentzVector>);
		std::vector<TLorentzVector> finalSuperJet1;
		std::vector<TLorentzVector> finalSuperJet2;
		int nCombinations = 0;
		//int nCombinationsArraySize = pow(nMiscJets,2);  //THIS DOES NOT WORK IN C++
		//double massCombinations[1028][2];    //defined 14 miscjet cases, this would be an array of length 16k+, don't want to lug this around ... 

};	
#endif
