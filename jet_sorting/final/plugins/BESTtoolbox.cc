//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// BESTtoolbox.cc ------------------------------------------------------------------------
//========================================================================================
// C++ file containing functions for use with CMS EDAnalyzer and EDProducer --------------
//////////////////////////////////////////////////////////////////////////////////////////

#include "/uscms_data/d3/cannaert/analysis/CMSSW_10_6_29/src/BESTtoolbox.h"

//========================================================================================
// Calculate Legendre Polynomials --------------------------------------------------------
//----------------------------------------------------------------------------------------
// Simple Legendre polynomial function that can calculate up to order 4 ------------------
// Inputs: argument of the polynomial and order desired ----------------------------------
//----------------------------------------------------------------------------------------

float LegendreP(float x, int order){
    if (order == 0) return 1;
    else if (order == 1) return x;
    else if (order == 2) return 0.5*(3*x*x - 1);
    else if (order == 3) return 0.5*(5*x*x*x - 3*x);
    else if (order == 4) return 0.125*(35*x*x*x*x - 30*x*x + 3);
    else return 0;
}

//========================================================================================
// Calculate Fox Wolfram Moments ---------------------------------------------------------
//----------------------------------------------------------------------------------------
// This function calculates the Fox Wolfram moments for jet constituents -----------------
// in various rest frames. ---------------------------------------------------------------
// Inputs: particles (jet constiuents boosted to rest frame) and empty array that --------
//         that will store the FW moments ------------------------------------------------
//----------------------------------------------------------------------------------------

int FWMoments(std::vector<TLorentzVector> particles, double (&outputs)[5] ){

    // get number of particles to loop over
    int numParticles = particles.size();

    // get energy normalization for the FW moments
    float s = 0.0;
    for(int i = 0; i < numParticles; i++){
        s += particles[i].E();
    }

    float H0 = 0.0;
    float H4 = 0.0;
    float H3 = 0.0;
    float H2 = 0.0;
    float H1 = 0.0;

    for (int i = 0; i < numParticles; i++){

        for (int j = i; j < numParticles; j++){

            // calculate cos of jet constituent angles
            float costh = ( particles[i].Px() * particles[j].Px() + particles[i].Py() * particles[j].Py()
                                        + particles[i].Pz() * particles[j].Pz() ) / ( particles[i].P() * particles[j].P() );
            float w1 = particles[i].P();
            float w2 = particles[j].P();

            // calculate legendre polynomials of jet constiteuent angles
            float fw0 = LegendreP(costh, 0);
            float fw1 = LegendreP(costh, 1);
            float fw2 = LegendreP(costh, 2);
            float fw3 = LegendreP(costh, 3);
            float fw4 = LegendreP(costh, 4);

            // calculate the Fox Wolfram moments
            H0 += w1 * w2 * fw0;
            H1 += w1 * w2 * fw1;
            H2 += w1 * w2 * fw2;
            H3 += w1 * w2 * fw3;
            H4 += w1 * w2 * fw4;
        }
    }

    // Normalize the Fox Wolfram moments
    if (H0 == 0) H0 += 0.001;      // to prevent dividing by zero
    outputs[0] = (H0);
    outputs[1] = (H1 / H0);
    outputs[2] = (H2 / H0);
    outputs[3] = (H3 / H0);
    outputs[4] = (H4 / H0);

    return 0;
}
