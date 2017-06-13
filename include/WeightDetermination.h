#ifndef SOFT_COMP_WEIGHT_DETERMINATION_H
#define SOFT_COMP_WEIGHT_DETERMINATION_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>

#include "TApplication.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TInterpreter.h"
#include "TLegend.h"
#include "TList.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"

#include "EventClass.h"

typedef std::vector<EventClass*> EventVector;

class SoftCompWeightDetermination
{
    private:
        std::string     m_rootFiles;
        int             m_numberOfEvents;
        float           m_energyDensityFinalBin;
        EventVector     m_eventVector;
        FloatVector     m_densityBins;
        IntVector       m_energies;

    public:
        /*
         * Default Constructor
         */
        SoftCompWeightDetermination(std::string rootFiles);

        /*
         * Default Destructor
         */
        ~SoftCompWeightDetermination();

        /*
         * Read data 
         */
        void LoadEvents();

        /*
         * Calculate mean bin positions for software compensation weight binning 
         */
        void MakeBinDensities();

        /*
         * Perform TMinuit fit to find software compensation weights
         */
        void Fit();

        /*
         * Work out chi squared for fit
         */
        double Chi2(const double *par);

        /*
         * Return the software compensated energy for the event given using the given parameters
         */
        float SoftwareCompensatedEnergy(EventClass *pEventClass, const double *par) const;

        /**
         * Find the energy density bin for the hit of given unbinned energy density 
         */
        float FindDensity(float energyDensity) const;


        template <class T>
        std::string NumberToString(T Number) const;
};

#endif
