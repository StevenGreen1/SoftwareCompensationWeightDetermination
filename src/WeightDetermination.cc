#include "WeightDetermination.h"

//==============================================================================

int main(int argc, char **argv)
{
    gStyle->SetOptStat(kFALSE); 
    gInterpreter->EnableAutoLoading();
    std::cout << "int main(int argc, char **argv)" << std::endl;
    std::string argument1 = argv[1];
    SoftCompWeightDetermination *pSoftCompWeightDetermination = new SoftCompWeightDetermination(argument1);
}

//==============================================================================

SoftCompWeightDetermination::SoftCompWeightDetermination(std::string rootFiles) : 
    m_rootFiles(rootFiles),
    m_numberOfEvents(0),
    m_energyDensityFinalBin(30.f)
{
    std::cout << "SoftCompWeightDetermination::SoftCompWeightDetermination" << std::endl;

    int energyArray[10] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
    IntVector energies(energyArray, energyArray + sizeof(energyArray)/sizeof(energyArray[0]));
    m_energies.insert(m_energies.end(), energies.begin(), energies.end());

    std::cout << "SoftCompWeightDetermination::MakeBinDensities" << std::endl;
    const float densityBins[10]  = {0.f, 2.f, 5.f, 7.5f, 9.5f, 13.f, 16.f, 20.f, 23.5f, 28.f};
    FloatVector densityBinsVec(densityBins, densityBins + sizeof(densityBins)/sizeof(densityBins[0]));
    m_densityBins.insert(m_densityBins.end(), densityBinsVec.begin(), densityBinsVec.end());

    this->LoadEvents();
    this->Fit();
}

//==============================================================================

SoftCompWeightDetermination::~SoftCompWeightDetermination()
{
    std::cout << "SoftCompWeightDetermination::~SoftCompWeightDetermination" << std::endl;
}

//==============================================================================

void SoftCompWeightDetermination::Fit()
{
    std::cout << "SoftCompWeightDetermination::Fit" << std::endl;
    const char *minName = "Minuit2";
    const char *algoName = "Minos";

    ROOT::Math::Minimizer* pMinimizer = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
    pMinimizer->SetMaxFunctionCalls(1000000);
    pMinimizer->SetTolerance(1.f);
    pMinimizer->SetPrintLevel(1);

    ROOT::Math::Functor functorChi2(this, &SoftCompWeightDetermination::Chi2,9);
    pMinimizer->SetFunction(functorChi2);

std::cout << "HERE" << std::endl;

    const float step = 0.000001;
    double fitPar0 = 1.0000000;
    pMinimizer->SetVariable(0, "p0", fitPar0, step);
    double fitPar1(-0.02999999);
    pMinimizer->SetVariable(1, "p1", fitPar1, step);
    double fitPar2(0.0002999999);
    pMinimizer->SetVariable(2, "p2", fitPar2, step);
    double fitPar3(-0.01999999);
    pMinimizer->SetVariable(3, "p3", fitPar3, step);
    double fitPar4(-0.0003999999);
    pMinimizer->SetVariable(4, "p4", fitPar4, step);
    double fitPar5(-0.000001999999);
    pMinimizer->SetVariable(5, "p5", fitPar5, step);
    double fitPar6(0.1999999);
    pMinimizer->SetVariable(6, "p6", fitPar6, step);
    double fitPar7(0.2999999);
    pMinimizer->SetVariable(7, "p7", fitPar7, step);
    double fitPar8(-0.0999999);
    pMinimizer->SetVariable(8, "p8", fitPar8, step);

    // do the minimization
    pMinimizer->Minimize();

    // Printouts
    const double *xs = pMinimizer->X();
    std::cout << "Minimum chi2 : " << pMinimizer->MinValue()  << std::endl;
    std::cout << "Target chi2  : " << m_numberOfEvents - 9 << std::endl;

    // expected minimum is m_numberOfEvents - #Params(9) 

    if ( this->Chi2(xs) < m_numberOfEvents - 9)
    {
        std::cout << "Minimizer " << minName << " - " << algoName << "   converged to a minimum" << std::endl;
    }

    else
    {
        std::cout << "Minimizer " << minName << " - " << algoName << "   failed to converge !!!" << std::endl;
        std::cout << "MinValue " << pMinimizer->MinValue() << std::endl;
    }

    for (int i = 0; i < 9; i++)
    {
        std::cout << "Parameter " << i << " is " << xs[i] << std::endl;
    }

    for (int i = 0; i < 9; i++)
    {
        std::cout << xs[i] << " ";
    }
    std::cout << std::endl;
}

//==============================================================================

void SoftCompWeightDetermination::LoadEvents()
{
    TApplication a("a", 0, 0);

    std::cout << "SoftCompWeightDetermination::LoadEvents" << std::endl;
    for (std::vector<int>::const_iterator it = m_energies.begin(); it != m_energies.end(); ++it)
    {
        const std::string energy = this->NumberToString(*it);
        const std::string energyString(energy + "_GeV_Kaon0L");

        std::string pathToFiles = "/r06/lc/sg568/JERDetailed/Calibration/Detector_Model_pseudo_clic_ild/Reco_Stage_76_ilcsoft_v01-17-10/SoftwareCompensationWeightTraining/RootFiles";

        TChain *pTChain = new TChain("SoftwareCompensationTrainingTree");
        TSystemDirectory directory(pathToFiles.c_str(), pathToFiles.c_str());
        TList *listOfFiles = directory.GetListOfFiles();

        if (listOfFiles)
        {
            TSystemFile *file;
            TString fileCandidate;
            TIter next(listOfFiles);

            while ((file=(TSystemFile*)next()))
            {
                fileCandidate = file->GetName();

                if (!file->IsDirectory() and fileCandidate.EndsWith("root") and fileCandidate.Contains(energyString.c_str()) and pTChain->GetEntries() < 10000)
                {
                    TString rootFileToAdd = pathToFiles + "/" + fileCandidate.Data();
                    pTChain->Add(rootFileToAdd.Data());
                }
            }
        }

        float mcEnergy(*it);
        float pfoEnergy(-1.f);
        float rawClusterEnergy(-1.f);
        FloatVector *pHitEnergies(NULL);
        FloatVector *pCellSize0(NULL);
        FloatVector *pCellSize1(NULL);
        FloatVector *pCellThickness(NULL);
        IntVector *pHitType(NULL);

        pTChain->SetBranchAddress("EnergyOfPfo",&pfoEnergy);
        pTChain->SetBranchAddress("RawEnergyOfCluster",&rawClusterEnergy);
        pTChain->SetBranchAddress("HitEnergies",&pHitEnergies);
        pTChain->SetBranchAddress("CellSize0",&pCellSize0);
        pTChain->SetBranchAddress("CellSize1",&pCellSize1);
        pTChain->SetBranchAddress("CellThickness",&pCellThickness);
        pTChain->SetBranchAddress("HitType",&pHitType);

        if (0 == pTChain->GetEntries()) continue;

        m_numberOfEvents += pTChain->GetEntries();

        for (unsigned int i = 0; i < pTChain->GetEntries(); i++)
        {
            pTChain->GetEntry(i);
//            std::cout << "mcEnergy, pfoEnergy, rawClusterEnergy, pHitEnergies, pCellSize0, pCellSize1, pCellThickness, pHitType" << std::endl;
//            std::cout << mcEnergy << " " << pfoEnergy << " " << rawClusterEnergy << " " << pHitEnergies->at(0) << " " << pCellSize0->at(0) << " " << pCellSize1->at(0) << " " << pCellThickness->at(0) << " " << pHitType->at(0) << std::endl;
            EventClass *pEventClass = new EventClass(mcEnergy, pfoEnergy, rawClusterEnergy, pHitEnergies, pCellSize0, pCellSize1, pCellThickness, pHitType);
            m_eventVector.push_back(pEventClass);
        }
    }
}

//==============================================================================

double SoftCompWeightDetermination::Chi2(const double *par)
{
//    std::cout << "SoftCompWeightDetermination::Chi2" << std::endl;
    double chi2(0.f);
   
    for (EventVector::iterator it = m_eventVector.begin(); it != m_eventVector.end(); it++)
    {
        EventClass *pEventClass = *it;
        const float softCompEnergy = this->SoftwareCompensatedEnergy(pEventClass, par);
        const float mcEnergy = pEventClass->GetMCEnergy();
        float diff = (softCompEnergy - mcEnergy)*(softCompEnergy - mcEnergy)/(0.5*mcEnergy);
        chi2 += diff;
    }
    return chi2;
}
//==============================================================================

float SoftCompWeightDetermination::SoftwareCompensatedEnergy(EventClass *pEventClass, const double *par) const 
{
//    std::cout << "SoftCompWeightDetermination::SoftwareCompensatedEnergy" << std::endl;
    float softCompEnergy(0.f);

    const float rawClusterEnergy(pEventClass->GetRawClusterEnergy());
    const float p1 = par[0] + par[1]*rawClusterEnergy + par[2]*rawClusterEnergy*rawClusterEnergy;
    const float p2 = par[3] + par[4]*rawClusterEnergy + par[5]*rawClusterEnergy*rawClusterEnergy;
    const float p3 = par[6] / (par[7] + exp(par[8]*rawClusterEnergy));

    softCompEnergy += pEventClass->GetECalEnergy();
    FloatVector hcalHitEnergy = pEventClass->GetHCalHitEnergies();
    FloatVector hcalHitEnergyDensities = pEventClass->GetHCalHitEnergyDensities();

    for (FloatVector::iterator it = hcalHitEnergyDensities.begin(); it != hcalHitEnergyDensities.end(); it++) 
    {
        const float hitEnergy(hcalHitEnergy.at(it-hcalHitEnergyDensities.begin()));
        const float binnedEnergyDensity(this->FindDensity(*it));
        const float weight = p1*exp(p2*binnedEnergyDensity) + p3;
        softCompEnergy += weight * hitEnergy;
    }
    return softCompEnergy; 
}

//==============================================================================

float SoftCompWeightDetermination::FindDensity(float hitEnergyDensity) const
{
    if (hitEnergyDensity >= m_densityBins.back())
    {
        return m_energyDensityFinalBin;
    }
    else
    {
        for (unsigned int iBin = 0; iBin < m_densityBins.size() - 1; iBin++)
        {
            const float lowBinContent = m_densityBins.at(iBin);
            const float highBinContent = m_densityBins.at(iBin+1);

            if (hitEnergyDensity >= lowBinContent && hitEnergyDensity < highBinContent)
            {
                return ((lowBinContent + highBinContent) * 0.5f);
            }
        }
    }
}

//==============================================================================

template <class T>
std::string SoftCompWeightDetermination::NumberToString ( T Number ) const 
{
    std::ostringstream ss;
    ss << Number;
    return ss.str();
}

//==============================================================================
