#include "Analysis/VVAnalysis/interface/SelectorBase.h"
#include <boost/algorithm/string.hpp>
#include <TStyle.h>
#include <regex>
#include "TParameter.h"

void SelectorBase::Begin(TTree * /*tree*/)
{
    TString option = GetOption();
}

void SelectorBase::SlaveBegin(TTree * /*tree*/)
{
    if (GetInputList() != nullptr) {
        TParameter<bool>* applyScaleFactors = (TParameter<bool>*) GetInputList()->FindObject("applyScaleFacs");
        if (applyScaleFactors != nullptr && applyScaleFactors->GetVal()) {
           SetScaleFactors();
        }
    }
}

void SelectorBase::Init(TTree *tree)
{
    if (!tree) return;
    fChain = tree;
    
    TString option = GetOption();

    if (GetInputList() != nullptr) {
	TNamed* ntupleType = (TNamed *) GetInputList()->FindObject("ntupleType");
    TNamed* name = (TNamed *) GetInputList()->FindObject("name");
    TNamed* chan = (TNamed *) GetInputList()->FindObject("channel");
    TNamed* selection = (TNamed *) GetInputList()->FindObject("selection");
	TNamed* year = (TNamed *) GetInputList()->FindObject("year");

	if (ntupleType != nullptr) {
	    std::string ntupleName = ntupleType->GetTitle();
	    if (ntupleName == "NanoAOD")
		    ntupleType_ = NanoAOD;
	    else if (ntupleName  == "UWVV")
		ntupleType_ = UWVV;
	    else if (ntupleName  == "Bacon")
		ntupleType_ = Bacon;
	    else
		    throw std::invalid_argument("Unsupported ntuple type!");
	}
	else {
	    std::cerr << "INFO: Assuming NanoAOD ntuples" << std::endl;
	    ntupleType_ = NanoAOD;
	}

        if (name != nullptr) {
            name_ = name->GetTitle();
        }
        else {
            name_ = GetNameFromFile();
        }
        if (name_.empty()){
            std::cerr << "INFO: Using default name \"Unknown\" for file" << std::endl;
            name_ = "Unknown";
        }
	if(year != nullptr) {
	    year_ = yearMap_[year->GetTitle()];
	}
	
	if (chan != nullptr) {
	    channelName_ = chan->GetTitle();
	}
	else if (ntupleType_ == UWVV)
            channelName_ = fChain->GetTree()->GetDirectory()->GetName();
        if (selection != nullptr) {
            selectionName_ = selection->GetTitle();
        }
    }
    
    if (selectionMap_.find(selectionName_) != selectionMap_.end()) {
	    selection_ = selectionMap_[selectionName_];
    }
    else
        throw std::invalid_argument(selectionName_ + " is not a valid selection!");
    
    isMC_ = false;
    if (name_.find("data") == std::string::npos){
        isMC_ = true;
    }
    if (doSystematics_ && isMC_ && !isNonprompt_)
        variations_.insert(systematics_.begin(), systematics_.end());

    if (channelMap_.find(channelName_) != channelMap_.end())
        channel_ = channelMap_[channelName_];
    else {
        std::string message = "Invalid channel choice! ";
        message += "Choice was " + channelName_ + "\n";
        message += "Valid choices are: ";
        for (const auto& chan : channelMap_)
            message += chan.first + ", " ;
        throw std::invalid_argument(message);
    }

    // only make the directory iff class isn't being run as a slave class /////
    TNamed* isSlaveClass = (TNamed *) GetInputList()->FindObject("isSlaveClass");
    if(isSlaveClass != nullptr) return;
    
    makeOutputDirs();
    SetBranches();
}

void SelectorBase::addSubprocesses(std::vector<std::string> processes) {
    subprocesses_ = processes;
    for (auto& proc : subprocesses_)
        std::cout << "Subproces " << proc << std::endl;
}

void SelectorBase::setSubprocesses(std::string process) {
    currentHistDir_ = dynamic_cast<TList*>(fOutput->FindObject(process.c_str()));
    if (currentHistDir_ == nullptr) {
        currentHistDir_ = new TList();
        currentHistDir_->SetName(process.c_str());
        std::cout << "Setting output dir name to " << process << std::endl;
        fOutput->Add(currentHistDir_);
    }
        //throw std::invalid_argument(process + " is not a valid subprocess for selector!");
}

void SelectorBase::makeOutputDirs() {
    std::vector<std::string> allProcesses = {name_};
    allProcesses.insert(allProcesses.begin(), subprocesses_.begin(), subprocesses_.end());
    
    for (auto& name : allProcesses) { 
        currentHistDir_ = dynamic_cast<TList*>(fOutput->FindObject(name.c_str()));
        
        if ( currentHistDir_ == nullptr ) {
            currentHistDir_ = new TList();
            currentHistDir_->SetName(name.c_str());
            fOutput->Add(currentHistDir_);
            if (name == name_)
                SetupNewDirectory();
        }
    }
    //TODO: Determine why this fails and if it's really needed
    //size_t existingObjectPtrsSize = allObjects_.size();
    // Watch for something that I hope never happens,
    //if ( existingObjectPtrsSize > 0 && allObjects_.size() != existingObjectPtrsSize ) {
    //    throw std::invalid_argument(Form("SelectorBase: Size of allObjects has changed!: %lu to %lu", existingObjectPtrsSize, allObjects_.size()));
    //}
    //UpdateDirectory();
}

void SelectorBase::SetScaleFactors() {
    std::invalid_argument("No scale factors defined for selector!");
}

void SelectorBase::SetBranches() {
    if (ntupleType_ == UWVV)
        SetBranchesUWVV();
    else if (ntupleType_ == NanoAOD)
        SetBranchesNanoAOD();
    else if (ntupleType_ == Bacon)
        SetBranchesBacon();
    else 
        throw std::invalid_argument("Undefined ntuple type!");
        
}

void SelectorBase::LoadBranches(Long64_t entry, SystPair variation) {
    if (ntupleType_ == NanoAOD)
        LoadBranchesNanoAOD(entry, variation);
    else if (ntupleType_ == Bacon)
        LoadBranchesBacon(entry, variation);
    else {
        LoadBranchesUWVV(entry, variation);
    }
}

Bool_t SelectorBase::Process(Long64_t entry)
{
    for (const auto& variation : variations_) {
        LoadBranches(entry, variation);
        FillHistograms(entry, variation);
    }
    return kTRUE;
}

Bool_t SelectorBase::Notify()
{
    return kTRUE;
}

float SelectorBase::GetPrefiringEfficiencyWeight(
        std::vector<float>* jetPt, std::vector<float>* jetEta) {
    float prefire_weight = 1;
    for (size_t i = 0; i < jetPt->size(); i++) {
        float jPt = jetPt->at(i);
        float jEta = std::abs(jetEta->at(i));
        prefire_weight *= (1-prefireEff_->GetEfficiency(prefireEff_->FindFixBin(jEta, jPt)));
    }
    return prefire_weight;
}

void SelectorBase::Terminate()
{
}
    
void SelectorBase::SlaveTerminate()
{
}
void SelectorBase::UpdateDirectory()
{
  for(TNamed** objPtrPtr : allObjects_) {
    if ( *objPtrPtr == nullptr ) std::invalid_argument("SelectorBase: Call to UpdateObject but existing pointer is null");
    *objPtrPtr = (TNamed *) currentHistDir_->FindObject((*objPtrPtr)->GetName());
    if ( *objPtrPtr == nullptr ) std::invalid_argument("SelectorBase: Call to UpdateObject but current directory has no instance");
  }
}

template<typename T>
void SelectorBase::InitializeHistMap(std::vector<std::string>& labels, std::unordered_map<HistLabel, T*>& histMap) {
    for (auto& label : labels) {
        if (channel_ != Inclusive) {
            HistLabel histlabel = {label, channel_, Central};
            histMap[histlabel] = {};
        }
        else {
            for (auto& chan : allChannels_) {
                HistLabel histlabel = {label, chan.first, Central};
                histMap[histlabel] = {};
            }
        }
    }
}

void SelectorBase::InitializeHistogramsFromConfig() {
    TList* histInfo = (TList *) GetInputList()->FindObject("histinfo");
    if (histInfo == nullptr ) 
        throw std::domain_error("Can't initialize histograms without passing histogram information to TSelector");

    InitializeHistMap(hists1D_,histMap1D_);
    InitializeHistMap(hists2D_,histMap2D_);
    InitializeHistMap(weighthists1D_, weighthistMap1D_);
    InitializeHistMap(weighthists2D_, weighthistMap2D_);

    std::vector<std::string> tempSystHistNames;
    for (auto& syst : systematics_) {
        for (auto hist : systHists_) {
            tempSystHistNames.push_back(hist + "_" + syst.second);
        }
        for (auto hist : systHists2D_) {
            tempSystHistNames.push_back(hist + "_" + syst.second);
        }
    }

    for (auto && entry : *histInfo) {  
        TNamed* currentHistInfo = dynamic_cast<TNamed*>(entry);
        const char* name = currentHistInfo->GetName();
        std::vector<std::string> histData = ReadHistDataFromConfig(currentHistInfo->GetTitle());
        
        std::vector<ChannelPair> channels = {{channel_, channelName_}};
        if (channel_ == Inclusive) {
            channels = allChannels_;
        }

        for (const auto& chan : channels) {
            HistLabel centralLabel = {name, chan.first, Central};
            if (histMap1D_.find(centralLabel) != histMap1D_.end() || histMap2D_.find(centralLabel) != histMap2D_.end()) { 
                InitializeHistogramFromConfig(name, chan, histData);
            }
            // May add a syst hist to plotobjects for plotting purposes, surpress errors in this case
            else if (std::find(tempSystHistNames.begin(), tempSystHistNames.end(), name) == tempSystHistNames.end())
                std::cerr << "Skipping invalid histogram '" << name << "'" << std::endl;
        }
    }

    for (auto& subprocess : subprocesses_) {
        setSubprocesses(subprocess);
        auto& subprocessMap = subprocessHistMaps1D_[subprocess];
        subprocessMap = {};
        for (auto& hist : histMap1D_) {
            subprocessMap[hist.first] = {};
            if (hist.second == nullptr)
                continue;
            AddObject<TH1D>(subprocessMap[hist.first], const_cast<TH1D&>(*hist.second));
            //subprocessWeightHistMaps1D_[subporcess] = weighthistMap1D_;
        }
    }

    currentHistDir_ = dynamic_cast<TList*>(fOutput->FindObject(name_.c_str()));
}

void SelectorBase::InitializeHistogramFromConfig(std::string name, ChannelPair channel, std::vector<std::string>& histData) {
    if (histData.size() != 4 && histData.size() != 7) {
        std::cerr << "Malformed data string for histogram '" << name
                    << ".' Must have form: 'Title; (optional info) $ nbins, xmin, xmax'"
                    << "\n   OR form: 'Title; (optional info) $ nbins, xmin, xmax nbinsy ymin ymax'"
                    << std::endl;
        exit(1);
    }
    std::string histName = getHistName(name, {}, channel.second);
    
    int nbins = std::stoi(histData[1]);
    float xmin = std::stof(histData[2]);
    float xmax = std::stof(histData[3]);

    HistLabel centralLabel = {name, channel.first, Central};
    if (histData.size() == 4) {
        AddObject<TH1D>(histMap1D_[centralLabel], histName.c_str(), histData[0].c_str(),nbins, xmin, xmax);
        if (doSystematics_ && !isNonprompt_ && std::find(systHists_.begin(), systHists_.end(), name) != systHists_.end()) {
            for (auto& syst : systematics_) {
                HistLabel systLabel = {name, channel.first, syst.first};
                std::string syst_histName = getHistName(name, syst.second, channel.second);
                histMap1D_[systLabel] = {};
                AddObject<TH1D>(histMap1D_[systLabel], syst_histName.c_str(), 
                    histData[0].c_str(),nbins, xmin, xmax);
                // TODO: Cleaner way to determine if you want to store systematics for weighted entries
                //if (isaQGC_ && doaQGC_ && (weighthistMap1D_.find(name) != weighthistMap1D_.end())) { 
                //    std::string weightsyst_histName = name+"_lheWeights_"+syst.second;
                //    AddObject<TH2D>(weighthistMap1D_[syst_histName], 
                //        (weightsyst_histName+"_"+channel).c_str(), histData[0].c_str(),
                //        nbins, xmin, xmax, 1000, 0, 1000);
                //}
            }
        }
        // Weight hists must be subset of 1D hists!
        if (isMC_ && !isNonprompt_ && (weighthistMap1D_.find(centralLabel) != weighthistMap1D_.end())) { 
            AddObject<TH2D>(weighthistMap1D_[centralLabel], 
                (name+"_lheWeights_"+channel.second).c_str(), histData[0].c_str(),
                nbins, xmin, xmax, 2000, 0, 2000);
        }
    }
    else {
        int nbinsy = std::stoi(histData[4]);
        float ymin = std::stof(histData[5]);
        float ymax = std::stof(histData[6]);
        AddObject<TH2D>(histMap2D_[centralLabel], histName.c_str(), histData[0].c_str(),nbins, xmin, xmax,
                nbinsy, ymin, ymax);
        if (doSystematics_ && std::find(systHists2D_.begin(), systHists2D_.end(), name) != systHists2D_.end()) {
            for (auto& syst : systematics_) {
                HistLabel systLabel = {name, channel.first, syst.first};
                std::string syst_histName = getHistName(name, syst.second, channel.second);
                histMap2D_[systLabel] = {};
                AddObject<TH2D>(histMap2D_[systLabel], syst_histName.c_str(), 
                    histData[0].c_str(),nbins, xmin, xmax, nbinsy, ymin, ymax);
            }
        }
        // 3D weight hists must be subset of 2D hists!
        if (isMC_ && std::find(weighthists2D_.begin(), weighthists2D_.end(), name) != weighthists2D_.end()) {
            AddObject<TH3D>(weighthistMap2D_[centralLabel], 
                (name+"_lheWeights_"+channel.second).c_str(), histData[0].c_str(),
                nbins, xmin, xmax, nbinsy, ymin, ymax, 2000, 0, 2000);
        }
    }
}

std::vector<std::string> SelectorBase::ReadHistDataFromConfig(std::string histDataString) {
    std::vector<std::string> histData;
    boost::split(histData, histDataString, boost::is_any_of("$"));
    std::vector<std::string> binInfo;
    if (histData.size() != 2)
        return {};
    
    boost::split(binInfo, histData[1], boost::is_any_of(","));
   
    histData.pop_back();
    for (const auto& x : binInfo) {
        histData.push_back(x);
    }
    
    return histData;
}

void SelectorBase::SetupNewDirectory()
{
}

std::string SelectorBase::concatenateNames(const std::string& baseName, std::string& toAppend) {
    return concatenateNames(baseName.c_str(), toAppend);
}

std::string SelectorBase::concatenateNames(const std::string& baseName, const char* toAppend) {
    std::string app = toAppend;
    return concatenateNames(baseName.c_str(), toAppend);
}

std::string SelectorBase::concatenateNames(const char* baseName, std::string& toAppend) {
    if (toAppend.empty())
        return baseName;
    const std::string delimit = "_";
    std::string name = baseName;
    name.reserve(name.size()+delimit.size()+toAppend.size());
    name.append(delimit);
    name.append(toAppend);
    return name;
}

std::string SelectorBase::getHistName(std::string histName, std::string variationName) {
    return getHistName(histName, variationName, channelName_);
}

std::string SelectorBase::getHistName(std::string histName, std::string variationName, std::string channel) {
    histName.append("_");
    if (!variationName.empty()) {
        histName.append(variationName);
        histName.append("_");
    }
    histName.append(channel);
    return histName;
}
