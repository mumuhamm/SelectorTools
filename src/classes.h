#include "Analysis/SelectorTools/interface/FakeRateSelector.h"
#include "Analysis/SelectorTools/interface/SelectorBase.h"
#include "Analysis/SelectorTools/interface/ZSelector.h"
#include "Analysis/SelectorTools/interface/LowPileupSelector.h"
#include "Analysis/SelectorTools/interface/LowPileupZSelector.h"
#include "Analysis/SelectorTools/interface/LowPileupWSelector.h"
#include "Analysis/SelectorTools/interface/LowPileupWBackgroundSelector.h"
#include "Analysis/SelectorTools/interface/ZZGenSelector.h"
#include "Analysis/SelectorTools/interface/WGenSelector.h"
#include "Analysis/SelectorTools/interface/ZGenSelector.h"
#include "Analysis/SelectorTools/interface/NanoGenSelectorBase.h"
#include "Analysis/SelectorTools/interface/WZSelector.h"
#include "Analysis/SelectorTools/interface/TTTSelector.h"
#include "Analysis/SelectorTools/interface/ThreeLepSelector.h"
#include "Analysis/SelectorTools/interface/WZSelectorBase.h"
#include "Analysis/SelectorTools/interface/WZBackgroundSelector.h"
#include "Analysis/SelectorTools/interface/ScaleFactor.h"
#include "Analysis/SelectorTools/interface/disambiguateFinalStates.h"
#include "Analysis/SelectorTools/interface/disambiguateFinalStatesZZ.h"
#include "Analysis/SelectorTools/interface/Efficiency.h"


namespace{
  namespace{
    FakeRateSelector pFakeRateSelector;
    WZSelectorBase pWZSelectorBase;
    SelectorBase pSelectorBase;
    ZSelector pZSelector;
    LowPileupZSelector pLowPileupZSelector;
    LowPileupWSelector pLowPileupWSelector;
    LowPileupWSelector pLowPileupWBackgroundSelector;
    LowPileupSelector pLowPileupSelector;
    WZSelector pWZSelector;
    NanoGenSelectorBase pNanoGenSelectorBase;
    ZZGenSelector pZZGenSelector;
    TTTSelector pTTTSelector;
    ThreeLepSelector pThreeLepSelector;
    WGenSelector pWGenSelector;
    ZGenSelector pZGenSelector;
    WZBackgroundSelector pWZBackgroundSelector;
    ScaleFactor pScaleFactor;
    disambiguateFinalStates pDisambiguator;
    disambiguateFinalStates pDisambiguatorZZ;
    Efficiency pEfficiency;
  }
}
