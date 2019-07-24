pwd
#For CMSSW_9_4_6
export SCRAM_ARCH=slc6_amd64_gcc630
cmsrel CMSSW_9_4_13
cd CMSSW_9_4_13/src
cmsenv

git cms-init

git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP

git cms-merge-topic guitargeek:ElectronID_MVA2017_940pre3

git cms-merge-topic cms-met:METFixEE2017_949_v2

#For BadMuon Filters (default code is reversed logic)
git cms-addpkg RecoMET/METFilters

#For DelPanj
git clone git@github.com:tiwariPC/DMAnaRun2.git DelPanj
cd DelPanj
git checkout 92X_2017data_deepCSV_genMet_lepID
cd -
cp -p DelPanj/tempfix/BadGlobalMuonTagger.cc RecoMET/METFilters/plugins/BadGlobalMuonTagger.cc

#For jetToolBox
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox
cd JMEAnalysis/JetToolbox
git checkout jetToolbox_94X_v3
cd -

#For DeepDoubleX
git cms-merge-topic 25371
git cms-addpkg RecoBTag/Combined
cd RecoBTag/Combined/
git clone -b V01-01-01 --depth 1 --no-checkout https://github.com/cms-data/RecoBTag-Combined.git data
cd data
git config core.sparseCheckout true
echo 'DeepDoubleX/94X/V01/' > .git/info/sparse-checkout
git checkout --
cd $CMSSW_BASE/src


#Compile (due to the external packages, will take about 15-20 mins)
scramv1 b clean
scramv1 b -j 10

# Add the area containing the MVA weights (from cms-data, to appear in "external").
# Note: the "external" area appears after "scram build" is run at least once, as above
#
cd $CMSSW_BASE/external
# below, you may have a different architecture, this is just one example from lxplus
cd slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data
cd data/RecoEgamma/PhotonIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/external
cd slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
# Go back to the src/
cd $CMSSW_BASE/src

scram b -j8

cd $CMSSW_BASE/external/slc6_amd64_gcc630/data/RecoEgamma/PhotonIdentification/data/
rm -rf Spring15/ Spring16/

cd $CMSSW_BASE/external/slc6_amd64_gcc630/data/RecoEgamma/ElectronIdentification/data/
rm -rf Spring15/ Spring16_GeneralPurpose_V1/ Spring16_HZZ_V1/

cd $CMSSW_BASE/src

## need to copy them so that the code can be tested in test directory where original config files are placed. 
cp DelPanj/CrabUtilities/effArea* DelPanj/TreeMaker/test/RunConfigTest/
