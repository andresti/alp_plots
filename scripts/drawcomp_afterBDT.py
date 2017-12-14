# to compare histograms - max samples with different options
# python scripts/drawcomp_afterBDT.py -w 0 -b 20170116-181049

import json
import os
from glob import glob
import importlib

# ROOT imports
from ROOT import TChain, TH1F, TFile, vector, TCanvas, gROOT

from Analysis.alp_analysis.samplelists import samlists
from Analysis.alp_analysis.alpSamples import samples
from Analysis.alp_plots.histOpt import hist_opt
import Analysis.alp_plots.UtilsDraw as UtilsDraw

TH1F.AddDirectory(0)
gROOT.SetBatch(True)

# parsing parameters
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-w", "--whichPlots", help="which plots to be produced", type=int, default='0')
parser.add_argument("-b", "--bdt", help="bdt version, equal to input file name (classifier report output)", required=True)
parser.add_argument("-o", "--oDir", help="output directory", default="plots")
parser.add_argument("--res", dest="plotResidual", help="to plot residuals (2==pulls)", type=int, default=0)
parser.add_argument("-r", "--clrebin", help="to rebin (classifier output)", type=int, default=-1)
parser.add_argument("-n", "--doNorm", help="do normalize", action='store_true')
parser.add_argument("-c", "--defaultCol", help="to use default colors", action='store_true')
parser.add_argument("-l", "--list", help="hist list", dest="hlist", type=int, default=1)
parser.add_argument("--lumi", help="int lumi to normalize to", dest="lumi", type=float, default=35.9*0.96) #0.96 trg eff for MC
parser.set_defaults(doNorm=False, defaultCol=False)
args = parser.parse_args()

iDir       = '../hh2bbbb_limit/'
filename = iDir+"/"+args.bdt+".root"
headerOpt = args.bdt
intLumi_fb = args.lumi # plots normalized to this
which = args.whichPlots
getChi = False
getVar = False
drawH2 = False
weights = [[],[]]
sf = [[],[]]

if args.hlist == 0:
    histList   = [ 'classifier']
elif args.hlist == 1: # exact list of BDT input variables
    histList   = [
                  'h_jets_ht', 'h_jets_ht_r',
                  'h_jet0_pt', 'h_jet1_pt', 'h_jet2_pt', 'h_jet3_pt',
                  'h_jet0_eta', 'h_jet1_eta', 'h_jet2_eta', 'h_jet3_eta',
                  'h_cmva3','h_cmva4',
                  'h_H0_mass','h_H0_pt','h_H0_csthst0_a','h_H0_dr','h_H0_dphi',
                  'h_H1_mass','h_H1_pt', 'h_H1_dr','h_H1_dphi',
                  'h_H0H1_mass', 'h_H0H1_pt', 'h_H0H1_csthst0_a',
                  'h_X_mass'
                  'classifier'
                 ]
elif args.hlist == 2:
    histList   = [
                  'h_csv3','h_csv4',
                  'h_H1_csthst2_a', 'h_H1_dr','h_H1_dphi','h_H0H1_dr',

                 ]

histList2  = ["DiJets[0].mass()-DiJets[1].mass()", "CSV_Jet2-CSV_Jet3", "CMVA_Jet2-CMVA_Jet3",] # -- not maintained

###############
if which == -3:
    samples = [['tt','qcd_m'], ['bkg']]
    fractions = ['','test']
    regions = ['','']
    weights = [[0.010760,17.635231,3.476259,0.509821,0.089584,0.046491,0.005935,0.002410],[]]
    legList = [['tt','qcd HT>200'], ["bkg (mixed data)"]]
    colorList = [604, 430]
    dofill = [True,True]
    isMC = True
    oname = 'comp_qcdttBkg_afterBDT'
    headerOpt = ""

elif which == -2:
    samples = [['HHTo4B_SM'], ['sig']]
    fractions = ['','train']
    regions = ['','']
    legList = [["ggHH4b - SM"], ["ggHH4b - PangeaSM"]]
    colorList = [[630], [633]]
    sf = [[7.335],[1.]] #7.335 2.439
    dofill = [True,True]
    isMC = True
    oname = 'comp_SMpangeaSM_afterBDT'
    headerOpt = ""

elif which == -21:
    samples = [['HHTo4B_BM2'], ['sig']]
    fractions = ['','appl']
    regions = ['','']
    legList = [["ggHH4b - BM2"], ["ggHH4b - PangeaBM2"]]
    colorList = [[630], [633]]
    sf = [[2.782],[1.]] #7.335 
    dofill = [True,True]
    isMC = True
    oname = 'comp_BM2pangeaBM2_afterBDT'
    headerOpt = ""

elif which == -1:
    samples = [['pan'], ['sig']]
    fractions = ['','']
    regions = ['','']
    legList = [["ggHH4b - Pangea"], ["ggHH4b - PangeaSM"]]
    colorList = [630, 633]
    dofill = [True,True]
    isMC = True
    oname = 'comp_SMpangea_afterBDT'
    headerOpt = ""

#elif which == -1:
 #   samples = [['sig'], ['sample']]
  ##  fractions = ['test','']
   # regions = ['sr','sr']
  #  legList = [["ggHH4b - pangeaSM"],["ggHH4b - SM"]]
  #  colorList = [633, 630]
  #  dofill = [True,True]
  #  isMC = True
  #  oname = 'comp_'+samples[0]+samples[1]+'_afterBDT'
  #  headerOpt = "SR"

elif which == 0:
    samples = [['sig'], ['bkg']]
    fractions = ['appl','appl']
    regions = ['','']
    legList = [["ggHH4b SM"], ["bkg (mixed data)"]]
    colorList = [[632], [430]]
    dofill = [True,True]
    sf = [[11.8756],[0.25]] #(33.53*0.5824*0.5824/(4172119.0*0.2))
    isMC = True
    oname = 'comp_sigBkg_afterBDT'
    headerOpt = "    appl samples" #{}".format()

elif which == 1:
    samples = [['sig'], ['bkg']]
    fractions = ['train','train']
    regions = ['','']
    legList = [["ggHH4b BM12"], ["bkg (mixed data)"]]
    colorList = [[632], [430]]
    dofill = [True,True]
    sf = [[1.],[0.25]] #(33.53*0.5824*0.5824/(4172119.0*0.2)) , 3.968109
    isMC = True
    oname = 'comp_sigBkg_afterBDT'
    headerOpt = "    train samples" #{}".format()

elif which == 2:
    samples = [['sig'], ['bkg']]
    fractions = ['test','test']
    regions = ['', '']
    legList = [["ggHH4b SM"], ["bkg (mixed data)"]]
    colorList = [[632], [430]]
    dofill = [True,True]
    sf = [[11.99561921],[0.25]] #(33.53*0.5824*0.5824/(4172119.0*0.2))
    isMC = True
    oname = 'comp_sigBkg_afterBDT'
    headerOpt = "    test samples" #{}".format()

elif which == 3:
    samples = [['sig'], ['bkg']]
    fractions = ['test','test']
    regions = ['cr', 'cr']
    legList = [["signal (HH4b SM) - CR"], ["bkg (mixed data) - CR"]]
    colorList = [632, 430]
    dofill = [True,True]
    isMC = False
    oname = 'comp_'+samples[0]+samples[1]+'_afterBDT'

elif which == 4:
    samples = [['sig'], ['bkg']]
    fractions = ['train','train']
    regions = ['sr', 'sr']
    legList = [["signal (HH4b SM) - SR, train fract"], ["bkg (mixed data) - SR, train fract"]]
    colorList = [632, 430]
    dofill = [True,True]
    oname = 'comp_'+samples[0]+samples[1]+'_afterBDT'

elif which == 5:
    samples = [['bkg'],['data']] #data always  second
    fractions = ['appl','']
    regions = ['ms','ms']
    legList = [["mixed data"], ["data"]]
    colorList = [[430], [1]]
    sf = [[0.25],[1.]]
    dofill = [True,False]
    isMC = False
    oname = 'comp_bkgdata_afterBDT'
    headerOpt = "   mass CR. appl sample"

elif which == 6:
    samples = [['bkg'],['bkg']]
    fractions = ['test','test']
    regions = ['','sr']
    legList = [["bkg (mixed data)"], ["data - SR"]]
    colorList = [430, 416]
    dofill = [True,True]
    isMC = False
    oname = 'comp_'+samples[0]+samples[1]+'_afterBDT'

elif which == 7:
    samples = [['sig'], ['sig']]
    fractions = ['train','test']
    regions = ['', '']
    legList = [["signal (HH4b SM) - train"], ["signal (HH4b SM) - test"]]
    colorList = [632, 600]
    dofill = [True,True]
    isMC = True
    oname = 'comp_'+samples[0]+samples[1]+'_afterBDT'

elif which == 8:
    samples = [['bkg'], ['bkg']]
    fractions = ['appl','test']
    regions = ['', '']
    legList = [["bkg (mixed data) - appl"], ["bkg (mixed data) - test"]]
    colorList = [[430], [400]]
    dofill = [True,True]
    isMC = False
    oname = 'comp_bkgAppl-bkgTest_afterBDT'

#elif which == 9:
#    samples = [['bkg'],['data']] #data always  second
#    fractions = ['test','']
#    regions = ['msbdt','msbdt'] 
#    legList = [["bkg (mixed data)"], ["data"]]
#    colorList = [430, 1]
#    sf = [[1.],[1.]]
#    dofill = [True,False]
#    isMC = False
#    oname = 'comp_bkgdata_afterBDT'
#    headerOpt = " BDT[0.8-1] & H1-H2 mass blinded" #h1-h2 mass cut

elif which == 9:
    samples = [['sig'],['bkg']] #data always  second
    fractions = ['test','test']
    regions = ['ms','ms'] 
    legList = [["ggHH4b SM"], ["mixed data"]]
    colorList = [[632], [430]]
    sf = [[(33.53*0.5824*0.5824/(4172119.0*0.2))],[0.25]]
    dofill = [True,True]
    isMC = True
    oname = 'comp_sigbkgms_afterBDT'
    headerOpt = "  H mass CR - test sam."

elif which == 10:
    samples = [['bkg'],['data']] #data always  second
    fractions = ['appl','']
    regions = ['ms','ms'] 
    legList = [["mixed data"], ["data"]]
    colorList = [[430], [1]]
    sf = [[0.25],[1.]]
    dofill = [True,False]
    isMC = False
    oname = 'comp_bkgdata_afterBDT'
    headerOpt = "  H mass CR - appl.sam." #

elif which == 101:
    samples = [['BTagCSVRun2016-mixed-extreme2'],['data']] #data always  second
    fractions = ['','']
    regions = ['ms2','ms2'] 
    legList = [["mixed data"], ["data"]]
    colorList = [[430], [1]]
    sf = [[0.25],[1.]]
    dofill = [True,False]
    isMC = False
    oname = 'comp_bkgExtr2data_afterBDT'
    headerOpt = "  H mass CR2 - extreme2 sam." #
    #drawH2 = True

elif which == 102:
    samples = [['BTagCSVRun2016-mixed-test-htall'],['data']] #data always  second --11nn
    fractions = ['','']
    regions = ['ms','ms'] 
    legList = [["mixed data"], ["data"]]
    colorList = [[430], [1]]
    sf = [[0.25],[1.]]
    dofill = [True,False]
    isMC = False
    oname = 'comp_bkgHtalldata_afterBDT'
    headerOpt = "  H mass CR. htall" #
   # drawH2 = True

elif which == 1002:
    samples = [['BTagCSVRun2016-mixed-train-pt1234'],['BTagCSVRun2016-mixed-11']] #data always  second
    fractions = ['','']
    regions = ['ms','ms'] 
    legList = [["mixed data - new"], ["mixed data"]]
    colorList = [[430], [435]]
    sf = [[1.],[1.]]
    dofill = [True,True]
    isMC = False
    oname = 'comp_bkgNewbkg_afterBDT'
    headerOpt = " H mass CR.  n-n 1-1" #

elif which == 103:
    samples = [['bkg'],['bkg']] #data always  second
    fractions = ['appl','train']
    regions = ['',''] 
    legList = [["mixed data - appl"], ["mixed data - train"]]
    colorList = [[430], [435]]
    sf = [[0.25],[0.25]]
    dofill = [True,True]
    isMC = False
    oname = 'comp_bkgbkgApplTrain_afterBDT'
    headerOpt = "  "

elif which == 104:
    samples = [['BTagCSVRun2016-11'],['BTagCSVRun2016-11']] #bkg
    fractions = ['train','appl']
    regions = ['ms','ms'] 
    legList = [["mixed mixed data - train"], ["mixed mixed data - appl"]]
    colorList = [[430], [435]]
    sf = [[0.25],[0.25]]
    dofill = [True,True]
    isMC = False
    oname = 'comp_mixedmixed_afterBDT'
    headerOpt = "  from n-n 1-1"

elif which == 105:
    samples = [['BTagCSVRun2016-11-mixedor-11'],['BTagCSVRun2016-11']]
    fractions = ['','']
    regions = ['ms','ms'] 
    legList = [["mixed mixed data - 11"], ["mixed data"]]
    colorList = [[597], [420]]
    sf = [[1.],[1.]]
    dofill = [True,True]
    isMC = False
    oname = 'comp_mixedmixedor1111_afterBDT'
    headerOpt = "  CR. from n-n 1-1. orig lib"

elif which == 106:
    samples = [['BTagCSVRun2016applTT-mixed-appl'],['BTagCSVRun2016applTT']]
    fractions = ['','']
    regions = ['',''] 
    legList = [["mixed mixed data - appl"], ["mixed data - appl"]]
    colorList = [[597], [420]]
    sf = [[0.25],[1.]]
    dofill = [True,True]
    isMC = False
    oname = 'comp_mixedmixeddatatt_afterBDT'
    headerOpt = "  mixed data + TT (4pc)"

elif which == 1006:
    samples = [['BTagCSVRun2016appl-mixedor-appl'],['BTagCSVRun2016appl']]
    fractions = ['','']
    regions = ['ms','ms'] 
    legList = [["mixed mixed data - appl"], ["mixed data - appl"]]
    colorList = [[597], [403]]
    sf = [[0.25],[1.]]
    dofill = [True,True]
    isMC = False
    oname = 'comp_mixedmixeddataapplorappl_afterBDT'
    headerOpt = "  CR. original hems"

elif which == 107:
    samples = [['BTagCSVRun2016appl-mixed-appl'],['bkg']]
    fractions = ['','appl']
    regions = ['',''] 
    legList = [["mixed mixed data - appl"], ["mixed data"]]
    colorList = [[597], [420]]
    sf = [[0.25],[1.]]
    dofill = [False,False]
    isMC = False
    oname = 'comp_mixedmixeddataapplappl_afterBDT'
    headerOpt = "   orig: mixed data appl"

elif which == 1007:
    samples = [['BTagCSVRun2016-11-mixed-11'],['BTagCSVRun2016-mixed-11']]
    fractions = ['','']
    regions = ['',''] 
    legList = [["mixed mixed data - 11"], ["mixed data"]]
    colorList = [[634], [403]]
    sf = [[1.],[1.]]
    dofill = [False,False]
    isMC = False
    oname = 'comp_mixedmixeddata1111_afterBDT'
    headerOpt = "  orig: mixed data n-n 1-1"

elif which == 108:
    samples = [['bkg'],['data']]
    fractions = ['appl','']
    regions = ['ms','ms']
    legList = [["mixed data - appl"], ["data"]]
    colorList = [[634], [1]]
    sf = [[0.25],[1.]]
    dofill = [True,False]
    isMC = False
    oname = 'comp_mixeddataappl_afterBDT'
    headerOpt = "  CR."

elif which == 1008:
    samples = [['BTagCSVRun2016-mixed-11'],['data']]
    fractions = ['','']
    regions = ['ms','ms']
    legList = [["mixed data - 11"], ["data"]]
    colorList = [[634], [1]]
    sf = [[1.],[1.]]
    dofill = [True,False]
    isMC = False
    oname = 'comp_mixeddata11_afterBDT'
    headerOpt = "  CR. orig sample: data"

elif which == 109:
    samples = [['BTagCSVRun2016-mixed-appl'],['data']]
    fractions = ['','']
    regions = ['ms','ms']
    legList = [["mixed data - appl"], ["data"]]
    colorList = [[634], [403]]
    sf = [[0.25],[1.]]
    dofill = [True,True]
    isMC = False
    oname = 'comp_mixeddataAppl_afterBDT'
    headerOpt = "  CR"

elif which == 1009:
    samples = [['BTagCSVRun2016-mixed-11'],['data']]
    fractions = ['','']
    regions = ['ms','ms']
    legList = [["mixed data - 11"], ["data"]]
    colorList = [[634], [403]]
    sf = [[1.],[1.]]
    dofill = [True,True]
    isMC = False
    oname = 'comp_mixeddata11_afterBDT'
    headerOpt = "  CR"

elif which == 1010:
    samples = [['BTagCSVRun2016-22-mixed-appl'],['BTagCSVRun2016-22']]
    fractions = ['','']
    regions = ['',''] 
    legList = [["mixed mixed data - appl"], ["mixed data"]]
    colorList = [[634], [403]]
    sf = [[0.25],[1.]]
    dofill = [True,True]
    isMC = False
    oname = 'comp_mixedmixeddata22_afterBDT'
    headerOpt = "  mixed data n-n 2-2"

#elif which == 10:
#    samples = [['bkg'],['bkg']] #data always  second
 #   fractions = ['test','test']
#    regions = ['ms',''] 
 #   legList = [["bkg (mixed data), mCut"], ["bkg (mixed data)"]]
 #   colorList = [404, 430]
 #   dofill = [True,True]
 #   isMC = False
 #   oname = 'comp_bkgbkg__mCR_afterBDT'
 #   headerOpt = "   H1-H2 mass blinded" #h1-h2 mass cut

#antitag ---
elif which == 11:
    samples = [['amixed-0'],['data']] #data always  second
    fractions = ['','']
    regions = ['',''] #    regions = ['cr','cr']
    legList = [["bkg (mixed data 5nn)"], ["data"]]
    colorList = [430, 1]
    dofill = [True,False]
    isMC = False
    sf = [[0.0625],[1.]] #1.
    oname = 'comp_antiTagdata5nn_afterBDT'
    headerOpt = "   5nn, 4thjetCMVA<=medWP" #h1-h2 mass cut

elif which == 12:
    samples = [['antiMixed'],['bkg']] #data always  second
    fractions = ['','test']
    regions = ['cr','cr'] #    regions = ['cr','cr']
    legList = [["mixed data - antiTag"], ["mixed data - 4CMVA"]]
    colorList = [416, 430]
    sf = [[1.],[1.]]
    dofill = [True,True]
    isMC = False
    oname = 'comp_antiTagMixed_afterBDT'
    headerOpt = "    " #h1-h2 mass cut

elif which == 13:
    samples = [[''],['data']] #data always  second  -- amixed-0-5-10-15
    fractions = ['','']
    regions = ['sr','sr'] # check on 0.8-1 region
    legList = [["bkg (mixed data)"], ["data"]] # -- nn 0-5-10-15
    colorList = [430, 1]
    dofill = [True,False]
    sf = [[0.25],[1.]]
    isMC = False
    oname = 'comp_antiTagdata_nn-0-5-10-15'
    headerOpt = "  BDT[0.8-1] & 4thjetCMVA<=medWP" #h1-h2 mass cut

elif which == 14: # -l 0 !!
    getVar = True
  #  samples = [['amixed-0','amixed-5','amixed-12','amixed-18','amixed-24'],[]] # bin variance
  #  samples = [['amixed-0','amixed-1','amixed-2','amixed-3','amixed-4','amixed-5','amixed-6','amixed-7','amixed-8','amixed-9','amixed-10','amixed-11','amixed-12','amixed-13','amixed-14','amixed-15','amixed-16','amixed-17','amixed-18','amixed-19','amixed-20','amixed-21','amixed-22','amixed-23','amixed-24'],[]] # bin variance
    samples = [['amixed-91','amixed-92','amixed-93','amixed-94','amixed-95','amixed-96','amixed-97','amixed-98','amixed-99'],[]]
#    samples = [['amixed-0','amixed-11','amixed-22','amixed-33','amixed-44','amixed-55','amixed-66','amixed-77','amixed-88','amixed-99'],[]]
#    samples = [['amixed-0','amixed-11','amixed-22','amixed-33','amixed-44','amixed-55','amixed-66','amixed-77','amixed-88','amixed-99'],[]]
    fractions = ['','']
    regions = ['','']
    legList = [["bkg (mixed data - 10th nn)"], []] #up to 5 nn
    colorList = [430, 1]
    dofill = [False,False]
    isMC = False
    oname = 'comp_antiTagdata_nnVar_10'
    headerOpt = "  4thjetCMVA<=medWP"

## 10nn
elif which == 15: # -l 0 !!
    getChi = True
    xbmin = 21 #41 21
    #labels = ['1-1','2-2','3-3','4-4','5-5','6-6','7-7','8-8','9-9','10-10']
    #samples = [['mixed-0','mixed-11','mixed-22','mixed-33','mixed-44','mixed-55','mixed-66','mixed-77','mixed-88','mixed-99'],['data']]
    #samples = [['mixed-0','mixed-1','mixed-2','mixed-3','mixed-4','mixed-5','mixed-6','mixed-7','mixed-8','mixed-9'],['data']]
    labels = ['9-10','10-9','8-7','7-8','6-5','5-6','4-3','3-4','2-1','1-2']
    samples = [['mixed-89','mixed-98','mixed-76','mixed-67','mixed-54','mixed-45','mixed-32','mixed-23','mixed-10','mixed-1'],['data']]
    fractions = ['','']
    regions = ['','']
    legList = [["bkg (mixed data - 10nn)"], ["data"]]
    colorList = [430, 1]
    dofill = [False,False]
    isMC = False
    #oname = 'comp_mixedData_10nn_chi_diag_30b'
    oname = 'comp_mixedData_10nn_chi_rdm_30b'
    #headerOpt = "  antitag selection - 10nn diagonal - 30bins"
    headerOpt = "  antitag selection - 10nn various comb. - 30bins"

## 10nn
elif which == 16: # -l 0 !!
    getVar = True
    samples = [['mixed-0','mixed-11','mixed-22','mixed-33','mixed-44','mixed-55','mixed-66','mixed-77','mixed-88','mixed-99'],[]]
    #samples = [['mixed-0','mixed-1','mixed-2','mixed-3','mixed-4','mixed-5','mixed-6','mixed-7','mixed-8','mixed-9'],[]]
    #samples = [['mixed-89','mixed-98','mixed-76','mixed-67','mixed-54','mixed-45','mixed-32','mixed-23','mixed-10','mixed-1'],[]]
    fractions = ['','']
    regions = ['','']
    legList = [["bkg (mixed data - 10nn)"], []]
    colorList = [430, 1]
    dofill = [False,False]
    isMC = False
    oname = 'comp_mixedData_10nn_var_diag'
    #oname = 'comp_mixedData_10nn_var_rdm'
    headerOpt = "  antitag selection - 10nn diagonal"
    #headerOpt = "  antitag selection - 10nn various comb."

## 20nn - short
elif which == 17: # -l 0 !!
    getVar = True
    samples = [['mixed-0','mixed-1','mixed-2','mixed-3','mixed-4','mixed-5','mixed-6','mixed-7','mixed-8','mixed-9',
                'mixed-10','mixed-11','mixed-12','mixed-13','mixed-14','mixed-15','mixed-16','mixed-17','mixed-18 mixed-38'],[]]
#    samples = [['mixed-19','mixed-20','mixed-21','mixed-22','mixed-23','mixed-24','mixed-25','mixed-26','mixed-27',
#                'mixed-28','mixed-29','mixed-30','mixed-31','mixed-32','mixed-33','mixed-34','mixed-35','mixed-36','mixed-37','mixed-38'],[]]
    fractions = ['','']
    regions = ['','']
    legList = [["bkg (mixed data - 20nn)"], []] #up to 5 nn
    colorList = [430, 1]
    dofill = [False,False]
    isMC = False
    oname = 'comp_mixedData_20nn_var_diag'
#    oname = 'comp_mixedData_20nn_var_20'
    headerOpt = "  def selection - 20nn diagonal"
#    headerOpt = "  def selection - 20nn comb. with 20th"

## 20nn - short
elif which == 18: # -l 0 !!
    getChi = True
    xbmin = 41 #21
#    labels = ['11','22','33','44','55','66','77','88','99','1010','1111','1212','1313','1414','1515','1616','1717','1818','1919','2020']
#    samples = [['mixed-0','mixed-1','mixed-2','mixed-3','mixed-4','mixed-5','mixed-6','mixed-7','mixed-8','mixed-9',
#                'mixed-10','mixed-11','mixed-12','mixed-13','mixed-14','mixed-15','mixed-16','mixed-17','mixed-18','mixed-38'],['data']]

    labels = ['201','202','203','204','205','206','207','208','209','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020']
    samples = [['mixed-19','mixed-20','mixed-21','mixed-22','mixed-23','mixed-24','mixed-25','mixed-26','mixed-27',
                'mixed-28','mixed-29','mixed-30','mixed-31','mixed-32','mixed-33','mixed-34','mixed-35','mixed-36','mixed-37','mixed-38'],['data']]
    fractions = ['','']
    regions = ['','']
    legList = [["bkg (mixed data - 20nn)"], ["data"]]
    colorList = [430, 1]
    dofill = [False,False]
    isMC = False
#    oname = 'comp_antiTagdata_20nn_Chi_diag'
    oname = 'comp_antiTagdata_20nn_Chi_20'
#    headerOpt = "  antitag selection - 20nn diagonal - 10bins"
    headerOpt = "  antitag selection - 20nn comb. with 20th - 10bins"

elif which == 20: # -l 0 !!
    samples = [['ggHbb'],['sig']] #data always  second
    fractions = ['','appl']
    regions = ['',''] #    regions = ['cr','cr']
    legList = [["ggHbb"], ["HH4b SM"]]
    colorList = [[602], [632]]
    dofill = [True,True]
    isMC = True
    sf = [[0.002879],[(33.53*0.5824*0.5824/(4172119.0*0.2))]] 
    oname = 'comp_gghsig_afterBDT'
    headerOpt = "    "

elif which == 21: # -l 0 !!
    samples = [['vbfHbb'],['sig']] #data always  second
    fractions = ['','appl']
    regions = ['',''] #    regions = ['cr','cr']
    legList = [["vbfHbb"], ["HH4b SM"]]
    colorList = [[434], [632]]
    dofill = [True,True]
    isMC = True
    sf = [[0.000720],[(33.53*0.5824*0.5824/(4172119.0*0.2))]]
    oname = 'comp_vbfhsig_afterBDT'
    headerOpt = "    "

elif which == 22: # -l 0 !!
    samples = [['ttHbb'],['sig']] #data always  second
    fractions = ['','appl']
    regions = ['',''] #    regions = ['cr','cr']
    legList = [["ttHbb"], ["HH4b SM"]]
    colorList = [[419], [632]]
    dofill = [True,True]
    isMC = True
    sf = [[0.00007621],[(33.53*0.5824*0.5824/(4172119.0*0.2))]]

    oname = 'comp_tthsig_afterBDT'
    headerOpt = "    "

elif which == 23: # -l 0 !!
    samples = [['ZHbbqq'],['sig']] #data always  second
    fractions = ['','appl']
    regions = ['',''] #    regions = ['cr','cr']
    legList = [["ZH to QQBB"], ["HH4b SM"]]
    colorList = [[398], [632]]
    dofill = [True,True]
    isMC = True
    sf = [[0.0007665],[(33.53*0.5824*0.5824/(4172119.0*0.2))]]
    oname = 'comp_zhsig_afterBDT'
    headerOpt = "    "

elif which == 24: # -l 0 !!
    samples = [['TTTT'],['sig']] #data always  second
    fractions = ['','appl']
    regions = ['',''] #    regions = ['cr','cr']
    legList = [["TTTT"], ["HH4b SM"]]
    colorList = [[413], [632]]
    dofill = [True,True]
    isMC = True
    sf = [[1.],[(33.53*0.5824*0.5824/(4172119.0*0.2))]] ## update xs!!
    oname = 'comp_ttttsig_afterBDT'
    headerOpt = "    "

elif which == 25: # -l 0 !!
    samples = [['ttbb'],['sig']] #data always  second
    fractions = ['','appl']
    regions = ['',''] #    regions = ['cr','cr']
    legList = [["ttbb"], ["HH4b SM"]]
    colorList = [[420], [632]]
    dofill = [True,True]
    isMC = True
    sf = [[1.],[(33.53*0.5824*0.5824/(4172119.0*0.2))]] ## update xs!!
    oname = 'comp_ttbbsig_afterBDT'
    headerOpt = "    "

elif which == 26: # -l 0 !!
    samples = [['TT'],['sig']] #data always  second
    fractions = ['','appl']
    regions = ['','']
    legList = [["tt"], ["HH4b SM"]]
    colorList = [[430], [632]]
    dofill = [True,True]
    isMC = True
    sf = [[0.01077],[(33.53*0.5824*0.5824/(4172119.0*0.2))]] ## update xs!!
    oname = 'comp_ttsig_afterBDT'
    headerOpt = "    "

elif which == 27: # -l 0 !!
    samples = [['TT'],['bkg']] #data always  second
    fractions = ['','appl']
    regions = ['',''] #    regions = ['cr','cr']
    legList = [["tt"], ["mixed data"]]
    colorList = [[425], [400]]
    dofill = [True,True]
    isMC = True
    sf = [[0.01077*35.9*0.96],[1.]] ## update xs!!
    oname = 'comp_ttbkg_afterBDT'
    headerOpt = "    "

elif which == 207: # -l 0 !!
    samples = [['TT-fix-00'],['TT']] #data always  second
    fractions = ['','']
    regions = ['','']
    legList = [["mixed data as TT (0-0)"], ["TT"]]
    colorList = [[603], [430]]
    dofill = [True,True]
    isMC = True
    sf = [[0.01077*35.9*0.96],[0.01077*35.9*0.96]]
    oname = 'tt-tt-fix-00'
    headerOpt = "   "

elif which == 208: # -l 0 !!
    samples = [['ttHbb-fix-00'],['ttHbb']] #data always  second
    fractions = ['','']
    regions = ['','']
    legList = [["mixed data as ttHbb (0-0)"], ["ttHbb"]]
    colorList = [[619], [419]]
    dofill = [True,True]
    isMC = True
    sf = [[0.00007621*35.9*0.96],[0.00007621*35.9*0.96]]
    oname = 'tth-tth-fix-00'
    headerOpt = "   "

elif which == 28: # -l 0 !!
    samples = [['TT-fix-large-appl'],['TT-fix-00-11']] #data always  second
    fractions = ['','']
    regions = ['','']
    legList = [["re-mixed data as tt -large-appl"], ["re-mixed data as tt - 00-11"]]
    colorList = [[603], [430]]
    dofill = [True,True]
    isMC = True
    sf = [[0.01077*35.9*0.96/96],[0.01077*35.9*0.96*4]]
    oname = 'comp_tt-fix-largeappl-tt-fix-0011_afterBDT'
    headerOpt = "   "

elif which == 29: # -l 0 !!
    samples = [['ttHbb-fix-00-11'],['ttHbb-fix-00']] #data always  second
    fractions = ['','']
    regions = ['','']
    legList = [["re-mixed data as ttHbb"], ["mixed data as ttHbb"]]
    colorList = [[619], [419]]
    dofill = [True,True]
    isMC = True
    sf = [[0.00007621],[0.00007621]]
    oname = 'comp_ttHbkgttH1100_afterBDT'
    headerOpt = "    "

elif which == 30: # -l 0 !!
    samples = [['ttbb-fix-00-11'],['ttbb']] #data always  second
    fractions = ['','']
    regions = ['','']
    legList = [["re-mixed data as ttbb"], ["ttbb"]]
    colorList = [[613], [420]]
    dofill = [True,True]
    isMC = True
    sf = [[1.],[1.]] ## update xs!!
    oname = 'comp_ttbbBkgttbb11_afterBDT'
    headerOpt = "    "

elif which == 31: # -l 0 !!
    samples = [['TTTT-fix-00-11'],['TTTT']] #data always  second
    fractions = ['','']
    regions = ['','']
    legList = [["re-mixed data as tttt"], ["tttt"]]
    colorList = [[613], [413]]
    dofill = [True,True]
    isMC = True
    sf = [[1.],[1.]] ## update xs!!
    oname = 'comp_ttttBkgtttt11_afterBDT'
    headerOpt = "    "

elif which == 32: # -l 0 !!
    samples = [['ggH-fix-00-11'],['ggH-fix-00']] #data always  second
    fractions = ['','']
    regions = ['','']
    legList = [["re-mixed data as ggHbb"], ["mixed data as ggHbb"]]
    colorList = [[613], [602]]
    dofill = [True,True]
    isMC = True
    sf = [[0.002879],[0.002879]]
    oname = 'comp_ggHbkgggH1100_afterBDT'
    headerOpt = "    "

elif which == 33: # -l 0 !!
    samples = [['vbfH-fix-00-11'],['vbfH-fix-00']] #data always  second
    fractions = ['','']
    regions = ['','']
    legList = [["re-mixed data as vbfHbb"], ["mixed data as vbfHbb"]]
    colorList = [[613], [434]]
    dofill = [True,True]
    isMC = True
    sf = [[0.000720],[0.000720]]
    oname = 'comp_vbfHbkgvbfH1100_afterBDT'
    headerOpt = "    "

elif which == 34: # -l 0 !!
    samples = [['ZH-fix-00-11'],['ZH-fix-00']] #data always  second
    fractions = ['','']
    regions = ['','']
    legList = [["re-mixed data as ZHbbqq"], ["mixed data as ZHbbqq"]]
    colorList = [[613], [398]]
    dofill = [True,True]
    isMC = True
    sf = [[0.0007665],[0.0007665]]
    oname = 'comp_ZHbkgZH1100_afterBDT'
    headerOpt = "    "

elif which == 35: # -l 0 !!
    samples = [['QCD-fix-fix'],['QCD-fix']] #data always  second
    fractions = ['','']
    regions = ['','']
    legList = [["re-mixed data as QCD"], ["mixed data as QCD"]]
    colorList = [[613,613,613,613,613,613,613], [603,603,603,603,603,603,603]]
    dofill = [True,False]
    isMC = True
    sf = [[17.635231, 3.476259, 0.509821, 0.089584, 0.046491, 0.005935, 0.002410,],[17.635231, 3.476259, 0.509821, 0.089584, 0.046491, 0.005935, 0.002410,]]
    oname = 'comp_qcdBkgqcd1100_afterBDT'
    headerOpt = "    "

elif which == 300: # -l 0 !!
    samples = [['QCD-fix','TT-fix-00','ttHbb-fix-00','ZH-fix-00','vbfH-fix-00','ggH-fix-00'],['bkg']] #data always  second 
    fractions = ['','appl']
    regions = ['','']
    legList = [["mixed data as QCD","mixed data as QCD","mixed data as QCD","mixed data as QCD","mixed data as QCD","mixed data as QCD","mixed data as QCD","mixed data as tt","mixed data as ttH","mixed data as ZH","mixed data as vbfHbb","mixed data as ggHbb"], ["mixed data"]] #debug!!! 
    colorList = [[425,425,425,425,425,425,425,634,419,398,434,613], [603]]
    dofill = [True,False]
    isMC = True
    sf = [[17.635231*35.9, 3.476259*35.9, 0.509821*35.9, 0.089584*35.9, 0.046491*35.9, 0.005935*35.9, 0.002410*35.9, 0.01077*35.9, 0.00007621*35.9, 0.0007665*35.9, 0.000720*35.9, 0.002879*35.9],[0.25]]
    oname = 'comp_allBkgBkg_afterBDT'
    headerOpt = "    appl"

elif which == 301: # -l 0 !!
    samples = [['vbfH-fix-00-11','ggH-fix-00-11','ZH-fix-00-11','ttHbb-fix-00-11','TT-fix-00-11',],['bkg']] #data always  second
    fractions = ['','appl']
    regions = ['','']
    legList = [["re-mixed data as vbfHbb","re-mixed data as ggHbb","re-mixed data as ZH","re-mixed data as ttH","re-mixed data as tt"], ["mixed data"]] #mixed data
    colorList = [[613,434,398,419,634], [603]]
    dofill = [True,False]
    isMC = True
    sf = [[0.002879*35.9*0.96, 0.000720*35.9*0.96, 0.0007665*35.9*0.96, 0.00007621*35.9*0.96, 0.01077*35.9*0.96],[0.25]]
    oname = 'comp_minBkgBkg11_afterBDT'
    headerOpt = "    appl"

elif which == 302: # -l 0 !!
    samples = [['vbfH-fix-00-11','ggH-fix-00-11','ZH-fix-00-11','ttHbb-fix-00-11','TT-fix-00-11'],['vbfHbb','ggHbb','ZHbbqq','ttHbb','TT']] #data always  second
    fractions = ['','']
    regions = ['','']
    legList = [["re-mixed data as vbfHbb","re-mixed data as ggHbb","re-mixed data as ZH","re-mixed data as ttH","re-mixed data as tt"], ["MC samples"]] #mixed data
    colorList = [[613,434,398,419,634], [603, 603,603,603,603,603]]
    dofill = [True,False]
    isMC = True
    sf = [[0.002879*35.9*0.96*4, 0.000720*35.9*0.96*4, 0.0007665*35.9*0.96*4, 0.00007621*35.9*0.96*4, 0.01077*35.9*0.96*4],[0.002879*35.9*0.96*4, 0.000720*35.9*0.96*4, 0.0007665*35.9*0.96*4, 0.00007621*35.9*0.96*4, 0.01077*35.9*0.96*4]]
    oname = 'comp_minBkgminBkg11_afterBDT'
    headerOpt = "    "  
else: 
    print "ERROR: wrong '-w' argument"
    exit()
###############

if args.defaultCol: colors = [0,0]
else: colors = colorList

snames1 = []
for s in samples[0]:
    if not s in samlists: 
        if not s in samples: 
            snames1.append(s)
        else:
            snames1.append(samples[s]['sam_name'])    
    else: 
        snames1.extend(samlists[s])
#print snames1

snames2 = []
for s in samples[1]:
    if not s in samlists: 
        if not s in samples: 
            snames2.append(s)
        else:
            snames2.append(samples[s]['sam_name'])    
    else: 
        snames2.extend(samlists[s])
#print snames2

plotDirs1 = []
for sam in snames1:
    option = ''
    if fractions[0]: 
        option += fractions[0]
        if regions[0]: option += "_"
    if regions[0]: option += regions[0]

    if option: plotDirs1.append(sam+'_'+option)
    else: plotDirs1.append(sam)
print "HISTS FROM FOLDER {}".format(plotDirs1) 

plotDirs2 = []
for sam in snames2:
    option = ''
    if fractions[1]: 
        option += fractions[1]
        if regions[1]: option += "_"
    if regions[1]: option += regions[1]

    if option: plotDirs2.append(sam+'_'+option)
    else: plotDirs2.append(sam)
print "HISTS FROM FOLDER {}".format(plotDirs2) 


oDir = args.oDir
oDir += "/"+args.bdt
if not os.path.exists(oDir): os.mkdir(oDir)
oDir += "/"+oname
if args.doNorm: oDir = oDir+"_norm/"
else: oDir = oDir+"/"
if not os.path.exists(oDir): os.mkdir(oDir)
oDir += option #keep the second sample options
if not os.path.exists(oDir): os.mkdir(oDir)

#ran1v = (1.5, 1.5, 1.18, 1.18, 1.32, 1.25, 1.25, 1.22, 1.3, 1.15, 1.4, 1.4, 1.06, 1.3, 1.11, 1.08, 1.5, 1.5, 1.2, 1.15, 1.25, 1.35, 1.06)
#ran2v = (0.88, 0.88, 0.94, 0.93, 0.93, 0.93, 0.9, 0.9, 0.93, 0.9, 0.65, 0.88, 0.94, 0.88, 0.92, 0.92, 0.87, 0.93, 0.89, 0.95, 0.78, 0.9, 0.94)

#----------------------------------
for n, h in enumerate(histList):
 #   ran1 = ran1v[n]
  #  ran2 = ran2v[n]
    hOpt = hist_opt[h]
    if h == 'classifier': 
        h+='-'+args.bdt
    hs1 = UtilsDraw.getHistos_bdt(h, filename, plotDirs1, weights[0], sf[0])
    hs2 = UtilsDraw.getHistos_bdt(h, filename, plotDirs2, weights[1], sf[1])

    if drawH2:
        UtilsDraw.drawH2(hs1, hs2, hist_opt["h2_bdt"], snames1, args.clrebin, oDir, legList)
    elif getVar: # variance check
        UtilsDraw.drawBinVar(hs1, snames1, legList[0], hOpt, oDir, args.clrebin, headerOpt, isMC)
    elif getChi: # chi square
        UtilsDraw.drawChiSquare(hs1, snames1, legList[0], hs2, hOpt, oDir, xbmin, headerOpt, isMC, labels)
    else: 
        if hs1 and hs2:
            n1,n1err,n2,n2err = UtilsDraw.drawH1(hs1, snames1, legList[0], hs2, snames2, legList[1], 
                         hOpt, args.plotResidual, args.doNorm, oDir, colors, dofill, args.clrebin, headerOpt, isMC)#, ran1,ran2)
        #if n2: 
        #   print "### n1/n2 numEvents: {} +- {} ###".format(n1/n2, UtilsDraw.getRelErr(n1,n1err,n2,n2err)*n1/n2) 
        #   print "### n1: {} +- {} ###".format(n1,n1err) 
        #   print "### n2: {} +- {} ### \n".format(n2,n2err) 

