import json
import os
import math
from glob import glob

import numpy as np
# ROOT imports
from ROOT import TChain, TPad, TH1D, TH2D, TFile, vector, TCanvas, TLatex, TLine, TLegend, THStack, gStyle, TGaxis, TPolyLine, TGraph
import ROOT
from rootpy.plotting import Hist, Hist2D
from array import array
from Analysis.alp_analysis.alpSamplesOptions  import sam_opt


# UTILS
#################
def getNames(samlist):
    snames = []
    for s in samlist:
        snames.extend(samlists[s]) 
    return snames
#------------

def getFiles(samlist):
    files = []
    samNames = getNames(samlist)
    return files
#------------

def getRelErr(a,erra,b,errb):
    return math.sqrt(math.pow(erra/a,2)+math.pow(errb/b,2))
#------------

def getHistos_bdt(hist, filename, plotDirs, lumi, normtolumi, weights, sf, bias_unc_only = False):
    hlist = [] 
    
    tf = TFile(filename)
    #tf = TFile("../hh2bbbb_limit//20171120-160644-bm0.root")
    if not tf: 
        print "## WARNING: files do not exist"  
    for i, Dir in enumerate(plotDirs):
        hname = hist+"_"+Dir 
        if(tf.Get(Dir+"/"+hname)):
            norm = 1
            w = 1.
            sf_ = 1.
            if len(weights)>0:
	            if weights[i]>0:
        	        w = weights[i]
            if len(sf)>0:
                if sf[i]>=0:
                    sf_ = sf[i]

            if normtolumi: norm = w*lumi*sf_
            else: norm = w*sf_
            print "scale[{}]: {}".format(i, norm)
            h = tf.Get(Dir+"/"+hname)
            h.Scale(norm)
            hlist.append(h)
        else:
            print "## WARNING: hist {} not found in {}".format(Dir+"/"+hname,tf)
        print hname, hlist
    print hlist    
    return hlist
#------------

def getHistos_tdr(hist, filelist, plotDir, lumi, normtolumi, weight):
    hlist = []    
    for i, f in enumerate(filelist):
        w = 1.
        tf = TFile(f)
        if not tf: 
            print "## WARNING: files do not exist"  

        if weight[i]>=0:
   	    w = weight[i]

        if(tf.Get(plotDir+"/"+hist)):
            h = tf.Get(plotDir+"/"+hist)
            h.Scale(w)
            hlist.append(h)
            #print h.Integral()
        else:
            print "## WARNING: hist {} not found in {}".format(hist,tf)

    return hlist
#------------


def getHistos(hist, filelist, plotDir, lumi, normtolumi, weight, sf):
    hlist = []    
    for i, f in enumerate(filelist):  #debug - not efficient to loop on file
        norm = 1
        w = 1.
        sf_ = 1.
        tf = TFile(f)
        if not tf: 
            print "## WARNING: files do not exist"  

        if "Run" in f:   #data
            if not "mixed" in f:
                print "is Data"
        else:
            if len(weight)>0:
                if weight[i]>=0:
                    w = weight[i]
            # h_w_oneInvFb plot is not maintained - commented out for the time being
            #else:
                #if not "mixed" in f: #only for MC
                    #if(tf.Get("h_w_oneInvFb")):
                    #    h = tf.Get("h_w_oneInvFb")
                    #    w = h.GetBinContent(1)
                    #else:
                    #    print "WARNING: 'h_w_oneInvFb' not found in {}".format(tf)
            if len(sf)>0:
                if sf[i]>=0:
                    sf_ = sf[i]
            if normtolumi: norm = w*lumi*sf_
            else: norm = w*sf_
        print "scale[{}]: {}".format(i, norm)
        if(tf.Get(plotDir+"/"+hist)):
            h = tf.Get(plotDir+"/"+hist)
            h.Scale(norm)
            hlist.append(h)
        else:
            print "## WARNING: hist {} not found in {}".format(hist,tf)

    return hlist
#------------

def setYmaxStack(h1, h2, herr, c, hname):
    if (h1.GetMaximum() > h2.GetMaximum()): ymax = (h1.GetMaximum()*c) 
    else: ymax = h2.GetMaximum()*c 
    h1.SetMaximum(ymax)
    h2.SetMaximum(ymax)
    herr.SetMaximum(ymax)
    if(hname == 'h_nevts'):
        h1.SetMinimum(0.)
        h2.SetMinimum(0.)
        herr.SetMinimum(0.)
#------------

def setYmaxStack_D(hD, h1, h2, herr, c, hname):
    if hD.GetMaximum() > h1.GetMaximum() and hD.GetMaximum() > h2.GetMaximum(): ymax = hD.GetMaximum()*c
    elif h1.GetMaximum() > h2.GetMaximum(): ymax = h1.GetMaximum()*c 
    else: ymax = h2.GetMaximum()*c 

    hD.SetMaximum(ymax)
    h1.SetMaximum(ymax)
    h2.SetMaximum(ymax)
    herr.SetMaximum(ymax)
    if(hname == 'h_nevts'):
        hD.SetMinimum(0.)
        h1.SetMinimum(0.)
        h2.SetMinimum(0.)
        herr.SetMinimum(0.)
#------------

def setYmax(hs, c, hname):
    ymax = 0.
    for h in hs:
        if (h.GetMaximum() > ymax): ymax = (h.GetMaximum()*c) 
    for h in hs:
        h.SetMaximum(ymax)
#------------

def getScale(hsOpt, hlist1, hlist2, skip1=-1, skip2=-1): #h2/h1

    for h in hlist1:         
        h.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
    for h in hlist2:         
        h.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])

    hsInt = 0
    hsSkip = 0
    norm = 0
    for i, h in enumerate(hlist1):
        if i==skip1: 
            hsSkip += h.Integral()
            print 'skip', hsSkip
            continue 
        hsInt += h.Integral()
        #print h.Integral()
    for i, h in enumerate(hlist2): 
        if isinstance(h,int):
           norm = h
           break
        elif i==skip2: 
            continue 
        else: norm += h.Integral()
    if  hsInt: scale = norm/(hsInt+hsSkip)
    else: scale = 1.
    #print hsInt
    return scale
#------------


def getStackH(histos, hsOpt, rebin, snames, color, scale, fill, postfit_file = None, residuals = None):

    if color: col = color[0]
    else: col = sam_opt[snames[0]]['fillcolor']

    if postfit_file:
        myfile = TFile.Open(postfit_file)
        fit = "postfit"
        bak = myfile.Get("shapes_fit_s").Get("hh_bbbb").Get("total_background")
        sig = myfile.Get("shapes_fit_s").Get("hh_bbbb").Get("total_signal")
        data = myfile.Get("shapes_fit_s").Get("hh_bbbb").Get("data")
        
        if snames[0] == "data":
            for i in range(0,66):
                x, y = array('d', [0]), array('d', [0])
                data.GetPoint(i, x, y)
            histos[0].SetBinContent(0, 0)
            histos[0].SetBinError(0, 0)
        for ibin in range(1, sig.GetNbinsX()+1):
            #To avoid histos going out of scope later            
            for i in range(len(histos)):
                if snames[i] == "bkg":
                    histos[i].SetBinContent(ibin, bak.GetBinContent(ibin))
                    histos[i].SetBinError(ibin, bak.GetBinError(ibin))
                elif snames[i] == "sig":
                    histos[i].SetBinContent(ibin, sig.GetBinContent(ibin))
                    histos[i].SetBinError(ibin, sig.GetBinError(ibin))
                elif snames[i] == "data":
                    x, y = array('d', [0]), array('d', [0])
                    data.GetPoint(ibin-1, x, y)
                    histos[i].SetBinContent(ibin, y[0])
                    histos[i].SetBinError(ibin, data.GetErrorY(ibin))
        scale = 1
    
    
    
    herr = histos[0].Clone("hs_error")
    herr.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
    herr.Reset()
    herr.Rebin(rebin)
    herr.GetXaxis().SetTitle(hsOpt['xname'])
    herr.SetFillStyle(3005)
    herr.SetFillColor(col)
    herr.SetLineColor(922)
    herr.SetLineWidth(0)         
    herr.SetMarkerSize(0)
    herr.SetMarkerColor(922)
    herr.SetMinimum(0.)

    hs   = THStack("hs","")
    for i, h in enumerate(histos):
	    #print histos
	    #print i, color
        if color: col = color[i]
        else: col = sam_opt[snames[i]]['fillcolor']
        #print sam_opt[snames[i]]['sam_name']
        
        h.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        h.SetMinimum(0.)
        
        #Just a quick hack to answer ARC question
        #Actually not needed as it is done in classifier report already
        """if "bkg_appl_ms" in h.GetName():
            print i, h
            #bkg_bias_fname = "/lustre/cmswork/dcastrom/projects/hh/april_2017/CMSSW_8_0_25/src/Analysis/hh2bbbb_limit/notebooks/bms_btagside_err_fixed/BM0/bias_correction_bigset_unscaled.json"
            bkg_bias_fname = "/lustre/cmswork/dcastrom/projects/hh/april_2017/CMSSW_8_0_25/src/Analysis/hh2bbbb_limit/notebooks/bias_22032018_with_weights_also_mass_cut/BM0/bias_correction_mass_cut_bigset_unscaled.json"
            with open(bkg_bias_fname,"r") as bkg_bias_file:
                json_dict = json.load(bkg_bias_file)
                print ("using bias file: ", bkg_bias_file)

            hbias = h.Clone("h_bias")
            for n in range(len(json_dict['bias_corr'])):
                #if not s_bin.overflow:
                #value = s_bin.value

                bias = json_dict['bias_corr'][n]
                var = json_dict['var'][n]
                bias_unc = json_dict['bias_corr_unc_bs'][n]
                bias_unc_stat = json_dict['bias_corr_unc_stat'][n]
                
                new_bkg_pred_stat = var
                                
                new_bkg_pred_tot_unc = np.sqrt(new_bkg_pred_stat**2 + bias_unc**2 + bias_unc_stat**2)
                hbias.SetBinContent(n+1, bias)
                hbias.SetBinError(n+1, new_bkg_pred_tot_unc)
            h.Add(hbias, -1)
        """
        
        h.Scale(scale)
        h.Rebin(rebin)
        h.SetMarkerStyle(8)
        h.SetMarkerSize(0.)
        print "h[{}]: {}".format(i, h.Integral())
        if fill and residuals == 5: 
            h.SetFillStyle(0)
            #h.SetFillColor(col)
            h.SetLineColorAlpha(col, 0.6)
            h.SetLineStyle(1)
            h.SetLineWidth(2)            
        elif fill: 
            h.SetFillColorAlpha(col,0.35)
            h.SetFillStyle(1001)
        if i==len(histos)-1 :
            h.SetLineStyle(1)
            #h.SetLineWidth(2)
            h.SetLineColor(col)
            if not fill: 
                h.SetMarkerSize(0.6)
                h.SetMarkerColor(col)
        else:     
            h.SetLineWidth(1)         
            h.SetLineColorAlpha(col,0.)
            h.SetLineStyle(1)
            #h.SetLineWidth(2)
            h.SetLineColor(col)
        if i==0: h_  = h.Clone("h_")
        else: h_.Add(h)
        hs.Add(h)
        herr.Add(h)
        
    if "sig" in snames:
        i = snames.index("sig")
        sig100 = histos[i].Clone("sig x5")
        sig100.Scale(5)
        sig100.SetLineStyle(2)
        sig100.SetFillStyle(0)
        histos.append(sig100)
        
        
    return hs, herr, h_

# DRAWING FUNCTIONS
#################
def drawChiSquare(hlist, snames, legstack, h2, hsOpt, oDir, xbmin, headerOpt, isMC, labels):  
    gStyle.SetOptStat(False)
    c1 = TCanvas("c1", hsOpt['hname'], 800, 800)       

    x_data =[]
    chisq_ =[]
    n_sample = 0
    nbins = hlist[0].GetNbinsX()
    #print nbins

    for k in range(0, nbins-xbmin+1):
        x_data.append(0.)

    for i, h in enumerate(h2):
        for b in range(xbmin, nbins+1):
            x_data[b-xbmin] = h.GetBinContent(b)
            print b, x_data[b-xbmin] 
    
    for i, h in enumerate(hlist):            
        n_sample +=1
        chisq_.append(0)
        for b in range(xbmin, nbins+1):
            x = h.GetBinContent(b)       
            #print b, x
            #print b, x_data[b-xbmin] 
            if (x+x_data[b-xbmin]) > 0: chisq_[i] += pow((x_data[b-xbmin]-x),2)/(x+x_data[b-xbmin])
            else: 
                print 'sample data bin {} : null bin content?'.format(b)       

    print "\n", n_sample  

    hchi = (TH1D)("h_chisq", "chi square", n_sample, 0., n_sample) 
    hchi.SetMarkerStyle(8)
    hchi.SetMarkerSize(1.)
    hchi.SetMarkerColor(1)

    for i, h in enumerate(hlist):         
        hchi.SetBinContent(i+1, chisq_[i])
       #hchi.SetBinError(j+1, dDiff)
        #print hchi.GetBinContent(i+1) #, hvar.GetBinError(j+1)       

    x_axis = hchi.GetXaxis()
    for b_n, label in enumerate(labels): 
        x_axis.SetBinLabel(b_n+1, label)        

    hchi.GetYaxis().SetTitle("chi2 data-mixed")
    hchi.Draw("PE")
    drawCMS(35.9, headerOpt)

    c1.Update()    
    c1.SaveAs(oDir+"/"+"bdtChisq_allnn.pdf")
    c1.SaveAs(oDir+"/"+"bdtChisq_allnn.png")
    c1.SaveAs(oDir+"/"+"bdtChisq_allnn.root")
#------------

def drawBinVar(hlist, snames, legstack, hsOpt, oDir, rebin, headerOpt, isMC):  
    gStyle.SetOptStat(False)
    c1 = TCanvas("c1", hsOpt['hname'], 800, 800)       
    
    x_=[]
    xsq_=[]
    xs_mean = []
    xs_sqmean = []
    sum_err = []
    n_sample = 0
    hvar = hlist[0].Clone("hvar")
    hvar.Reset()

    #once for all
    nbins = hlist[0].GetNbinsX()
    #print nbins
    for k in range(0, nbins):
        x_.append(0.)
        xsq_.append(0.)
        sum_err.append(0.)

    for i, h in enumerate(hlist):         
        n_sample +=1
        for b in range(1, nbins+1):
            x = h.GetBinContent(b)
            #if b==nbins: print b, x
            x_[b-1] += x
            xsq_[b-1] += pow(x,2)
    #print "\n", n_sample

    for k in range(0, nbins):
        xs_mean.append(x_[k]/n_sample)
        xs_sqmean.append(xsq_[k]/n_sample)

    for i, h in enumerate(hlist):         
        for b in range(1, nbins+1):
            x = h.GetBinContent(b)
           #sum_err[b-1] +=  (pow(x,2) + pow(xs_mean[b-1],2) + 2*x*xs_mean[b-1])*x            
            sum_err[b-1] += pow((2*x - 2*xs_mean[b-1] -1),2)*x

    for j, m in enumerate(xs_mean):
        var = (xs_sqmean[j] - pow(m,2))
       #dVar = math.sqrt(4./pow(n_sample,2)*sum_err[j])
       #dVar = 2*pow(var,2)/(n_sample-1)
       #print var #, dVar
       #print m, math.sqrt(m)

       #if x>0. : 
       #r = var / m
       #dRat = math.sqrt(pow(dVar/var,2)+pow(math.sqrt(m)/m,2))*r
        diff = var - m        
        dDiff = math.sqrt(1./pow(n_sample,2)*sum_err[j])
       #else: 
        #    print 'null mean - r set to zero'
         #   diff = 0.
          #  dDiff = 0.

        hvar.SetBinContent(j+1, diff)
        hvar.SetBinError(j+1, dDiff)
        #print hvar.GetBinContent(j+1), hvar.GetBinError(j+1)       

    hvar.GetYaxis().SetTitle("mean-Variance")
    hvar.Draw("PE")
    legend = setLegend(1,1)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hvar,"mixed data")
    drawCMS(35.9, headerOpt)
    legend.Draw("same")
    l = TLine(0.,0.,1.,0.);
    l.SetLineStyle(3)
    l.Draw("same")
    c1.Update()    
    c1.SaveAs(oDir+"/"+"bdtVar_nn.pdf")
    c1.SaveAs(oDir+"/"+"bdtVar_nn.png")
    c1.SaveAs(oDir+"/"+"bdtVar_nn.root")
#------------

def drawH1(hlist, snames, legstack, hsOpt, residuals, norm, oDir, colors, dofill, rebin, headerOpt, isMC, fit_results = None, postfit_file = None, bm = 0):
    H_ref = 800     
    W_ref = 800
    if residuals == 5: H_ref = 1600
    W = W_ref
    H = H_ref
    
    iPos = 11
    iPeriod = 4
    
    c1 = TCanvas("c1", hsOpt['hname'], H_ref, W_ref)
    setTDRStyle()
    
    T = 0.08*H_ref
    B = 0.12*H_ref 
    L = 0.12*W_ref
    R = 0.04*W_ref
    
    c1.SetFillColor(0)
    c1.SetBorderMode(0)
    c1.SetFrameFillStyle(0)
    c1.SetFrameBorderMode(0)
    c1.SetLeftMargin( L/W )
    c1.SetRightMargin( R/W )
    c1.SetTopMargin( T/H )
    c1.SetBottomMargin( B/H )
    #c1.SetBottomMargin( 0 )
    
    if residuals == -1 or residuals == -2 or residuals == -4 or residuals == -12:
        pad1 = TPad("pad1", "pad1", 0, 0.4, 1, 1.0)
        pad1.SetTopMargin(0.1) 
        pad1.SetBottomMargin(0.03) 
        pad1.Draw()             
        pad1.cd()
        #if residuals == -2: pad1.SetLogy()        
    elif residuals == 5: 
        pad1 = TPad("pad1", "pad1", 0, 0., 0.5, 1.0)
        #pad1.SetTopMargin(0) 
        pad1.SetBottomMargin(0.09)
        #pad1.SetLeftMargin(0.)
        #pad1.SetRightMargin(0.) 
        pad1.Draw()             
        pad1.cd()
        #if residuals == -2: pad1.SetLogy() """       
    else:
        pad1 = c1
        
    if hsOpt["ylog"] == True:
        pad1.SetLogy()
    pad1.SetLogy()
    
    if rebin > 0: rb = rebin
    else: rb =  hsOpt['rebin']

    isNevts=False
    if hsOpt['hname']=="h_nevts": isNevts=True
    ymax = 0.00001
    
    scales = []
    for i in range(len(hlist) - 1):
        if(norm):
            scales.append(getScale(hsOpt, hlist[i], hlist[-1]))
            print "sc_to_norm%d: %d" % (i, scales[i])
        else: scales.append(1.)
    scales.append(1)
    hs, herr, h = [], [], []
    
    if residuals == 5:
        h_bkg_bias_corr = get_bias_corrected_histo(hlist[0][0], region = "ms")
        #print hlist, hlist[0][0].Integral()*4
        #for n in range(1, hlist[1][0].GetNbinsX()+1):
        #    print n, hlist[1][0].GetBinContent(n), hlist[0][0].GetBinContent(n) * 4, h_bkg_bias_corr.GetBinContent(n)*4, h_bkg_bias_corr.GetBinError(n)*2
    
    for i in range(len(hlist)):
        hs_tmp, herr_tmp, h_tmp = getStackH(hlist[i], hsOpt, rb, snames[i], colors[i], scales[i], dofill[i], postfit_file, residuals)
        if hsOpt["ylog"] == True:
            if hs_tmp.GetMaximum() > ymax: ymax = hs_tmp.GetMaximum()*10
            if herr_tmp.GetMaximum() > ymax: ymax = herr_tmp.GetMaximum()*10
        else:
            #if hs_tmp.GetMaximum() > ymax: ymax = hs_tmp.GetMaximum()*1.20
            #if herr_tmp.GetMaximum() > ymax: ymax = herr_tmp.GetMaximum()*1.20
            if hs_tmp.GetMaximum() > ymax: ymax = hs_tmp.GetMaximum()*40
            if herr_tmp.GetMaximum() > ymax: ymax = herr_tmp.GetMaximum()*40
        print "ymax", ymax
        #if isNevts:  print "h1Int ",  h1.GetBinContent(1), h1.GetBinError(1)
        hs.append(hs_tmp)
        herr.append(herr_tmp)
        h.append(h_tmp)

    #Only do KS if 2 sets of histos
    ks = None
    if len(hlist) == 2:
        if len(hlist[0]) == 1 and len(hlist[1]) == 1:
            ks = hlist[0][0].KolmogorovTest(hlist[1][0])
            print("KS: ", ks)
            print("Chi2: ", hlist[0][0].Chi2Test(hlist[1][0], "UU NORM"))
        else:
            kol1 = hlist[0][0].Clone()
            kol2 = hlist[1][0].Clone()
            maxn = len(hlist[0])
            if "sig" in snames[0]: maxn -= 1
            for ihist in range(1, maxn):
                #for i in range(1, kol1.GetNbinsX()+1):
                #    print ihist, i, hlist1[ihist].GetBinContent(i), kol1.GetBinContent(i)
                kol1.Add(hlist[0][ihist])
            negBin = False
            for i in range(1, kol1.GetNbinsX()+1):
                if kol1.GetBinContent(i) < 0:
                    negBin = True
            for ihist in range(1, len(hlist[1])):
                kol2.Add(hlist[1][ihist])
            ks = 0.        
            if not negBin:
                ks = kol1.KolmogorovTest(kol2, "NX")
            print("KS: ", ks)
        
    #debug -- needed before drawing hs
    herr[0].GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
    #herr[0].GetXaxis().SetNdivisions(000) 
    herr[0].SetMinimum(1)
    if residuals == 5:
        ymax = ymax / 2
    herr[0].SetMaximum(ymax)
    herr[0].GetYaxis().SetTitleSize(20)
    herr[0].GetYaxis().SetTitleFont(43)
    herr[0].GetYaxis().SetTitleOffset(1.40)
    herr[0].GetYaxis().SetLabelFont(43)
    herr[0].GetYaxis().SetLabelSize(18)
    herr[0].GetYaxis().SetTitle("Events")
    pad1.SetTickx(1)
    if len(hlist[0])>1: 
        herr[0].Draw("E2")
        herr[0].SetFillColor(922)
    else: herr[0].Draw("E")
    #--
    #Legend & headers
    #if isMC: drawCMS(-1, headerOpt)
    #else: drawCMS(35.9, headerOpt)
    CMS_lumi(pad1, iPeriod, iPos)
    
    if residuals == 5:
        legend = setLegend(0.4,0.7,0.90,0.9)
    else:
        legend = setLegend(0.7,0.60,0.90,0.85)    
    for i in range(1, len(hlist)):
        legend.AddEntry(hlist[i][len(hlist[i])-1], legstack[i][0]) #debug - one leg for second samplelist
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextAlign(33)   
    #if ks and residuals == -4:
    #    latex.DrawLatex(0.5, 0.78, "KS p-val: %.3f" % ks)
    nskip = 0
    match = False 
    for n, sam in enumerate(snames[0]):
        #print n, same
        m = len(snames[0]) - n - 1
        # to get one legend for all HT bins - to be implemented for other samples         
        if sam.find("QCD")>=0 and match: 
            nskip+=1 
            continue
        if len(legstack[0]) > n-nskip:
            if sam.find("QCD")>=0: match = True
            legend.AddEntry(hlist[0][n], legstack[0][n-nskip])            
    
    #-------------
    
    if residuals == 5:
        h_bkg_bias_corr.SetLineColor(ROOT.kBlue)
        h_bkg_bias_corr.SetFillStyle(1001)
        h_bkg_bias_corr.SetMarkerStyle(0)
        h_bkg_bias_corr.SetFillColorAlpha(ROOT.kBlue, 0.3)
        h_bkg_bias_corr.Draw("hist same")
        h_bkg_bias_corr_err = h_bkg_bias_corr.Clone("bkg_bias_corrected_unc")
        #herr.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        #herr.Reset()
        #herr.Rebin(rebin)
        #herr.GetXaxis().SetTitle(hsOpt['xname'])
        h_bkg_bias_corr_err.SetFillStyle(3005)
        h_bkg_bias_corr_err.SetFillColor(ROOT.kBlue)
        h_bkg_bias_corr_err.SetLineColor(922)
        h_bkg_bias_corr_err.SetLineWidth(0)         
        h_bkg_bias_corr_err.SetMarkerSize(0)
        h_bkg_bias_corr_err.SetMarkerColor(922)
        h_bkg_bias_corr_err.Draw("e2 same")
        legend.AddEntry(h_bkg_bias_corr, "Mixed data with bias correction")
        legend.AddEntry(h_bkg_bias_corr_err, "Background uncertainty")
        
    if len(hlist[0])>1: legend.AddEntry(herr[0], 'Total uncertainty')
    legend.Draw("same")

    if(ymax > 1000): TGaxis.SetMaxDigits(3)
    for i in range(len(hs)):
        hs[i].SetMaximum(ymax)
        herr[i].SetMaximum(ymax)
        plotH(hlist[i], hs[i], herr[i], dofill[i], residuals)
        if i == len(hs) - 1:
            herr[i].Draw("Esameaxis")
    
    
    #hlist1[-1].Draw("same hist")

    nevs = []
    neverrs = []
    if isNevts: 
        for i in range(len(herr)):
            nevs.append(herr[i].GetBinContent(1))
            neverrs.append(herr[i].GetBinError(1))
    
    """if len(hlist[0])>1: 
        herr[0].Draw("E2 same")
        herr[0].SetFillColor(922)
    else: herr[0].Draw("E same")"""

    if not (residuals == -1 or residuals == -2 or residuals == -4 or residuals == -12):
        c1.Update()    
        c1.SaveAs(oDir+"/"+hsOpt['hname']+".pdf")
        c1.SaveAs(oDir+"/"+hsOpt['hname']+".png")            
        #c1.SaveAs(oDir+"/"+hsOpt['hname']+".root") 
        herr[0].GetXaxis().SetLabelSize(0.025)
        herr[0].SetNdivisions(510, "X")
    else:
        herr[0].GetXaxis().SetLabelSize(0.)

    if residuals==-3: # division --- utility not clear...
        c1.cd()
        pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.4)
        pad2.SetTopMargin(0.)
        pad2.SetBottomMargin(0.2)
        pad2.Draw()
        pad2.cd()
        hrat = h[-1].Clone("h_div")
        hb = h[-1].Clone("h_bias")

        for i in range(0, len(bias)):
            hb.SetBinContent(i+1, bias[i])
            #print(hb.GetBinContent(i+1))

        for ibin in range(1, h[1].GetNbinsX()+1):
            hrat.SetBinContent(ibin, h[0].GetBinContent(ibin)-h[1].GetBinContent(ibin) )
            #print(hrat.GetBinContent(ibin))
            if(h2.GetBinContent(ibin)>0): 
                hrat.SetBinError(ibin, h[1].GetBinError(ibin)/h[1].GetBinContent(ibin) )
            else: 
                hrat.SetBinError(ibin, 0.)

        # MC uncertainy shadow plot
        h_error = herr[0].Clone("h_err")
        h_error.Reset()
        for ibin in range(1, herr[0].GetNbinsX()+1):
            if(herr[0].GetBinContent(ibin)>0): 
                h_error.SetBinContent (ibin, 1.)
                h_error.SetBinError   (ibin, herr[0].GetBinError(ibin)/herr[0].GetBinContent(ibin))
            else: 
                h_error.SetBinContent(ibin, 1.)
                h_error.SetBinError   (ibin, 0.)

        hrat.SetTitle("")
        hrat.GetXaxis().SetTitleSize(20)
        hrat.GetXaxis().SetTitleFont(43)
        hrat.GetXaxis().SetTitleOffset(4.)
        hrat.GetXaxis().SetLabelFont(43)
        hrat.GetXaxis().SetLabelSize(20)
        hrat.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        minbin = hrat.GetXaxis().GetFirst()
        maxbin = hrat.GetXaxis().GetLast()
        #print minbin, maxbin        
        ymax_ = 1000 #1.5 
        ymax = 800.
        ymin_ = 0.5
        ymin = -200.

       # print 'rb', rb      
        hrat.GetYaxis().SetRangeUser(ymin,ymax)
        hrat.GetYaxis().SetTitleSize(20)
        hrat.GetYaxis().SetTitleFont(43)
        hrat.GetYaxis().SetTitleOffset(1.40)
        hrat.GetYaxis().SetLabelFont(43)
        hrat.GetYaxis().SetLabelSize(18)
        hrat.SetMarkerStyle(8)
        hrat.SetMarkerSize(0.8)
        hrat.SetMarkerColor(1)
        hrat.SetLineColor(1)
        hrat.GetXaxis().SetTitle(hsOpt['xname'])
        hrat.GetYaxis().SetTitle('exp - obs')


        c00 = TCanvas("c0", "", 800, 800)    
        c00.cd()
        hrat.Draw("X0")
        hb.SetMarkerColor(2)
        hb.SetMarkerSize(0.8)
        hb.Draw("same X0")
        c00.SaveAs(oDir+"/"+hsOpt['hname']+"_div2.root")
        c00.SaveAs(oDir+"/"+hsOpt['hname']+"_div2.png")

        l = TLine(hsOpt['xmin'],1.5,hsOpt['xmax'],1.5);
        l0 = TLine(hsOpt['xmin'],1.4,hsOpt['xmax'],1.4);
        l00 = TLine(hsOpt['xmin'],1.3,hsOpt['xmax'],1.3);
        l000 = TLine(hsOpt['xmin'],1.2,hsOpt['xmax'],1.2);
        l1 = TLine(hsOpt['xmin'],1.1,hsOpt['xmax'],1.1);
        l2 = TLine(hsOpt['xmin'],1.,hsOpt['xmax'],1.);
        l3 = TLine(hsOpt['xmin'],0.9,hsOpt['xmax'],0.9);
        l4 = TLine(hsOpt['xmin'],0.8,hsOpt['xmax'],0.8);
        l.SetLineStyle(3)
        l0.SetLineStyle(3)
        l00.SetLineStyle(3)
        l000.SetLineStyle(3)
        l1.SetLineStyle(3)
        l2.SetLineStyle(3)
        l3.SetLineStyle(3)
        l4.SetLineStyle(3)
        l2.Draw("same")

        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_div.pdf")
        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_div.png")            
        #c1.SaveAs(oDir+"/"+hsOpt['hname']+"_div.root") 
        c1.Clear()

    if residuals==-1: # simple ratio -- data as h2!!!
        c1.cd()
        pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0.)
        pad2.SetBottomMargin(0.2)
        pad2.Draw()
        pad2.cd()
        hrat = h[1].Clone("h_rat")
        hrat.Divide(h[0])
        # consider only data error in the ratio plot
        for ibin in range(1, h[1].GetNbinsX()+1):
            #print(h2.GetBinContent(ibin))
            if(h[1].GetBinContent(ibin)>0): 
                hrat.SetBinError(ibin, h[1].GetBinError(ibin)/h[1].GetBinContent(ibin) )
            else: 
                hrat.SetBinError(ibin, 0.)

        # MC uncertainy shadow plot
        h_error = herr[0].Clone("h_err")
        h_error.Reset()
        for ibin in range(1, herr[0].GetNbinsX()+1):
            if(herr[0].GetBinContent(ibin)>0): 
                h_error.SetBinContent (ibin, 1.)
                h_error.SetBinError   (ibin, herr[0].GetBinError(ibin)/herr[0].GetBinContent(ibin))
            else: 
                h_error.SetBinContent(ibin, 1.)
                h_error.SetBinError   (ibin, 0.)

        hrat.SetTitle("")
        hrat.GetXaxis().SetTitleSize(20)
        hrat.GetXaxis().SetTitleFont(43)
        hrat.GetXaxis().SetTitleOffset(4.)
        hrat.GetXaxis().SetLabelFont(43)
        hrat.GetXaxis().SetLabelSize(20)
        hrat.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        minbin = hrat.GetXaxis().GetFirst()
        maxbin = hrat.GetXaxis().GetLast()
        ymax_ = 1.5
        ymax = 1.
        ymin_ = 0.5
        ymin = 1.
    
        for ibin in range(minbin, maxbin+1):       
            binc = hrat.GetBinContent(ibin) 
            if binc == 0. or (binc is None): 
                continue
            if binc*1.05 > ymax: 
                ymax = binc*1.05
            if binc*0.95 < ymin: 
                ymin = binc*0.95

        hrat.GetYaxis().SetRangeUser(ymin,ymax)
        hrat.GetYaxis().SetTitleSize(20)
        hrat.GetYaxis().SetTitleFont(43)
        hrat.GetYaxis().SetTitleOffset(1.40)
        hrat.GetYaxis().SetLabelFont(43)
        hrat.GetYaxis().SetLabelSize(18)
        hrat.SetMarkerStyle(8)
        hrat.SetMarkerSize(0.8)
        hrat.SetMarkerColor(1)
        hrat.SetLineColor(1)
        hrat.GetXaxis().SetTitle(hsOpt['xname'])
        hrat.GetYaxis().SetTitle('data/bkg')
        hrat.Draw("E X0")
        h_error.SetFillColor(430)
        h_error.Draw("E2same")

        hrat.Fit("pol1")
        myfunc = hrat.GetFunction("pol1")
        print myfunc
        fithist = myfunc.CreateHistogram()
        for ibin in range(1,fithist.GetNbinsX()+1):
            print "%.5f," % (fithist.GetBinContent(ibin)),
        print
        l = TLine(hsOpt['xmin'],1.5,hsOpt['xmax'],1.5);
        l0 = TLine(hsOpt['xmin'],1.4,hsOpt['xmax'],1.4);
        l00 = TLine(hsOpt['xmin'],1.3,hsOpt['xmax'],1.3);
        l000 = TLine(hsOpt['xmin'],1.2,hsOpt['xmax'],1.2);
        l1 = TLine(hsOpt['xmin'],1.1,hsOpt['xmax'],1.1);
        l2 = TLine(hsOpt['xmin'],1.,hsOpt['xmax'],1.);
        l3 = TLine(hsOpt['xmin'],0.9,hsOpt['xmax'],0.9);
        l4 = TLine(hsOpt['xmin'],0.8,hsOpt['xmax'],0.8);
        l.SetLineStyle(3)
        l0.SetLineStyle(3)
        l00.SetLineStyle(3)
        l000.SetLineStyle(3)
        l1.SetLineStyle(3)
        l2.SetLineStyle(3)
        l3.SetLineStyle(3)
        l4.SetLineStyle(3)
        #l.Draw("same")
        #l0.Draw("same")
        #l00.Draw("same")
        #l000.Draw("same")
        #l1.Draw("same")
        l2.Draw("same")
        #if ymin<0.9: l3.Draw("same")
        #if ymin<0.8: l4.Draw("same")

        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.pdf")
        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.png")            
        c1.Clear()


    elif residuals == -2: # plot data - fitted (sig+bkg) in residual plot
        c1.cd()
        pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.4)
        pad2.SetTopMargin(0.)
        pad2.SetBottomMargin(0.2)
        pad2.Draw()
        pad2.cd()
        
        #hlist1, snames1, legstack1, hlist2, snames2, legstack2, hsOpt, residuals, norm, oDir, colors, dofill, rebin, headerOpt, isMC
        histos = {}
        histos["data"] = hlist[1][0]
        histos["sig"] = hlist[0][1]
        histos["bkg"] = hlist[0][0]
        #print hlist
        #for binn in range(1, histos["data"].GetNbinsX()+1):
        #    print binn, histos["data"].GetBinContent(binn), histos["sig"].GetBinContent(binn), histos["bkg"].GetBinContent(binn)
        (h_data_bkg, h_sig, h_error) = getHistosPostFit(histos, hsOpt, snames[0], colors, fit_results, postfit_file)
        #print h_data_bkg.Integral(), h_sig.Integral(), h_err.Integral()
        
        h_data_bkg.SetTitle("")
        h_data_bkg.GetXaxis().SetTitleSize(20)
        h_data_bkg.GetXaxis().SetTitleFont(43)
        h_data_bkg.GetXaxis().SetTitleOffset(4.)
        h_data_bkg.GetXaxis().SetLabelFont(43)
        h_data_bkg.GetXaxis().SetLabelSize(20)
        h_data_bkg.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        #h_data_bkg.GetYaxis().SetRangeUser(-20,20)
        minbin = h_data_bkg.GetXaxis().GetFirst()
        maxbin = h_data_bkg.GetXaxis().GetLast()
        y_max = max(h_data_bkg.GetMaximum(), hlist[0][-1].GetMaximum(), h_error.GetMaximum())*1.25
        y_min = min(h_data_bkg.GetMinimum(), hlist[0][-1].GetMinimum(), h_error.GetMinimum())*1.25

        h_data_bkg.GetYaxis().SetRangeUser(y_min,y_max)
        h_data_bkg.GetYaxis().SetTitleSize(20)
        h_data_bkg.GetYaxis().SetTitleFont(43)
        h_data_bkg.GetYaxis().SetTitleOffset(1.40)
        h_data_bkg.GetYaxis().SetLabelFont(43)
        h_data_bkg.GetYaxis().SetLabelSize(18)
        #h_data_bkg.SetMarkerStyle(33)
        h_data_bkg.SetMarkerStyle(20)
        h_data_bkg.SetMarkerSize(0.7)
        h_data_bkg.SetMarkerColor(1)
        h_data_bkg.SetLineColor(1)
        h_data_bkg.GetXaxis().SetTitle(hsOpt['xname'])
        #h_data_bkg.GetYaxis().SetTitle('data/bkg')
        h_data_bkg.Draw("e3 x0")
        hlist[0][-1].Draw("hist same")
        
        h_sig.Draw("hist same")
        
        h_error.SetFillColor(632)
        h_error.Draw("E2same")
        h_data_bkg.Draw("same e3 x0")
        
        #Draw best fit slope from MS CR
        if True:
            slope = [0.98610, 0.98652, 0.98695, 0.98738, 0.98781, 0.98823, 0.98866, 0.98909, 0.98951, 0.98994, 0.99037, 0.99080, 0.99122, 0.99165, 0.99208, 0.99250, 0.99293, 0.99336, 0.99379, 0.99421, 0.99464, 0.99507, 0.99549, 0.99592, 0.99635, 0.99678, 0.99720, 0.99763, 0.99806, 0.99848, 0.99891, 0.99934, 0.99977, 1.00019, 1.00062, 1.00105, 1.00147, 1.00190, 1.00233, 1.00276, 1.00318, 1.00361, 1.00404, 1.00446, 1.00489, 1.00532, 1.00575, 1.00617, 1.00660, 1.00703, 1.00745, 1.00788, 1.00831, 1.00874, 1.00916, 1.00959, 1.01002, 1.01044, 1.01087, 1.01130, 1.01173, 1.01215, 1.01258, 1.01301, 1.01343, 1.01386, 1.01429, 1.01472, 1.01514, 1.01557, 1.01600, 1.01643, 1.01685, 1.01728, 1.01771, 1.01813, 1.01856, 1.01899, 1.01942, 1.01984, 1.02027, 1.02070, 1.02112, 1.02155, 1.02198, 1.02241, 1.02283, 1.02326, 1.02369, 1.02411, 1.02454, 1.02497, 1.02540, 1.02582, 1.02625, 1.02668, 1.02710, 1.02753, 1.02796, 1.02839,]
            print slope
            bkg_slope = histos["data"].Clone("bkg_slope")
            for b in range(1, bkg_slope.GetNbinsX()+1):
                bkg_slope.SetBinContent(b, bkg_slope.GetBinContent(b) / slope[b-1])
            bkg_slope.Add(histos["data"], -1)
            bkg_slope.SetLineColor(ROOT.kBlue)
            bkg_slope.Draw("hist same")
    
        #h_data_bkg.GetYaxis().SetRangeUser(y_min-150,y_max)
        
        leg_coords = 0.65,0.2,0.9,0.4
        if "legpos" in hsOpt:
            if hsOpt["legpos"] == "top":
                leg_coords = 0.65,0.78,0.9,1.
            elif hsOpt["legpos"] == "left" or hsOpt["legpos"] == "topleft":
                leg_coords = 0.1,0.78,0.35,1.
            elif hsOpt["legpos"] == "middle":
                leg_coords = 0.47,0.0,0.63,0.25
        leg = TLegend(*leg_coords)
        leg.SetTextSize(0.05)
        leg.AddEntry(h_data_bkg, "Data - fitted background", "p")
        leg.AddEntry(h_sig, "HH4b fitted")
        leg.AddEntry(hlist[0][-1], "HH4b fitted x5")
        leg.AddEntry(h_error, "Total uncertainty")
        leg.Draw("same")
        
        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.pdf")
        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.png")            
        c1.Clear()
        
    elif residuals == -12: # plot data - fitted (sig+bkg) in residual plot
        c1.cd()
        pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.4)
        pad2.SetTopMargin(0.)
        pad2.SetBottomMargin(0.2)
        pad2.Draw()
        pad2.cd()
        
        #hlist1, snames1, legstack1, hlist2, snames2, legstack2, hsOpt, residuals, norm, oDir, colors, dofill, rebin, headerOpt, isMC
        histos = {}
        histos["data"] = hlist[1][0]
        histos["sig"] = hlist[0][0]
        histos["bkg"] = hlist[0][1]
        #print hlist
        #for binn in range(1, histos["data"].GetNbinsX()+1):
        #    print binn, histos["data"].GetBinContent(binn), histos["sig"].GetBinContent(binn), histos["bkg"].GetBinContent(binn)
        (h_data_bkg, h_sig, h_error) = getHistosPostFitTemp(histos, hsOpt, snames[0], colors, fit_results, postfit_file)
        
        
        if True:
            slope = [0.98610, 0.98652, 0.98695, 0.98738, 0.98781, 0.98823, 0.98866, 0.98909, 0.98951, 0.98994, 0.99037, 0.99080, 0.99122, 0.99165, 0.99208, 0.99250, 0.99293, 0.99336, 0.99379, 0.99421, 0.99464, 0.99507, 0.99549, 0.99592, 0.99635, 0.99678, 0.99720, 0.99763, 0.99806, 0.99848, 0.99891, 0.99934, 0.99977, 1.00019, 1.00062, 1.00105, 1.00147, 1.00190, 1.00233, 1.00276, 1.00318, 1.00361, 1.00404, 1.00446, 1.00489, 1.00532, 1.00575, 1.00617, 1.00660, 1.00703, 1.00745, 1.00788, 1.00831, 1.00874, 1.00916, 1.00959, 1.01002, 1.01044, 1.01087, 1.01130, 1.01173, 1.01215, 1.01258, 1.01301, 1.01343, 1.01386, 1.01429, 1.01472, 1.01514, 1.01557, 1.01600, 1.01643, 1.01685, 1.01728, 1.01771, 1.01813, 1.01856, 1.01899, 1.01942, 1.01984, 1.02027, 1.02070, 1.02112, 1.02155, 1.02198, 1.02241, 1.02283, 1.02326, 1.02369, 1.02411, 1.02454, 1.02497, 1.02540, 1.02582, 1.02625, 1.02668, 1.02710, 1.02753, 1.02796, 1.02839,]
            print slope
            bkg_slope = histos["data"].Clone("bkg_slope")
            for b in range(1, bkg_slope.GetNbinsX()+1):
                bkg_slope.SetBinContent(b, bkg_slope.GetBinContent(b) / slope[b-1])
            bkg_slope.Add(histos["data"], -1)
            bkg_slope.SetLineColor(ROOT.kBlue)
            
        
        h_data_bkg.SetTitle("")
        h_data_bkg.GetXaxis().SetTitleSize(20)
        h_data_bkg.GetXaxis().SetTitleFont(43)
        h_data_bkg.GetXaxis().SetTitleOffset(4.)
        h_data_bkg.GetXaxis().SetLabelFont(43)
        h_data_bkg.GetXaxis().SetLabelSize(20)
        h_data_bkg.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        #h_data_bkg.GetYaxis().SetRangeUser(-20,20)
        minbin = h_data_bkg.GetXaxis().GetFirst()
        maxbin = h_data_bkg.GetXaxis().GetLast()
        y_max = max(h_data_bkg.GetMaximum(), hlist[0][-1].GetMaximum(), h_error.GetMaximum())*1.25
        y_min = min(h_data_bkg.GetMinimum(), hlist[0][-1].GetMinimum(), h_error.GetMinimum())*1.25

        h_data_bkg.GetYaxis().SetRangeUser(-60, 80)
        h_data_bkg.GetYaxis().SetTitleSize(20)
        h_data_bkg.GetYaxis().SetTitleFont(43)
        h_data_bkg.GetYaxis().SetTitleOffset(1.40)
        h_data_bkg.GetYaxis().SetLabelFont(43)
        h_data_bkg.GetYaxis().SetLabelSize(18)
        #h_data_bkg.SetMarkerStyle(33)
        h_data_bkg.SetMarkerStyle(20)
        h_data_bkg.SetMarkerSize(0.7)
        h_data_bkg.SetMarkerColor(1)
        h_data_bkg.SetLineColor(1)
        h_data_bkg.GetXaxis().SetTitle(hsOpt['xname'])
        #h_data_bkg.GetYaxis().SetTitle('data/bkg')
        h_data_bkg.Draw("e3 x0")
        #hlist[0][-1].Draw("hist same")
        
        h_sig.Draw("hist same")
        
        h_error.SetFillColor(632)
        h_error.Draw("E2same")
        h_data_bkg.Draw("same e3 x0")
        bkg_slope.Draw("hist same")
        #    
        #h_data_bkg.GetYaxis().SetRangeUser(y_min-150,y_max)
        
        leg_coords = 0.65,0.2,0.9,0.4
        if "legpos" in hsOpt:
            if hsOpt["legpos"] == "top":
                leg_coords = 0.65,0.78,0.9,1.
            elif hsOpt["legpos"] == "left" or hsOpt["legpos"] == "topleft":
                leg_coords = 0.1,0.78,0.35,1.
            elif hsOpt["legpos"] == "middle":
                leg_coords = 0.47,0.0,0.63,0.25
        leg = TLegend(*leg_coords)
        leg.SetTextSize(0.05)
        leg.AddEntry(h_data_bkg, "Data - fitted background", "p")
        leg.AddEntry(h_sig, "HH4b fitted")
        #leg.AddEntry(hlist[0][-1], "HH4b fitted x5")
        leg.AddEntry(h_error, "Only bias uncertainty")
        leg.Draw("same")
        
        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.pdf")
        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.png")            
        c1.Clear()
        
    elif residuals == -4 or residuals == -14: # plot data - fitted (sig+bkg) in residual plot
        c1.cd()
        pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.4)
        pad2.SetTopMargin(0.)
        pad2.SetBottomMargin(0.2)
        pad2.Draw()
        pad2.cd()
        pad2.SetTicky(0)
        
        #hlist1, snames1, legstack1, hlist2, snames2, legstack2, hsOpt, residuals, norm, oDir, colors, dofill, rebin, headerOpt, isMC
        histos = {}
        histos["data"] = hlist[1][0]
        histos["sig"] = hlist[0][0]
        histos["bkg"] = hlist[0][1]
        #print hlist
        #for binn in range(1, histos["data"].GetNbinsX()+1):
        #    print binn, histos["data"].GetBinContent(binn), histos["sig"].GetBinContent(binn), histos["bkg"].GetBinContent(binn)
        (h_data_bkg, h_sig, h_error) = getHistosPostFitRatio(histos, hsOpt, snames[0], colors, fit_results, postfit_file, 
                only_bias_unc = residuals == -14)
        #print h_data_bkg.Integral(), h_sig.Integral(), h_err.Integral()
        
        h_data_bkg.SetTitle("")
        h_data_bkg.GetXaxis().SetTitleSize(20)
        h_data_bkg.GetXaxis().SetTitleFont(43)
        h_data_bkg.GetXaxis().SetTitleOffset(4.)
        h_data_bkg.GetXaxis().SetLabelFont(43)
        h_data_bkg.GetXaxis().SetLabelSize(20)
        h_data_bkg.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        #h_data_bkg.GetYaxis().SetRangeUser(-20,20)
        minbin = h_data_bkg.GetXaxis().GetFirst()
        maxbin = h_data_bkg.GetXaxis().GetLast()
        
        max_error = 0.
        min_error = 0.
        hlist[0][-1].Divide(histos["bkg"])
        hlist[0][-1].SetLineWidth(2)
        
        for i in range(1, h_error.GetNbinsX()+1):
            if h_error.GetBinLowEdge(i+1) > hsOpt["xmax"]: continue
            if h_error.GetBinLowEdge(i+1) <= hsOpt["xmin"]: continue
            if h_error.GetBinContent(i) + h_error.GetBinError(i) > max_error: max_error = h_error.GetBinContent(i) + h_error.GetBinError(i)
            if h_error.GetBinContent(i) - h_error.GetBinError(i) < min_error: min_error = h_error.GetBinContent(i) - h_error.GetBinError(i)
            #print h_error.GetBinContent(i) + h_error.GetBinError(i), max_error, h_error.GetBinContent(i) - h_error.GetBinError(i), min_error
            
        #y_max = max(h_data_bkg.GetMaximum(), hlist[0][-1].GetMaximum(), h_error.GetMaximum())*1.5
        #y_min = min(h_data_bkg.GetMinimum(), hlist[0][-1].GetMinimum(), h_error.GetMinimum())*1.5
        y_max = max(max(h_data_bkg.GetMaximum(), h_sig.GetMaximum(), max_error)*1.6, hlist[0][-1].GetMaximum() * 1.05) 
        y_min = min(min(h_data_bkg.GetMinimum(), h_sig.GetMinimum(), min_error)*1.6, hlist[0][-1].GetMinimum() * 1.05) 

        if "classifier" in hsOpt["hname"]:
            y_max = max(h_data_bkg.GetMaximum(), h_sig.GetMaximum(), max_error)*1.1 
            y_min = min(h_data_bkg.GetMinimum(), h_sig.GetMinimum(), min_error)*1.1
            
        print h_data_bkg.GetMaximum(), h_sig.GetMaximum(), max_error, hlist[0][-1].GetMaximum()
        h_data_bkg.GetYaxis().SetRangeUser(y_min,y_max)
        h_data_bkg.GetYaxis().SetTitleSize(20)
        h_data_bkg.GetYaxis().SetTitleFont(43)
        h_data_bkg.GetYaxis().SetTitleOffset(1.40)
        h_data_bkg.GetYaxis().SetLabelFont(43)
        h_data_bkg.GetYaxis().SetLabelSize(18)
        #h_data_bkg.SetMarkerStyle(33)
        h_data_bkg.SetMarkerStyle(20)
        h_data_bkg.SetMarkerSize(0.7)
        h_data_bkg.SetMarkerColor(1)
        #h_data_bkg.SetLineColor(1)
        h_data_bkg.GetXaxis().SetTitle(hsOpt['xname'])
        #h_data_bkg.GetYaxis().SetTitle('data/bkg')
        h_data_bkg.SetLineWidth(0)         
    
        h_data_bkg.Draw("e1")
        h_data_bkg.GetYaxis().SetTitle("Normalized residuals")
        if not residuals == -14 and not "classifier" in hsOpt["hname"]:
            hlist[0][-1].Draw("hist same")
        
        h_sig.Draw("hist same")
        
        h_error.SetFillColor(600)
        h_error.Draw("E2same")
        h_data_bkg.Draw("same e1")
        
        
        leg_coords = 0.65,0.2,0.9,0.4
        #if "legpos" in hsOpt:
        #    if hsOpt["legpos"] == "top":
        #        leg_coords = 0.65,0.78,0.9,1.
        #    elif hsOpt["legpos"] == "left" or hsOpt["legpos"] == "topleft":
        leg_coords = 0.2,0.7,0.35,0.96
        if "legpos" in hsOpt:
            if hsOpt["legpos"] == "middle":
                leg_coords = 0.47,0.7,0.63,0.96
        leg = TLegend(*leg_coords)
        leg.SetTextSize(0.05)
        leg.AddEntry(h_data_bkg, "(Data - background) / background", "p")
        leg.AddEntry(h_sig, "(HH to 4b signal) / background")
        if not residuals == -14:
            if not "classifier" in hsOpt["hname"]:
                leg.AddEntry(hlist[0][-1], "HH4b fitted x5")
            leg.AddEntry(h_error, "Total uncertainty")
        else:
            leg.AddEntry(bkg_slope, "Thingie")
            leg.AddEntry(h_error, "Bias uncertainty")        
        leg.Draw("same")
        
        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.pdf")
        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.png")            
        c1.Clear()
        

    elif residuals == 5: # plot residuals before and after bias correction
        bias_range = (-5, 5)
        points = 10000
        #hlist1, snames1, legstack1, hlist2, snames2, legstack2, hsOpt, residuals, norm, oDir, colors, dofill, rebin, headerOpt, isMC
        histos = {}
        histos["data"] = hlist[1][0]
        histos["bkg"] = hlist[0][0]
        #print hlist
        #for binn in range(1, histos["data"].GetNbinsX()+1):
        #    print binn, histos["data"].GetBinContent(binn), histos["sig"].GetBinContent(binn), histos["bkg"].GetBinContent(binn)
        #(h_data_bkg, h_sig, h_error) = getHistosPostFit(histos, hsOpt, snames[0], colors, fit_results, postfit_file)
        #print h_data_bkg.Integral(), h_sig.Integral(), h_err.Integral()
        
        """h_data_bkg = histos["data"].Clone()
        h_data_bkg.Add(histos["bkg"], -1)
        
        h_data_bkg.SetTitle("")
        h_data_bkg.GetXaxis().SetTitleSize(20)
        h_data_bkg.GetXaxis().SetTitleFont(43)
        h_data_bkg.GetXaxis().SetTitleOffset(4.)
        h_data_bkg.GetXaxis().SetLabelFont(43)
        h_data_bkg.GetXaxis().SetLabelSize(20)
        h_data_bkg.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])"""
        #h_data_bkg.GetYaxis().SetRangeUser(-20,20)
        
        hres_nocorr = histos["data"].Clone("h_res_nocorr")
        hres_nocorr.Reset()
        hres_nocorr.SetMinimum(bias_range[0])
        hres_nocorr.SetMaximum(bias_range[1])
        
        residuals_nocorr = []
        
        h_bkg_norm = histos["bkg"].Clone("h_bkg_norm")
        h_bkg_norm.Scale(histos["data"].Integral()/histos["bkg"].Integral())
        
        for i in range(1, hres_nocorr.GetXaxis().GetNbins()+1):
            n1 = h_bkg_norm.GetBinContent(i)
            n2 = histos["data"].GetBinContent(i)
            e1 = h_bkg_norm.GetBinError(i)
            e2 = histos["data"].GetBinError(i)
            #print  "NoCorr: ",i, n1, n2, e1, e2
            if n1 and e1: 
                hres_nocorr.SetBinContent(i,(n1-n2)/math.sqrt(e1*e1+e2*e2)) #sign is fine!!!!!!!!
                err = (pow(n1,3) + 15*pow(n1,2)*n2+15*pow(n2,2)*n1 + pow(n2,3))/(4*pow((n1+n2),3))
                hres_nocorr.SetBinError(i, err)
            if hres_nocorr.IsBinOverflow(i) and hres_nocorr.GetBinContent(i) == 0: 
                print "zero"
                residuals_nocorr.append(0)
                continue
            else: 
                residuals_nocorr.append(hres_nocorr.GetBinContent(i))
            print hres_nocorr.GetBinError(i)
        #residuals_nocorr.append(0)
                
        #hres_nocorr.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        hres_nocorr.SetMarkerStyle(8)
        hres_nocorr.SetMarkerSize(0.8)
        hres_nocorr.SetMarkerColor(1)
        hres_nocorr.SetLineColor(1)
        hres_nocorr.GetXaxis().SetTitle(hsOpt['xname'])
        hres_nocorr.SetNdivisions(520, "X")
        hres_nocorr.GetXaxis().SetRangeUser(hsOpt['xmin'],1.02)
        #f1 =new TF1("f1","-x",-10,10);
        #TAxis* a = h->GetXaxis();
        #TGaxis *A1 = new TGaxis(0,2,10,2,"f1",510,"-");
   
        hres_nocorr.GetYaxis().SetTitle('Residuals (sigma units)')
        
        hres = hres_nocorr.Clone("h_res")
        residuals_corr = []
        h_bkg_bias_corr_norm = h_bkg_bias_corr.Clone("h_bkg_bias_corr_norm")
        
        vall = []
        for i in range(1, hres.GetXaxis().GetNbins()+1):
            print "yyy", i, h_bkg_bias_corr_norm.GetBinContent(i)
            vall.append(h_bkg_bias_corr_norm.GetBinContent(i))
            
        print "vall", vall
        h_bkg_bias_corr_norm.Scale(histos["data"].Integral() / h_bkg_bias_corr_norm.Integral())
        
        #for i in range(1, hres.GetXaxis().GetNbins()+1):
        #    print "yyy", i, h_bkg_bias_corr_norm.GetBinContent(i)
        
        for i in range(1, hres.GetXaxis().GetNbins()+1):
            n1 = h_bkg_bias_corr_norm.GetBinContent(i)
            n2 = histos["data"].GetBinContent(i)
            e1 = h_bkg_bias_corr_norm.GetBinError(i)
            e2 = histos["data"].GetBinError(i)
            #print  i, n1, n2, e1, e2, (n1-n2)/math.sqrt(e1*e1+e2*e2)
            if e1 > 0 or e2 > 0: 
                hres.SetBinContent(i,(n1-n2)/math.sqrt(e1*e1+e2*e2)) #sign is fine!!!!!!!!
                err = (pow(n1,3) + 15*pow(n1,2)*n2+15*pow(n2,2)*n1 + pow(n2,3))/(4*pow((n1+n2),3))
                hres.SetBinError(i, err)
            else:
                hres.SetBinContent(i, 0)
                hres.SetBinError(i, 0)
                #residuals_corr.append(0)
            if hres.IsBinOverflow(i) and hres.GetBinContent(i) == 0:
                print "zero"
                residuals_corr.append(0)
                continue            
            else: 
                residuals_corr.append(hres.GetBinContent(i))
            print "res", i, residuals_corr[-1]
            
        #residuals_corr.append(0)
        
        print residuals_corr
        c5 = TCanvas("tc5", "tc5", 200, 200)
        c5.cd()
        res_pull_nocorr = Hist(40, bias_range[0], bias_range[1], name = "h_res_pull")
        for v in residuals_nocorr: res_pull_nocorr.fill(v)
        _vals = []
        #for i in h_test: l_vals.append(i.value)
        #    print(l_vals)
        res_pull_nocorr.Fit("gaus", "ILLS")
        gaus_uncorr = res_pull_nocorr.GetFunction("gaus")
        #gaus_nocorr.Draw()
        res_pull_nocorr.Draw()
        
        
        res_pull = Hist(40, -5., 5, name = "h_res_pull")
        #res_pull.GetXaxis().SetTitle("residuals (sigma units)")
        for v in residuals_corr: res_pull.fill(v)
        res_pull.Fit("gaus", "ILLS")

        for i in range(0, res_pull.GetNbinsX()+1):
            _vals.append(res_pull.GetBinContent(i))
            #print i, res_pull.GetBinContent(i)
        print _vals, sum(_vals)        
        
        
        gaus_corr = res_pull.GetFunction("gaus")
        #test.FillRandom("gaus",2000)
        #test.SetFillColor(2)"""
        
        
        mean_uncorr = (gaus_uncorr.GetParameter(1), gaus_uncorr.GetParError(1))
        sigma_uncorr = (gaus_uncorr.GetParameter(2), gaus_uncorr.GetParError(2))
        
        mean_corr = (gaus_corr.GetParameter(1), gaus_corr.GetParError(1))
        sigma_corr = (gaus_corr.GetParameter(2), gaus_corr.GetParError(2))
        
        c1.cd()
        pad2 = TPad("pad2", "pad2", 0.5, 0.5, 0.9, 1)
        pad2.SetTopMargin(0.1)
        #pad2.SetBottomMargin(0.1)
        pad2.SetRightMargin(0.)
        pad2.Draw()
        pad2.cd()
        
        hres_nocorr.Draw("E X0")
        line = TLine(hsOpt['xmin'],0.,hsOpt['xmax'],0.);
        line.Draw()
        
        #test.Draw("hbar same")
        pad2.RedrawAxis()
        
        
        c1.cd()
        pad3 = TPad("pad3", "pad3", 0.5, 0., 0.9, 0.5)
        pad3.SetRightMargin(0.)
        #pad3.SetTopMargin(0.1)
        pad3.SetBottomMargin(0.18)
        pad3.Draw()
        pad3.cd()
        
        hres.Draw("E X0")
        line2 = TLine(hsOpt['xmin'],0.,hsOpt['xmax'],0.);
        line2.Draw()
        
        
        c1.cd()
        pad_g1 = TPad("pad4", "pad4", 0.9, 0.5, 1, 1.)
        pad_g1.SetLeftMargin(0.)
        pad_g1.SetTopMargin(0.1)
        #pad_g1.SetBottomMargin(0.1)
        pad_g1.SetTickx(0)
        pad_g1.Draw()
        pad_g1.cd()
        
        empty0 = Hist(10, 0., 10, name = "empty0")
        empty0.SetMinimum(bias_range[0])
        empty0.SetMaximum(bias_range[1])
        empty0.GetXaxis().SetLabelSize(0.)
        empty0.GetYaxis().SetLabelSize(0.)
        empty0.Draw()
        
        xs_nocorr = []
        
        x_range = np.linspace(-8,8, points)
        for y in x_range:
            xs_nocorr.append(gaus_uncorr.Eval(y))
            
        
        #pline = TPolyLine(1000,np.array(xs), x_range)
        pline_nocorr = TGraph(10000, np.array(xs_nocorr), x_range)
        #pline.SetFillColor(38)
        pline_nocorr.SetFillStyle(0)
        pline_nocorr.SetLineColor(2)
        pline_nocorr.SetLineWidth(1)
        #pline_nocorr.Draw("f same")
        pline_nocorr.Draw("same")
        
        text_mean = "Mean: %.2f #pm %.2f" % mean_uncorr
        text_sigma = "Sigma: %.2f #pm %.2f" % sigma_uncorr
        latex.SetTextFont(42)
        latex.SetTextAlign(12)
        latex.SetTextSize(0.08)
        latex.DrawLatex(0.25, 0.85, "#bf{Before correction:}")
        latex.DrawLatex(0.25, 0.78, text_mean)
        latex.DrawLatex(0.25, 0.71, text_sigma)
        
        c1.cd()
        pad_g2 = TPad("pad5", "pad5", 0.9, 0., 1, 0.5)
        pad_g2.SetLeftMargin(0.)
        #pad_g2.SetTopMargin(0.1)
        pad_g2.SetBottomMargin(0.18)
        pad_g2.SetTickx(0)
        pad_g2.Draw()
        pad_g2.cd()
        
        empty = Hist(10, 0., 10, name = "empty")
        empty.SetMinimum(bias_range[0])
        empty.SetMaximum(bias_range[1])
        empty.GetXaxis().SetLabelSize(0.)
        empty.GetYaxis().SetLabelSize(0.)
        empty.Draw()
        
        xs = []
        
        for y in x_range:
            xs.append(gaus_corr.Eval(y))
            
        
        pline = TGraph(10000, np.array(xs), x_range)
        pline.SetFillStyle(0)
        pline.SetLineColor(2)
        pline.SetLineWidth(1)
        #pline.Draw("f same")
        pline.Draw("same")
        
        text_mean = "Mean: %.2f #pm %.2f" % mean_corr
        text_sigma = "Sigma: %.2f #pm %.2f" % sigma_corr
        #latex = TLatex()
        #latex.SetNDC()
        #latex.SetTextSize(0.035)
        #latex.SetTextColor(1)
        latex.SetTextFont(42)
        latex.SetTextAlign(12)
        latex.SetTextSize(0.08)
        latex.DrawLatex(0.25, 0.9, "#bf{After correction:}")
        latex.DrawLatex(0.25, 0.83, text_mean)
        latex.DrawLatex(0.25, 0.76, text_sigma)
        
        
        #gStyle.SetOptStat(111)
        #gStyle.SetOptFit(1)
        
        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.svg")
        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.pdf")
        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.png")
        #c1.Clear()
        
    elif residuals == 10:  # to get comparison of hres with default bias
        c3 = TCanvas("c3", "comp_res"+hsOpt['hname'], 1000, 1000)
        pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0.0) 
        pad1.Draw()             
        pad1.cd()
        #c3.Divide(1,2)
        #c3.cd(1)        
        hres = h[1].Clone("h_res")        
        checkbin = False
        for i in range(1, hres.GetXaxis().GetNbins()+1):
            n1 = h[0].GetBinContent(i)
            n2 = h[1].GetBinContent(i)
            e1 = h[0].GetBinError(i)
            e2 = h[1].GetBinError(i)
            if n1 :#and e1: 
                hres.SetBinContent(i,(n2-n1)) # order is correct!! bkg - truth
                hres.SetBinError(i, math.sqrt(pow(e1,2)+pow(e2,2)))
        hres.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        hres.SetMarkerStyle(8)
        hres.SetMarkerSize(0.8)
        hres.SetMarkerColor(1)
        hres.SetLineColor(1)
        hres.GetXaxis().SetTitle(hsOpt['xname'])
        hres.GetYaxis().SetTitle('Events')

        bkg_bias_fname = "/lustre/cmswork/dcastrom/projects/hh/april_2017/CMSSW_8_0_25/src/Analysis/hh2bbbb_limit/notebooks/bias_21022018_with_weights/BM%d/bias_correction_bigset_unscaled.json" % bm
        with open(bkg_bias_fname,"r") as bkg_bias_file:
            json_dict = json.load(bkg_bias_file)
            print ("using bias file: ", bkg_bias_file)

        hbias = hres.Clone("h_bias")
        for n in range(len(json_dict['bias_corr'])):
                #if not s_bin.overflow:
                #value = s_bin.value

                bias = json_dict['bias_corr'][n]
                var = json_dict['var'][n]
                bias_unc = json_dict['bias_corr_unc_bs'][n]
                bias_unc_stat = json_dict['bias_corr_unc_stat'][n]
                
                new_bkg_pred_stat = var
                                
                new_bkg_pred_tot_unc = np.sqrt(new_bkg_pred_stat**2 + bias_unc**2 + bias_unc_stat**2)
                hbias.SetBinContent(n+1, bias)
                hbias.SetBinError(n+1, new_bkg_pred_tot_unc)

        
        for i in range(1, hres.GetXaxis().GetNbins()+1):
            n1 = hres.GetBinContent(i)
            n2 = hbias.GetBinContent(i)
            e1 = hres.GetBinError(i)
            e2 = hbias.GetBinError(i)
        
        ymax = 650.
        ymin = -200.
        if hbias.GetMaximum() > ymax: ymax = hbias.GetMaximum()*1.15
        if hres.GetMaximum() > ymax: ymax = hres.GetMaximum()*1.15
        if hbias.GetMinimum() < ymin: ymin = hbias.GetMinimum()*1.15
        if hres.GetMinimum() < ymin: ymin = hres.GetMinimum()*1.15
        hres.GetYaxis().SetRangeUser(ymin,ymax)
        hbias.GetYaxis().SetRangeUser(ymin,ymax)

        #print "res", hres.Integral(), "bias", hbias.Integral()
        hbias.SetMarkerColor(ROOT.kRed)
        hbias.SetLineColor(2)
        hbias.Draw("E X0")
        hres.Draw("E X0 same")
        
        legend = setLegend(1,1)
        legend.AddEntry(hres, "tt mixed - tt MC", "p")
        legend.AddEntry(hbias, "Bias", "p")
        legend.Draw("same")
        
        #Residual panel
        c3.cd()
        pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0.)
        pad2.SetBottomMargin(0.2)
        pad2.Draw()
        pad2.cd()

        hresidual = hres.Clone("h_residual")
        hresidual.Add(hbias, -1)
        # consider only data error in the ratio plot
        for ibin in range(1, hresidual.GetNbinsX()+1):
            #print "A", ibin, hresidual.GetBinContent(ibin), hres.GetBinContent(ibin), hbias.GetBinContent(ibin)
            hresidual.SetBinContent(ibin, hresidual.GetBinContent(ibin) / math.sqrt(hres.GetBinError(ibin)**2 + hbias.GetBinError(ibin)**2) )
            hresidual.SetBinError(ibin, 0.01)
            #print ibin, hresidual.GetBinContent(ibin), hres.GetBinError(ibin), hbias.GetBinError(ibin)
            #if(h2.GetBinContent(ibin)>0): 
            #    #hrat.SetBinError(ibin, h2.GetBinError(ibin)/h2.GetBinContent(ibin) )
            #else: 
            #    hrat.SetBinError(ibin, 0.)"""

        # MC uncertainy shadow plot
        h_error = herr[0].Clone("h_err")
        h_error.Reset()
        for ibin in range(1, herr[0].GetNbinsX()+1):
                h_error.SetBinContent (ibin, 0.)
                h_error.SetBinError   (ibin, 1.)

        hresidual.SetTitle("")
        hresidual.GetXaxis().SetTitleSize(20)
        hresidual.GetXaxis().SetTitleFont(43)
        hresidual.GetXaxis().SetTitleOffset(4.)
        hresidual.GetXaxis().SetLabelFont(43)
        hresidual.GetXaxis().SetLabelSize(20)
        hresidual.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])

        #hresidual.GetYaxis().SetRangeUser(min(h_error.GetMinimum(), hresidual.GetMinimum())*1.15,max(h_error.GetMaximum(), hresidual.GetMaximum())*1.15)
        hresidual.GetYaxis().SetRangeUser(-3, 3)
        hresidual.GetYaxis().SetTitleSize(20)
        hresidual.GetYaxis().SetTitleFont(43)
        hresidual.GetYaxis().SetTitleOffset(1.40)
        hresidual.GetYaxis().SetLabelFont(43)
        hresidual.GetYaxis().SetLabelSize(18)
        hresidual.SetMarkerStyle(8)
        hresidual.SetMarkerSize(0.8)
        hresidual.SetMarkerColor(1)
        hresidual.SetLineColor(1)
        hresidual.GetXaxis().SetTitle("BDT output")
        hresidual.GetYaxis().SetTitle('difference in #sigma units')
        hresidual.Draw("E1 X0")
        h_error.GetXaxis().SetTitle("")
        h_error.SetFillColor(430)
        h_error.Draw("E2same")
    
        hresidual.Fit("pol5")
        myfunc = hresidual.GetFunction("pol5")
        print myfunc
        fithist = myfunc.CreateHistogram()
        #for binn in range(1, fithist.GetNbinsX()):
        #    print (binn, fithist.GetBinContent(binn))
        #print myfunc.Integral(0., 74/80.), myfunc.Integral(76./80, 1.)*80/4., myfunc.Integral(0., 1.), myfunc.Integral(0., 2.)
        squares = 0.
        squares_0_2 = 0.
        squares_70p = 0.
        for binn in range(1, 81):
            #print (binn, myfunc.Eval((binn-0.5)/80.), myfunc.Integral((binn - 1.) / 80., binn / 80.))
            #print myfunc.Eval((binn-0.5)/80.), ", ",
            """if abs(hres.GetBinContent(binn)) > abs(hbias.GetBinContent(binn)):
                print hres.GetBinContent(binn), ", ",
            else:     
                print hbias.GetBinContent(binn), ", ","""
            if binn > 0:
                squares += myfunc.Eval((binn-0.5)/80.)**2
            if binn > 16:
                squares_0_2 += myfunc.Eval((binn-0.5)/80.)**2
            if binn > 72:
                squares_70p += myfunc.Eval((binn-0.5)/80.)**2
        print
        print "Tabel: BM%d &    & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f" % (bm, myfunc.Integral(0., 1.), myfunc.Integral(0.2, 1.), myfunc.Integral(0.9, 1.), squares, squares_0_2, squares_70p)
        
        """for binn in range(1, 81):
            #print (binn, myfunc.Eval((binn-0.5)/80.), myfunc.Integral((binn - 1.) / 80., binn / 80.))
            #print myfunc.Eval((binn-0.5)/80.), ", ",
            print max(abs(hres.GetBinContent(binn)), abs(hbias.GetBinContent(binn)), hbias.GetBinError(binn)), ", ",
        print """

        c3.SaveAs(oDir+"/"+hsOpt['hname']+"_resc.pdf")
        c3.SaveAs(oDir+"/"+hsOpt['hname']+"_resc.png")            
        c3.SaveAs(oDir+"/"+hsOpt['hname']+"_resc.root")            

    elif residuals >=1: # to get residuals and pulls -- data as h2!!!
        if residuals==2:
            c2 = TCanvas("c2", "res"+hsOpt['hname'], 800, 800)    
            c2.Divide(1,2)
            c2.cd(1)
        else:
            c2 = TCanvas("c2", "res"+hsOpt['hname'], 800, 400)
        hres = h[1].Clone("h_res")
        checkbin = False
        for i in range(1, hres.GetXaxis().GetNbins()+1):
            n1 = h[0].GetBinContent(i)
            n2 = h[1].GetBinContent(i)
            e1 = h[0].GetBinError(i)
            e2 = h[1].GetBinError(i)
            #print  i, n1, n2, e1, e2
            if n1 and e1: 
                hres.SetBinContent(i,(n1-n2)/math.sqrt(e1*e1+e2*e2)) #sign is fine!!!!!!!!
                err = (pow(n1,3) + 15*pow(n1,2)*n2+15*pow(n2,2)*n1 + pow(n2,3))/(4*pow((n1+n2),3))
                hres.SetBinError(i, err)


        hres.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        hres.SetMarkerStyle(8)
        hres.SetMarkerSize(0.8)
        hres.SetMarkerColor(1)
        hres.SetLineColor(1)
        hres.GetXaxis().SetTitle(hsOpt['xname'])
        hres.GetYaxis().SetTitle('residuals (sigma units)')
        hres.Draw("E X0")
        line = TLine(hsOpt['xmin'],0.,hsOpt['xmax'],0.);
        line.Draw()
        if residuals>1: #draw also pull with fit
            c2.cd(2)
            gStyle.SetOptStat(111)
            gStyle.SetOptFit(1)
            res_a = []
            for i in range(1, hres.GetXaxis().GetNbins()+1):                
                if (hres.IsBinOverflow(i) and hres.GetBinContent(i) == 0): continue
                else:  res_a.append(hres.GetBinContent(i))
            res_pull = Hist(40, -5., 5, name = "h_res_pull")
            res_pull.GetXaxis().SetTitle("residuals (sigma units)")
            for v in res_a: res_pull.fill(v)
            res_pull.Fit("gaus", "LL")
            res_pull.Draw("E X0")
            text = "BM%d" % bm
            text = text.replace("BM0", "SM").replace("BM13", "Box")
            latex.DrawLatex(0.18, 0.88, text)
            latex.SetTextFont(52)
            latex.SetTextSize(0.04)

        c2.Update()    
        c2.SaveAs(oDir+"/"+hsOpt['hname']+"_res.pdf")
        c2.SaveAs(oDir+"/"+hsOpt['hname']+"_res.png")            
        c2.SaveAs(oDir+"/"+hsOpt['hname']+"_res.svg")
        #c2.SaveAs(oDir+"/"+hsOpt['hname']+"_res.root")   

    return [nevs, neverrs]
#-----------

def drawH1tdrAcc(hs, snames, leg,
                   oDir, colors, headerOpt, isMC):
    gStyle.SetOptStat(False)
    c1 = TCanvas("c1", 'acc', 800, 800)
    markers = [23,24,25]

    legend = TLegend(0.50,0.7,0.85,0.85)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)

    for h_ in hs[0]:
        hd =  h_.Clone('hd')    
    for i, h_ in enumerate(hs):
        if i==0: 
            continue
        for h__ in h_:
            h =  h__.Clone('h')
        h.SetTitle("")
        h.Divide(hd) 
        h.GetXaxis().SetTitle('generator p_{T} (GeV)')
        h.GetYaxis().SetTitle('acceptance')
        h.GetYaxis().SetLabelOffset(0.010)
        h.GetYaxis().SetTitleOffset(1.45)
        h.GetXaxis().SetLabelOffset(0.010)
        h.GetXaxis().SetTitleOffset(1.45)
        h.SetMaximum(1.15)
        h.GetXaxis().SetRangeUser(0.,100.)
        h.SetMarkerStyle(markers[i-1])
        h.SetMarkerSize(1.)
        h.SetMarkerColor(colors[i-1])
        legend.AddEntry(h, leg[i-1], "P")
        if i==1: 
            h.Draw("PE")
            #Legend & headers
            drawCMStdr(headerOpt)

        h.Draw("PE same") #draw first twice
    legend.Draw("same")
    c1.Update()
    c1.SaveAs(oDir+"/acc.pdf")
    c1.SaveAs(oDir+"/acc.png")


def drawH1tdr(hlist1, snames1, leg1, hsOpt, residuals, norm, oDir, colors, dofill, rebin, headerOpt, isMC):
    gStyle.SetOptStat(False)
    c1 = TCanvas("c1", hsOpt['hname'], 800, 800)       

    if rebin > 0: rb = rebin
    else: rb =  hsOpt['rebin']

    ymax = 0.
    nh1 = hlist1[0].Integral()  

    herr = hlist1[0].Clone("hs_error")  
    herr.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])

    herr.Rebin(rb)
    herr.GetXaxis().SetTitle(hsOpt['xname'])
    herr.GetYaxis().SetTitle(hsOpt['yname'])
    herr.SetMinimum(0.)
    herr.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
    herr.SetMinimum(0.)
    herr.Draw()

    #Legend & headers
    if isMC: drawCMS(-1, headerOpt)
    else: drawCMS(35.9, headerOpt)
    legend = setLegend(1,1)
    if len(leg1) > 1:
        for n, leg in enumerate(leg1):
            legend.AddEntry(hlist1[n], leg)
    else:
        legend.AddEntry(hlist1[len(hlist1)-1], leg1[0])
    legend.Draw("same")
    #-------------

    scale = 1.
    for i, h in enumerate(hlist1):         
        if(norm): scale = nh1/h.Integral()
        col = colors[i]
        h.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        h.SetMinimum(0.)
        h.Scale(scale)
        h.Rebin(rb)
        h.SetMarkerStyle(8)
        h.SetMarkerSize(0.)
        if dofill[i]: 
            h.SetFillColorAlpha(col,0.1)
            #h.SetFillColor(col)
            h.SetFillStyle(1)
        else:     
            h.SetLineWidth(1)         
            h.SetLineColorAlpha(col,0.)
            #h.SetLineColor(col)
        if i == 0: ymax = h.GetMaximum()*1.15
        h.SetMaximum(ymax)
        if i ==0: h.Draw("HISTEsame") 
        else: h.Draw("HISTEsame")

    nev1 = 0
    nev2 = 0
    nev1err = 0
    nev2err = 0

    c1.Update()    
    c1.SaveAs(oDir+"/"+hsOpt['hname']+".pdf")
    c1.SaveAs(oDir+"/"+hsOpt['hname']+".png")            
   #c1.SaveAs(oDir+"/"+hsOpt['hname']+".root") 

    return [nev1,nev1err,nev2,nev2err]
#------------


def drawH1only(hlist1, snames1, legstack1, hsOpt, oDir, colors, dofill, rebin, headerOpt, isMC):
    gStyle.SetOptStat(False)
    gStyle.SetOptTitle(0);
    c1 = TCanvas("c1", hsOpt['hname'], 800, 800)       

    if rebin > 0: rb = rebin
    else: rb =  hsOpt['rebin']

    ymax = 0.
    hs1, herr1, h1 =  getStackH(hlist1, hsOpt, rb, snames1, colors[0], 1., dofill[0])
    ymax = hs1.GetMaximum()*1.15

   #debug -- needed before drawing hs
    herr1.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
    herr1.SetMinimum(0.)
    herr1.SetMaximum(ymax)
    if len(hlist1)>1: 
        herr1.Draw("E2")
        herr1.SetFillColor(922)
    else: herr1.Draw("E")
   #--

    hs1.SetMaximum(ymax)
    herr1.SetMaximum(ymax)

    #Legend & headers
    if isMC: drawCMS(-1, headerOpt)
    else: drawCMS(35.9, headerOpt)
    legend = setLegend(1,1)
    nskip = 0
    match = False
    for n, sam in enumerate(snames1):
        # to get one legend for all HT bins - to be implemented for other samples         
        if sam.find("QCD")>=0 and match: 
            nskip+=1 
            continue
        if len(legstack1) > n-nskip:
            if sam.find("QCD")>=0: match = True
            legend.AddEntry(hlist1[n], legstack1[n-nskip])            
    if len(hlist1)>1: legend.AddEntry(herr1, 'bkg. unc. (stat.only)')
    legend.Draw("same")
    #-------------

    plotH(hlist1, hs1, herr1, dofill[0])
    nev1 = 0
    nev1err = 0
    if hsOpt['hname']=='h_nevts': 
        nev1 = herr1.GetBinContent(1)
        nev1err = herr1.GetBinError(1)

    c1.Update()    
    c1.SaveAs(oDir+"/"+hsOpt['hname']+".pdf")
    c1.SaveAs(oDir+"/"+hsOpt['hname']+".png")            
   #c1.SaveAs(oDir+"/"+hsOpt['hname']+".root") 

    return [nev1,nev1err]
#-------------------------


def fillH2(h,h1):
    h2 = (TH2D)("h2_bdt","h2_bdt",h.GetNbinsX(),0.,1.,h.GetNbinsX(),0.,1.) 
    for i in range(1,h.GetNbinsX()):
        for j in range(1,h1.GetNbinsX()):
         h2.SetBinContent(i, j, h.GetBinContent(i)+h1.GetBinContent(j))
    return h2

def drawH2(hs, hs1, hsOpt, sname, rebin, oDir, legs):
    gStyle.SetOptStat(False)
    legend = setLegend(1,1)

    if len(hs1):
        hs2 = []
        for h in hs:
            for h1 in hs1:
                hs2.append(fillH2(h1,h))
    else:
        hs2 = hs

    for i, h2 in enumerate(hs2):
        name = hsOpt['hname']
        print name
        c2 = TCanvas(name, name, 800, 800)
        h2.Rebin2D(rebin, rebin)
        h2.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        h2.GetYaxis().SetRangeUser(hsOpt['ymin'],hsOpt['ymax'])
        h2.GetXaxis().SetTitle(hsOpt['xname']+" "+str(legs[1]))
        h2.GetYaxis().SetTitle(hsOpt['yname']+" "+str(legs[0]))

        h2.GetXaxis().SetLabelFont(43)
        h2.GetXaxis().SetLabelSize(20)
        h2.GetXaxis().SetTitleFont(43)
        h2.GetXaxis().SetTitleOffset(1.)
        h2.GetXaxis().SetTitleSize(23) 
        h2.GetYaxis().SetLabelFont(43)
        h2.GetYaxis().SetLabelSize(20)
        h2.GetYaxis().SetTitleFont(43)
        h2.GetYaxis().SetTitleOffset(1.2)
        h2.GetYaxis().SetTitleSize(23)    
 
        h2.Draw("COLZ")
        c2.Update() 
        palette = h2.FindObject("palette")
        c2.Update() 

        c2.SaveAs(oDir+"/"+name+".pdf")
        c2.SaveAs(oDir+"/"+name+".png")            
        c2.SaveAs(oDir+"/"+name+".root")  

    return True
#------------

# DRAWING TOOLS
#################
def plotH(hlist, h, herr, fill, res = None):
    if fill:         
        if len(hlist)>1: 
            h.Draw("HISTsame")
            herr.Draw("E2same")
        elif res == 5:  #No error bars
            h.Draw("HISTsame")
        else:
            h.Draw("HISTEsame")  
    else:
        if len(hlist)>1: 
            h.Draw("HISTsame")   
            herr.Draw("E2same")
        else: 
            h.Draw("Esame")

def drawCMS(lumi, text, onTop=False ):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextAlign(31)     
    if (type(lumi) is float or type(lumi) is int) and float(lumi) > 0: latex.DrawLatex(0.90, 0.94, "%.1f fb^{-1}  (13 TeV)" % (float(lumi)))
    else: latex.DrawLatex(0.90, 0.94, "simulation (13 TeV)")
    if not onTop: latex.SetTextAlign(11)
    latex.SetTextFont(61)
    latex.SetTextSize(0.03 if len(text)>0 else 0.035)
    if not onTop: 
        latex.DrawLatex(0.2, 0.855, "CMS   "+text)
        latex.SetTextFont(52)
        latex.SetTextSize(0.02)
        #latex.DrawLatex(0.15, 0.830, "preliminary")
    else: 
        latex.DrawLatex(0.15, 0.88, "CMS   "+text)
        latex.SetTextFont(52)
        latex.SetTextSize(0.02)
        #latex.DrawLatex(0.15, 0.850, "preliminary")
#------------

def drawCMStdr(text, onTop=False):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextAlign(33)
    latex.DrawLatex(0.90, 0.94, "14 TeV")
    if not onTop: latex.SetTextAlign(11)
    latex.SetTextFont(62)
    latex.SetTextSize(0.03 if len(text)>0 else 0.035)
    if not onTop:
        latex.DrawLatex(0.15, 0.85, "CMS Phase-2   "+text)
        latex.SetTextFont(52)
        latex.SetTextSize(0.030)
        latex.DrawLatex(0.15, 0.820, "Simulation. Preliminary")
        latex.SetTextFont(52);
        latex.SetTextSize(0.030);
        latex.DrawLatex(0.15, 0.780, "pp #rightarrow HH #rightarrow b#bar{b}b#bar{b}");

def setLegend(x0 = 0.55, y0 = 0.6, x1 = 0.9, y1 = 0.85):
    leg = TLegend(x0, y0, x1, y1)
    leg.SetLineWidth(0)
    leg.SetTextSize(0.032)
    return leg
#------------

def drawPostFitH1(hlist1, snames1, legstack1, hlist2, snames2, legstack2, hsOpt, residuals, norm, oDir, colors, dofill, rebin, headerOpt, isMC):
    gStyle.SetOptStat(False)
    gStyle.SetOptTitle(0);

    c1 = TCanvas("c1", hsOpt['hname'], 800, 800)       
    if residuals==-1:
        pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0.03) 
        pad1.Draw()             
        pad1.cd()               
    
    isNevts=False
    if hsOpt['hname']=="h_nevts": isNevts=True
    ymax = 0.

    hs1, herr1, h1 =  getStackH(hlist1, hsOpt, rb, snames1, colors[0], scale1, dofill[0])

    if hs1.GetMaximum() > ymax: ymax = hs1.GetMaximum()*1.15
    if herr1.GetMaximum() > ymax: ymax = herr1.GetMaximum()*1.15
    if isNevts:  print "h1Int ",  h1.GetBinContent(1), h1.GetBinError(1)

    if(norm): scale2 = getScale(hsOpt,hlist2, hlist2)
    else: scale2 = 1.
    if scale2 != 1.: print "sc_to_norm2: ",scale2
    hs2, herr2, h2 =  getStackH(hlist2, hsOpt, rb, snames2, colors[1], scale2, dofill[1])
    if hs2.GetMaximum() > ymax: ymax = hs2.GetMaximum()*1.15
    if herr2.GetMaximum() > ymax: ymax = herr2.GetMaximum()*1.15
    if isNevts:  print "h2Int ",  h2.GetBinContent(1), h2.GetBinError(1)

    if len(hlist1) == 1 and len(hlist2) == 1:
        print("KS: ", hlist1[0].KolmogorovTest(hlist2[0]))
        print("Chi2: ", hlist1[0].Chi2Test(hlist2[0], "UU NORM"))
    else:
        print("Not doing KS test, stacks as input")
   #debug -- needed before drawing hs
    herr1.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
    herr1.SetMinimum(0.)
    herr1.SetMaximum(ymax)
    if len(hlist1)>1: 
        herr1.Draw("E2")
        herr1.SetFillColor(922)
    else: herr1.Draw("E")
   #--

    #Legend & headers
    if isMC: drawCMS(-1, headerOpt)
    else: drawCMS(35.9, headerOpt)
    legend = setLegend(1,1)
    legend.AddEntry(hlist2[len(hlist2)-1], legstack2[0]) #debug - one leg for second samplelist
    nskip = 0
    match = False 
    for n, sam in enumerate(snames1):
        # to get one legend for all HT bins - to be implemented for other samples         
        if sam.find("QCD")>=0 and match: 
            nskip+=1 
            continue
        if len(legstack1) > n-nskip:
            if sam.find("QCD")>=0: match = True
            legend.AddEntry(hlist1[n], legstack1[n-nskip])            
    if len(hlist1)>1: legend.AddEntry(herr1, 'bkg. unc. (stat.only)')
    legend.Draw("same")
    #-------------

    if(ymax > 1000): TGaxis.SetMaxDigits(3)
    hs1.SetMaximum(ymax)
    herr1.SetMaximum(ymax)
    plotH(hlist1, hs1, herr1, dofill[0])
    hs2.SetMaximum(ymax)
    herr2.SetMaximum(ymax)
    plotH(hlist2, hs2, herr2, dofill[1])
    herr2.Draw("Esameaxis")

    nev1 = 0
    nev2 = 0
    nev1err = 0
    nev2err = 0
    if isNevts: 
        nev1 = herr1.GetBinContent(1)
        nev2 = herr2.GetBinContent(1)
        nev1err = herr1.GetBinError(1)
        nev2err = herr2.GetBinError(1)

    herr1.GetXaxis().SetLabelSize(0.)

    c1.cd()
    pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(0.2)
    pad2.Draw()
    pad2.cd()
    hrat = h2.Clone("h_rat")
    hrat.Divide(h1)
    # consider only data error in the ratio plot
    for ibin in range(1, h2.GetNbinsX()+1):
	    #print(h2.GetBinContent(ibin))
	    if(h2.GetBinContent(ibin)>0): 
		hrat.SetBinError(ibin, h2.GetBinError(ibin)/h2.GetBinContent(ibin) )
	    else: 
		hrat.SetBinError(ibin, 0.)

    # MC uncertainy shadow plot
    herr = herr1.Clone("h_err")
    herr.Reset()
    for ibin in range(1, herr1.GetNbinsX()+1):
	    if(herr1.GetBinContent(ibin)>0): 
		herr.SetBinContent (ibin, 1.)
		herr.SetBinError   (ibin, herr1.GetBinError(ibin)/herr1.GetBinContent(ibin))
	    else: 
		herr.SetBinContent(ibin, 1.)
		herr.SetBinError   (ibin, 0.)
    
    hrat.SetTitle("")
    hrat.GetXaxis().SetTitleSize(20)
    hrat.GetXaxis().SetTitleFont(43)
    hrat.GetXaxis().SetTitleOffset(4.)
    hrat.GetXaxis().SetLabelFont(43)
    hrat.GetXaxis().SetLabelSize(20)
    hrat.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
    minbin = hrat.GetXaxis().GetFirst()
    maxbin = hrat.GetXaxis().GetLast()
    ymax_ = 1.5
    ymax = 1.
    ymin_ = 0.5
    ymin = 1.

    for ibin in range(minbin, maxbin+1):       
	    binc = hrat.GetBinContent(ibin) 
	    if binc == 0. or (binc is None): 
		continue
	    if binc*1.05 > ymax: 
		ymax = binc*1.05
	    if binc*0.95 < ymin: 
		ymin = binc*0.95

    hrat.GetYaxis().SetRangeUser(ymin,ymax)
    hrat.GetYaxis().SetTitleSize(20)
    hrat.GetYaxis().SetTitleFont(43)
    hrat.GetYaxis().SetTitleOffset(1.40)
    hrat.GetYaxis().SetLabelFont(43)
    hrat.GetYaxis().SetLabelSize(18)
    hrat.SetMarkerStyle(8)
    hrat.SetMarkerSize(0.8)
    hrat.SetMarkerColor(1)
    hrat.SetLineColor(1)
    hrat.GetXaxis().SetTitle(hsOpt['xname'])
    hrat.GetYaxis().SetTitle('data/bkg')
    hrat.Draw("E X0")
    herr.SetFillColor(430)
    herr.Draw("E2same")

    l = TLine(hsOpt['xmin'],1.5,hsOpt['xmax'],1.5)
    l0 = TLine(hsOpt['xmin'],1.4,hsOpt['xmax'],1.4)
    l00 = TLine(hsOpt['xmin'],1.3,hsOpt['xmax'],1.3)
    l000 = TLine(hsOpt['xmin'],1.2,hsOpt['xmax'],1.2)
    l1 = TLine(hsOpt['xmin'],1.1,hsOpt['xmax'],1.1)
    l2 = TLine(hsOpt['xmin'],1.,hsOpt['xmax'],1.)
    l3 = TLine(hsOpt['xmin'],0.9,hsOpt['xmax'],0.9)
    l4 = TLine(hsOpt['xmin'],0.8,hsOpt['xmax'],0.8)
    l.SetLineStyle(3)
    l0.SetLineStyle(3)
    l00.SetLineStyle(3)
    l000.SetLineStyle(3)
    l1.SetLineStyle(3)
    l2.SetLineStyle(3)
    l3.SetLineStyle(3)
    l4.SetLineStyle(3)
    #l.Draw("same")
    #l0.Draw("same")
    #l00.Draw("same")
    #l000.Draw("same")
    #l1.Draw("same")
    l2.Draw("same")
    #if ymin<0.9: l3.Draw("same")
    #if ymin<0.8: l4.Draw("same")
    
    c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.pdf")
    c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.png")
    c1.Clear()

    return [nev1,nev1err,nev2,nev2err]


#def getHistosPostFit(histos, hsOpt, rebin, snames, color, scale, fill):
def getHistosPostFit(histos, hsOpt, snames, color, fit_results, postfit_file = None):
    
    print hsOpt
    h_data_bkg = histos["data"].Clone("data-bkg")
    h_data_bkg.Add(histos["bkg"], -1)

    for ibin in range(1, h_data_bkg.GetNbinsX()+1):
        #print h_data_bkg.GetBinError(ibin), histos["bkg"].GetBinError(ibin), histos["bkg"].GetBinError(ibin)/histos["bkg"].GetBinContent(ibin)
        #print h_data_bkg.GetBinError(ibin), histos["bkg"].GetBinError(ibin), histos["data"].GetBinError(ibin), histos["sig"].GetBinError(ibin), histos["sig"].GetBinContent(ibin)
        #print h_data_bkg.GetBinError(ibin), histos["bkg"].GetBinError(ibin), histos["data"].GetBinError(ibin), histos["sig"].GetBinError(ibin)#, histos["bkg"].GetBinError(ibin)/histos["bkg"].GetBinContent(ibin)
        pass

    h_sig = histos["sig"].Clone("signal")
    h_err = histos["sig"].Clone("error_bar")    
    
    #if color: col = color[0][0]
    #else: col = sam_opt[snames[0]]['fillcolor']
    #print col, color, sam_opt[snames[0]]
    
    h_err.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
    #h_err.Reset()
    #herr.Rebin(rebin)
    h_err.GetXaxis().SetTitle(hsOpt['xname'])
    h_err.SetFillStyle(3005)
    h_err.SetFillColor(sam_opt["sig"]['fillcolor'])
    h_err.SetLineColor(922)
    h_err.SetLineWidth(0)         
    h_err.SetMarkerSize(0)
    h_err.SetMarkerColor(922)
    #h_err.SetMinimum(0.)
    
    
    h_sig.SetLineStyle(1)
    h_sig.SetLineWidth(2)
    h_sig.SetLineColor(sam_opt["sig"]['linecolor'])
    
    """for i, h in enumerate(histos):         
        if color: col = color[i]
        else: col = sam_opt[snames[i]]['fillcolor']
        #print sam_opt[snames[i]]['sam_name']
        
        h.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        h.SetMinimum(0.)
        h.Scale(scale)
        h.Rebin(rebin)
        h.SetMarkerStyle(8)
        h.SetMarkerSize(0.)
        print "h[{}]: {}".format(i, h.Integral())
        if fill: 
            h.SetFillColorAlpha(col,0.2)
            h.SetFillStyle(1)
        if i==len(histos)-1 :
            h.SetLineStyle(1)
            h.SetLineWidth(2)
            h.SetLineColor(col)
            if not fill: 
                h.SetMarkerSize(0.6)
                h.SetMarkerColor(col)
        else:     
            h.SetLineWidth(1)         
            h.SetLineColorAlpha(col,0.)
    """
    
    
    
    err = max((fit_results["sig"][0] - fit_results["sig"][1]) / fit_results["sig"][0], (fit_results["sig"][2] - fit_results["sig"][0]) / fit_results["sig"][0])

    #Set error centered at zero as requested by ARC
    for ibin in range(1, h_err.GetNbinsX()+1):
        h_err.SetBinContent(ibin, 0. )

    #If not loading already morphed fit results
    if postfit_file == None:
        for ibin in range(1, h_err.GetNbinsX()+1):
            #print ibin, err, h_err.GetBinContent(ibin), err * h_err.GetBinContent(ibin), h_err.GetBinError(ibin)
            h_err.SetBinError(ibin, math.sqrt((err * h_err.GetBinContent(ibin))**2 + h_data_bkg.GetBinError(ibin)**2) )
    else:
        for ibin in range(1, h_err.GetNbinsX()+1):
            #print ibin, err, h_err.GetBinContent(ibin), err * h_err.GetBinContent(ibin), h_err.GetBinError(ibin)
            h_err.SetBinError(ibin, math.sqrt(h_sig.GetBinError(ibin)**2 + h_data_bkg.GetBinError(ibin)**2) )
            #h_err.SetBinError(ibin, math.sqrt(h_data_bkg.GetBinError(ibin)**2) )
            #print ibin, h_err.GetBinContent(ibin), h_err.GetBinError(ibin)
            
    return h_data_bkg, h_sig, h_err


def getHistosPostFitRatio(histos, hsOpt, snames, color, fit_results, postfit_file = None, only_bias_unc = False):
    
    print hsOpt
    h_data_bkg = histos["data"].Clone("data_over_bkg")
    h_data_bkg.Add(histos["bkg"], -1)
    h_data_bkg.Divide(histos["bkg"])

    for ibin in range(1, h_data_bkg.GetNbinsX()+1):
        #print h_data_bkg.GetBinError(ibin), histos["bkg"].GetBinError(ibin), histos["bkg"].GetBinError(ibin)/histos["bkg"].GetBinContent(ibin)
        #print h_data_bkg.GetBinError(ibin), histos["bkg"].GetBinError(ibin), histos["data"].GetBinError(ibin), histos["sig"].GetBinError(ibin), histos["sig"].GetBinContent(ibin)
        #print h_data_bkg.GetBinError(ibin), histos["bkg"].GetBinError(ibin), histos["data"].GetBinError(ibin), histos["sig"].GetBinError(ibin)#, histos["bkg"].GetBinError(ibin)/histos["bkg"].GetBinContent(ibin)
        pass

    h_sig = histos["sig"].Clone("signal")
    h_sig.Divide(histos["bkg"])
    h_err = histos["sig"].Clone("error_bar")    
    
    #if color: col = color[0][0]
    #else: col = sam_opt[snames[0]]['fillcolor']
    #print col, color, sam_opt[snames[0]]
    
    h_err.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
    #h_err.Reset()
    #herr.Rebin(rebin)
    h_err.GetXaxis().SetTitle(hsOpt['xname'])
    h_err.SetFillStyle(3005)
    h_err.SetFillColor(sam_opt["sig"]['fillcolor'])
    h_err.SetLineColor(922)
    h_err.SetLineWidth(0)         
    h_err.SetMarkerSize(0)
    h_err.SetMarkerColor(922)
    #h_err.SetMinimum(0.)
    
    
    h_sig.SetLineStyle(1)
    h_sig.SetLineWidth(2)
    h_sig.SetLineColor(sam_opt["sig"]['linecolor'])
    
    """for i, h in enumerate(histos):         
        if color: col = color[i]
        else: col = sam_opt[snames[i]]['fillcolor']
        #print sam_opt[snames[i]]['sam_name']
        
        h.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        h.SetMinimum(0.)
        h.Scale(scale)
        h.Rebin(rebin)
        h.SetMarkerStyle(8)
        h.SetMarkerSize(0.)
        print "h[{}]: {}".format(i, h.Integral())
        if fill: 
            h.SetFillColorAlpha(col,0.2)
            h.SetFillStyle(1)
        if i==len(histos)-1 :
            h.SetLineStyle(1)
            h.SetLineWidth(2)
            h.SetLineColor(col)
            if not fill: 
                h.SetMarkerSize(0.6)
                h.SetMarkerColor(col)
        else:     
            h.SetLineWidth(1)         
            h.SetLineColorAlpha(col,0.)
    """
    
    
    
    err = max((fit_results["sig"][0] - fit_results["sig"][1]) / fit_results["sig"][0], (fit_results["sig"][2] - fit_results["sig"][0]) / fit_results["sig"][0])

    #Set error centered at zero as requested by ARC
    for ibin in range(1, h_err.GetNbinsX()+1):
        h_err.SetBinContent(ibin, 0. )

    #If not loading already morphed fit results
    if postfit_file == None:
        for ibin in range(1, h_err.GetNbinsX()+1):
            #print ibin, err, h_err.GetBinContent(ibin), err * h_err.GetBinContent(ibin), h_err.GetBinError(ibin)
            h_err.SetBinError(ibin, math.sqrt((err * h_err.GetBinContent(ibin))**2 + h_data_bkg.GetBinError(ibin)**2) )
    else:
        for ibin in range(1, h_err.GetNbinsX()+1):
            #print ibin, err, h_err.GetBinContent(ibin), err * h_err.GetBinContent(ibin), h_err.GetBinError(ibin)
            if not only_bias_unc:
                h_err.SetBinError(ibin, math.sqrt(h_sig.GetBinError(ibin)**2 + h_data_bkg.GetBinError(ibin)**2) )
            else:
                h_err.SetBinError(ibin, math.sqrt(h_data_bkg.GetBinError(ibin)**2) )
            #print ibin, h_err.GetBinContent(ibin), h_err.GetBinError(ibin)
            
    return h_data_bkg, h_sig, h_err



"""def make_residual_comparison(hlist1, hsOpt, rb, snames1, colors[0], scale1, dofill[0], hlist2, hsOpt, rb, snames2, colors[1], scale2, dofill[1])
    
    hs1, herr1, h1 =  getStackH(hlist1, hsOpt, rb, snames1, colors[0], scale1, dofill[0])
    if hs1.GetMaximum() > ymax: ymax = hs1.GetMaximum()*1.15
    if herr1.GetMaximum() > ymax: ymax = herr1.GetMaximum()*1.15
    if isNevts:  print "h1Int ",  h1.GetBinContent(1), h1.GetBinError(1)

    if(norm): scale2 = getScale(hsOpt,hlist2, hlist2)
    else: scale2 = 1.
    if scale2 != 1.: print "sc_to_norm2: ",scale2
    hs2, herr2, h2 =  getStackH(hlist2, hsOpt, rb, snames2, colors[1], scale2, dofill[1])"""

def get_odir(args, oname, option=""):
    oDir = args.oDir
    if not os.path.exists(oDir): os.mkdir(oDir)
    oDir += "/"+args.report_dir
    if not os.path.exists(oDir): os.mkdir(oDir)
    
    #bm = args.bdt.split("-")[2]
    #oDir += "/"+bm
    
    if not os.path.exists(oDir): os.mkdir(oDir)
    oDir += "/"+oname
    if args.doNorm: oDir = oDir+"_norm"
    if args.clrebin > 1: 
        oDir += "_rebin_" + str(args.clrebin)
    oDir = oDir+"/"
    if not os.path.exists(oDir): os.mkdir(oDir)
    oDir += option #keep the second sample options
    if not os.path.exists(oDir): os.mkdir(oDir)
    return oDir
    

def tdrGrid( gridOn):
  tdrStyle.SetPadGridX(gridOn)
  tdrStyle.SetPadGridY(gridOn)

#fixOverlay: Redraws the axis
def fixOverlay(): gPad.RedrawAxis()
    
def setTDRStyle():
  tdrStyle =  ROOT.TStyle("tdrStyle","Style for P-TDR")

   #for the canvas:
  tdrStyle.SetCanvasBorderMode(0)
  tdrStyle.SetCanvasColor(ROOT.kWhite)
  tdrStyle.SetCanvasDefH(800) #Height of canvas
  tdrStyle.SetCanvasDefW(800) #Width of canvas
  tdrStyle.SetCanvasDefX(0)   #POsition on screen
  tdrStyle.SetCanvasDefY(0)


  tdrStyle.SetPadBorderMode(0)
  #tdrStyle.SetPadBorderSize(Width_t size = 1)
  tdrStyle.SetPadColor(ROOT.kWhite)
  tdrStyle.SetPadGridX(False)
  tdrStyle.SetPadGridY(False)
  tdrStyle.SetGridColor(0)
  tdrStyle.SetGridStyle(3)
  tdrStyle.SetGridWidth(1)

#For the frame:
  tdrStyle.SetFrameBorderMode(0)
  tdrStyle.SetFrameBorderSize(1)
  tdrStyle.SetFrameFillColor(0)
  tdrStyle.SetFrameFillStyle(0)
  tdrStyle.SetFrameLineColor(1)
  tdrStyle.SetFrameLineStyle(1)
  tdrStyle.SetFrameLineWidth(1)
  
#For the histo:
  #tdrStyle.SetHistFillColor(1)
  tdrStyle.SetHistFillStyle(1001)
  tdrStyle.SetHistLineColor(1)
  tdrStyle.SetHistLineStyle(0)
  tdrStyle.SetHistLineWidth(1)
  #tdrStyle.SetLegoInnerR(Float_t rad = 0.5)
  #tdrStyle.SetNumberContours(Int_t number = 20)

  tdrStyle.SetEndErrorSize(2)
  #tdrStyle.SetErrorMarker(20)
  #tdrStyle.SetErrorX(0.)
  
  tdrStyle.SetMarkerStyle(20)
  
#For the fit/function:
  tdrStyle.SetOptFit(1)
  tdrStyle.SetFitFormat("5.4g")
  tdrStyle.SetFuncColor(2)
  tdrStyle.SetFuncStyle(1)
  tdrStyle.SetFuncWidth(1)

#For the date:
  tdrStyle.SetOptDate(0)
  # tdrStyle.SetDateX(Float_t x = 0.01)
  # tdrStyle.SetDateY(Float_t y = 0.01)

# For the statistics box:
  tdrStyle.SetOptFile(0)
  tdrStyle.SetOptStat(0) # To display the mean and RMS:   SetOptStat("mr")
  tdrStyle.SetStatColor(ROOT.kWhite)
  tdrStyle.SetStatFont(42)
  tdrStyle.SetStatFontSize(0.025)
  tdrStyle.SetStatTextColor(1)
  tdrStyle.SetStatFormat("6.4g")
  tdrStyle.SetStatBorderSize(1)
  tdrStyle.SetStatH(0.1)
  tdrStyle.SetStatW(0.15)
  # tdrStyle.SetStatStyle(Style_t style = 1001)
  # tdrStyle.SetStatX(Float_t x = 0)
  # tdrStyle.SetStatY(Float_t y = 0)

# Margins:
  tdrStyle.SetPadTopMargin(0.05)
  tdrStyle.SetPadBottomMargin(0.13)
  tdrStyle.SetPadLeftMargin(0.16)
  tdrStyle.SetPadRightMargin(0.02)

# For the Global title:

  tdrStyle.SetOptTitle(0)
  tdrStyle.SetTitleFont(42)
  tdrStyle.SetTitleColor(1)
  tdrStyle.SetTitleTextColor(1)
  tdrStyle.SetTitleFillColor(10)
  tdrStyle.SetTitleFontSize(0.05)
  # tdrStyle.SetTitleH(0) # Set the height of the title box
  # tdrStyle.SetTitleW(0) # Set the width of the title box
  # tdrStyle.SetTitleX(0) # Set the position of the title box
  # tdrStyle.SetTitleY(0.985) # Set the position of the title box
  # tdrStyle.SetTitleStyle(Style_t style = 1001)
  # tdrStyle.SetTitleBorderSize(2)

# For the axis titles:

  tdrStyle.SetTitleColor(1, "XYZ")
  tdrStyle.SetTitleFont(42, "XYZ")
  tdrStyle.SetTitleSize(0.06, "XYZ")
  # tdrStyle.SetTitleXSize(Float_t size = 0.02) # Another way to set the size?
  # tdrStyle.SetTitleYSize(Float_t size = 0.02)
  tdrStyle.SetTitleXOffset(0.9)
  tdrStyle.SetTitleYOffset(1.25)
  # tdrStyle.SetTitleOffset(1.1, "Y") # Another way to set the Offset

# For the axis labels:

  tdrStyle.SetLabelColor(1, "XYZ")
  tdrStyle.SetLabelFont(42, "XYZ")
  tdrStyle.SetLabelOffset(0.007, "XYZ")
  tdrStyle.SetLabelSize(0.05, "XYZ")

# For the axis:

  tdrStyle.SetAxisColor(1, "XYZ")
  tdrStyle.SetStripDecimals(True)
  tdrStyle.SetTickLength(0.03, "XYZ")
  tdrStyle.SetNdivisions(510, "XYZ")
  tdrStyle.SetPadTickX(1)  # To get tick marks on the opposite side of the frame
  tdrStyle.SetPadTickY(1)

# Change for log plots:
  #tdrStyle.SetOptLogx(0)
  #tdrStyle.SetOptLogy(0)
  #tdrStyle.SetOptLogz(0)

# Postscript options:
  #tdrStyle.SetPaperSize(20.,20.)
  # tdrStyle.SetLineScalePS(Float_t scale = 3)
  # tdrStyle.SetLineStyleString(Int_t i, const char* text)
  # tdrStyle.SetHeaderPS(const char* header)
  # tdrStyle.SetTitlePS(const char* pstitle)

  # tdrStyle.SetBarOffset(Float_t baroff = 0.5)
  # tdrStyle.SetBarWidth(Float_t barwidth = 0.5)
  # tdrStyle.SetPaintTextFormat(const char* format = "g")
  # tdrStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0)
  # tdrStyle.SetTimeOffset(Double_t toffset)
  # tdrStyle.SetHistMinimumZero(kTRUE)

  tdrStyle.SetHatchesLineWidth(5)
  tdrStyle.SetHatchesSpacing(0.05)

  tdrStyle.cd()
  
  
# CMS_lumi
#   Initiated by: Gautier Hamel de Monchenault (Saclay)
#   Translated in Python by: Joshua Hardenbrook (Princeton)
#   Updated by:   Dinko Ferencek (Rutgers)
#


def CMS_lumi(pad,  iPeriod,  iPosX ):
    cmsText     = "CMS";
    cmsTextFont   = 61  

    writeExtraText = True
    extraText   = "Preliminary"
    extraTextFont = 52 

    lumiTextSize     = 0.6
    lumiTextOffset   = 0.2

    cmsTextSize      = 0.75
    cmsTextOffset    = 0.1

    relPosX    = 0.045
    relPosY    = 0.035
    relExtraDY = 1.2

    extraOverCmsTextSize  = 0.76

    lumi_13TeV = "35.9 fb^{-1}"
    lumi_8TeV  = "19.7 fb^{-1}" 
    lumi_7TeV  = "5.1 fb^{-1}"
    lumi_sqrtS = ""

    drawLogo      = False

    outOfFrame    = False
    if(iPosX/10==0 ): outOfFrame = True

    alignY_=3
    alignX_=2
    if( iPosX/10==0 ): alignX_=1
    if( iPosX==0    ): alignY_=1
    if( iPosX/10==1 ): alignX_=1
    if( iPosX/10==2 ): alignX_=2
    if( iPosX/10==3 ): alignX_=3
    align_ = 10*alignX_ + alignY_

    H = pad.GetWh()
    W = pad.GetWw()
    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()
    e = 0.025

    pad.cd()

    lumiText = ""
    if ( iPeriod==4 ):
        lumiText += lumi_13TeV
        lumiText += " (13 TeV)"
    elif ( iPeriod==7 ):
        if( outOfFrame ):lumiText += "#scale[0.85]{"
        lumiText += lumi_13TeV 
        lumiText += " (13 TeV)"
        lumiText += " + "
        lumiText += lumi_8TeV 
        lumiText += " (8 TeV)"
        lumiText += " + "
        lumiText += lumi_7TeV
        lumiText += " (7 TeV)"
        if( outOfFrame): lumiText += "}"
    elif ( iPeriod==12 ):
        lumiText += "8 TeV"
    elif ( iPeriod==0 ):
        lumiText += lumi_sqrtS
            
    print lumiText

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)    
    
    extraTextSize = extraOverCmsTextSize*cmsTextSize
    
    latex.SetTextFont(42)
    latex.SetTextAlign(31) 
    latex.SetTextSize(lumiTextSize*t)    

    latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText)

    if( outOfFrame ):
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11) 
        latex.SetTextSize(cmsTextSize*t)    
        latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText)
  
    pad.cd()

    posX_ = 0
    if( iPosX%10<=1 ):
        posX_ =   l + relPosX*(1-l-r)
    elif( iPosX%10==2 ):
        posX_ =  l + 0.5*(1-l-r)
    elif( iPosX%10==3 ):
        posX_ =  1-r - relPosX*(1-l-r)

    posY_ = 1-t - relPosY*(1-t-b)

    if( not outOfFrame ):
        if( drawLogo ):
            posX_ =   l + 0.045*(1-l-r)*W/H
            posY_ = 1-t - 0.045*(1-t-b)
            xl_0 = posX_
            yl_0 = posY_ - 0.15
            xl_1 = posX_ + 0.15*H/W
            yl_1 = posY_
            CMS_logo = ROOT.TASImage("CMS-BW-label.png")
            pad_logo =  ROOT.TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 )
            pad_logo.Draw()
            pad_logo.cd()
            CMS_logo.Draw("X")
            pad_logo.Modified()
            pad.cd()          
        else:
            latex.SetTextFont(cmsTextFont)
            latex.SetTextSize(cmsTextSize*t)
            latex.SetTextAlign(align_)
            latex.DrawLatex(posX_, posY_, cmsText)
            if( writeExtraText ) :
                latex.SetTextFont(extraTextFont)
                latex.SetTextAlign(align_)
                latex.SetTextSize(extraTextSize*t)
                latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText)
    elif( writeExtraText ):
        if( iPosX==0):
            posX_ =   l +  relPosX*(1-l-r)
            posY_ =   1-t+lumiTextOffset*t

        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(extraTextSize*t)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, extraText)      

    pad.Update()
    
    
def get_bias_corrected_histo(histo, region = "", only_bias_unc = False):
    if region == "ms":
        bkg_bias_fname = "/lustre/cmswork/dcastrom/projects/hh/april_2017/CMSSW_8_0_25/src/Analysis/hh2bbbb_limit/notebooks/bias_22032018_with_weights_also_mass_cut/BM0/bias_correction_mass_cut_bigset_unscaled.json"
        bkg_bias_fname = "/lustre/cmswork/dcastrom/projects/hh/april_2017/CMSSW_8_0_25/src/Analysis/hh2bbbb_limit/bias_08062018BM0/bias_correction_mass_cut_bigset_unscaled.json"
    elif region == "btag":
        bkg_bias_fname = "/lustre/cmswork/dcastrom/projects/hh/april_2017/CMSSW_8_0_25/src/Analysis/hh2bbbb_limit/notebooks/bms_btagsideBM0/bias_correction_bigset_unscaled.json"
    elif region == "sig_unfixed":
        bkg_bias_fname = "/lustre/cmswork/dcastrom/projects/hh/april_2017/CMSSW_8_0_25/src/Analysis/hh2bbbb_limit/notebooks/bias_22032018_with_weights_also_mass_cut/BM0/bias_correction_bigset_unscaled.json"
    else:
        bkg_bias_fname = "/lustre/cmswork/dcastrom/projects/hh/april_2017/CMSSW_8_0_25/src/Analysis/hh2bbbb_limit/bias_08062018BM0/bias_correction_bigset_unscaled.json"
        
    with open(bkg_bias_fname,"r") as bkg_bias_file:
        json_dict = json.load(bkg_bias_file)
        print ("using bias file: ", bkg_bias_file)

    hcorr = histo.Clone("h_bias_corrected")
    hcorr.Scale(4)
    hbias = histo.Clone("h_bias")
    for n in range(len(json_dict['bias_corr'])):
        #if region == "btag" and n == len(json_dict['bias_corr']) - 1: continue
        #if not s_bin.overflow:
        #value = s_bin.value

        bias = json_dict['bias_corr'][n]
        var = json_dict['var'][n]
        bias_unc = json_dict['bias_corr_unc_bs'][n]
        bias_unc_stat = json_dict['bias_corr_unc_stat'][n]
        
        bkg_pred_initial = hcorr.GetBinContent(n+1)
        if var > np.sqrt(bkg_pred_initial):
          new_bkg_pred_stat = var
        else:
          new_bkg_pred_stat = np.sqrt(bkg_pred_initial)

                        
        new_bkg_pred_tot_unc = np.sqrt(new_bkg_pred_stat**2 + bias_unc**2 + bias_unc_stat**2) * 1.93
        
        if only_bias_unc == True:
            new_bkg_pred_tot_unc = np.sqrt(bias_unc**2 + bias_unc_stat**2) * 1.93
        
        #hbias.SetBinContent(n+1, bias)
        #hbias.SetBinError(n+1, new_bkg_pred_tot_unc)
        print "Corr", n, bias, hcorr.GetBinContent(n+1), var, bkg_pred_initial, new_bkg_pred_tot_unc, new_bkg_pred_stat, bias_unc, bias_unc_stat
        hcorr.SetBinContent(n+1, hcorr.GetBinContent(n+1) - bias)
        hcorr.SetBinError(n+1, new_bkg_pred_tot_unc)
        print "Corr", n, bias, hcorr.GetBinContent(n+1)#, new_bkg_pred_tot_unc
        
    hcorr.Scale(0.25)
    #hcorr.Add(hbias, -1)
    return hcorr
    
    
    
def getHistosPostFitTemp(histos, hsOpt, snames, color, fit_results, postfit_file = None):
    
    print hsOpt
    h_bkg_bias_corr = get_bias_corrected_histo(histos["bkg"], region = "", only_bias_unc = True)
    h_data_bkg = histos["data"].Clone("data-bkg")
    h_data_bkg.Add(h_bkg_bias_corr, -1)

    print histos
    for ibin in range(1, h_data_bkg.GetNbinsX()+1):
        #print "asi", ibin, h_bkg_bias_corr.GetBinContent(ibin), histos["bkg"].GetBinContent(ibin), histos["data"].GetBinContent(ibin)
        #h_data_bkg.GetBinError(ibin), histos["bkg"].GetBinError(ibin), histos["bkg"].GetBinError(ibin)/histos["bkg"].GetBinContent(ibin)
        #print h_data_bkg.GetBinError(ibin), histos["bkg"].GetBinError(ibin), histos["bkg"].GetBinError(ibin)/histos["bkg"].GetBinContent(ibin)
        #print h_data_bkg.GetBinError(ibin), histos["bkg"].GetBinError(ibin), histos["data"].GetBinError(ibin), histos["sig"].GetBinError(ibin), histos["sig"].GetBinContent(ibin)
        #print h_data_bkg.GetBinError(ibin), histos["bkg"].GetBinError(ibin), histos["data"].GetBinError(ibin), histos["sig"].GetBinError(ibin)#, histos["bkg"].GetBinError(ibin)/histos["bkg"].GetBinContent(ibin)
        pass

    h_sig = histos["sig"].Clone("signal")
    h_err = histos["sig"].Clone("error_bar")    
    
    #if color: col = color[0][0]
    #else: col = sam_opt[snames[0]]['fillcolor']
    #print col, color, sam_opt[snames[0]]
    
    h_err.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
    #h_err.Reset()
    #herr.Rebin(rebin)
    h_err.GetXaxis().SetTitle(hsOpt['xname'])
    h_err.SetFillStyle(3005)
    h_err.SetFillColor(sam_opt["sig"]['fillcolor'])
    h_err.SetLineColor(922)
    h_err.SetLineWidth(0)         
    h_err.SetMarkerSize(0)
    h_err.SetMarkerColor(922)
    #h_err.SetMinimum(0.)
    
    
    h_sig.SetLineStyle(1)
    h_sig.SetLineWidth(2)
    h_sig.SetLineColor(sam_opt["sig"]['linecolor'])
    
    #Set error centered at zero as requested by ARC
    for ibin in range(1, h_err.GetNbinsX()+1):
        h_err.SetBinContent(ibin, 0. )

    #If not loading already morphed fit results
    for ibin in range(1, h_err.GetNbinsX()+1):
        #print ibin, err, h_err.GetBinContent(ibin), err * h_err.GetBinContent(ibin), h_err.GetBinError(ibin)
        h_err.SetBinError(ibin, h_data_bkg.GetBinError(ibin))
            
    return h_data_bkg, h_sig, h_err

