import json
import os
import math
from glob import glob

# ROOT imports
from ROOT import TChain, TPad, TH1D, TH2D, TFile, vector, TCanvas, TLatex, TLine, TLegend, THStack, gStyle, TGaxis
from rootpy.plotting import Hist, Hist2D

from Analysis.alp_analysis.alpSamplesOptions  import sam_opt

##
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

def getHistos_bdt(hist, filename, plotDirs, weights, sf):
    hlist = [] 
    tf = TFile(filename)
    if not tf: 
        print "WARNING: files do not exist"  
    #print 'aa'
    for i, Dir in enumerate(plotDirs):
        hname = hist+"_"+Dir            
        #print hname
        if(tf.Get(Dir+"/"+hname)):
            w = 1.
            sf_ = 1.
            if len(weights)>0:
	            if weights[i]>0:
        	        w = weights[i]
            if len(sf)>0:
                if sf[i]>=0:
                    sf_ = sf[i]
            print w*sf_
            h = tf.Get(Dir+"/"+hname)
            #print h.GetName()
            #print h.Integral()
            h.Scale(w*sf_)
            print h.Integral()
            hlist.append(h)
        else:
            print "WARNING: hist {} not found in {}/{}".format(hist,tf,Dir)

    return hlist
#------------

def getHistos_tdr(hist, filelist, plotDir, lumi, normtolumi, weight):
    hlist = []    
    for i, f in enumerate(filelist):  #debug - not efficient to loop on file
        w = 1.
        tf = TFile(f)
        #print tf
        if not tf: 
            print "WARNING: files do not exist"  

        if weight[i]>=0:
   	    w = weight[i]

        #print w
        if(tf.Get(plotDir+"/"+hist)):
            h = tf.Get(plotDir+"/"+hist)
            # print plotDir
            #print "h1Int ",  h.GetBinContent(1), h.GetBinError(1)
            h.Scale(w)
            hlist.append(h)
            print h.Integral()
        else:
            print "WARNING: hist {} not found in {}".format(hist,tf)

    return hlist
#------------


def getHistos(hist, filelist, plotDir, lumi, normtolumi, weight, sf):
    hlist = []    
    for i, f in enumerate(filelist):  #debug - not efficient to loop on file
        norm = 1
        w = 1.
        sf_ = 1.
        tf = TFile(f)
       # print f
        if not tf: 
            print "WARNING: files do not exist"  

        if "Run" in f:   #data
            if not "mixed" in f:
                print "is Data"
        else:
            if len(weight)>0:
                if weight[i]>=0:
                    w = weight[i]
            else:
                if not "mixed" in f: #only for MC
                    if(tf.Get("h_w_oneInvFb")):
                        h = tf.Get("h_w_oneInvFb")
                        w = h.GetBinContent(1)
                    else:
                        print "WARNING: 'h_w_oneInvFb' not found in {}".format(tf)
            if len(sf)>0:
                if sf[i]>=0:
                    sf_ = sf[i]
            if normtolumi: norm = w*lumi*sf_
            else: norm = w*sf_

        print norm    
        if(tf.Get(plotDir+"/"+hist)):
            h = tf.Get(plotDir+"/"+hist)
            #print "h1Int ",  h.GetBinContent(1), h.GetBinError(1)
            h.Scale(norm)
            hlist.append(h)
        else:
            print "WARNING: hist {} not found in {}/{}".format(hist,tf,plotdir)

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
        print h.Integral()
    for i, h in enumerate(hlist2): 
        if isinstance(h,int):
           norm = h
           break
        elif i==skip2: 
            continue 
        else: norm += h.Integral()
    if  hsInt: scale = norm/(hsInt+hsSkip)
    else: scale = 1.
    print hsInt
    return scale
#------------

def getStackH(histos, hsOpt, rebin, snames, color, scale, fill, hl2):

    #if color: col = color
    #else: col = sam_opt[snames[0]]['fillcolor']

    herr = histos[0].Clone("hs_error")  
    herr.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
    herr.Reset()
    herr.Rebin(rebin)
    herr.GetXaxis().SetTitle(hsOpt['xname'])
    herr.GetYaxis().SetTitle('Events') #hsOpt['yname']
    #herr.GetYaxis().SetTitleSize(30)
    #herr.GetYaxis().SetTitleFont(43)
    #herr.GetYaxis().SetTitleOffset(1.40)
    #herr.GetYaxis().SetLabelFont(43)
    #herr.GetYaxis().SetLabelSize(25)
    herr.SetFillStyle(3005)
   # herr.SetFillColor(col)
    herr.SetLineColor(922)
    herr.SetLineWidth(0)         
    herr.SetMarkerSize(0)
    herr.SetMarkerColor(922)
    herr.SetMinimum(0.)

    hs   = THStack("hs","")
    for i, h in enumerate(histos):         
        if color: col = color[i]
        else: col = sam_opt[snames[i]]['fillcolor']
        #print sam_opt[snames[i]]['sam_name']

        h.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        h.SetMinimum(0.)
        h.Scale(scale)
        h.Rebin(rebin)
        h.SetMarkerStyle(8)
        h.SetMarkerSize(0.)
        print h.Integral()
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
        if i==0: h_  = h.Clone("h_")
        else: h_.Add(h)
        hs.Add(h)
        herr.Add(h)
    return hs, herr, h_

##
# DRAWING FUNCTIONS
#################
def drawEffH(histos, leg, fld, snames1, oDir, colors):
    gStyle.SetOptStat(False)
    ceff = TCanvas("ceff", "eff", 800, 800)

    for i, h in enumerate(histos):
        if i == 0: continue
        h.SetMarkerStyle(8)
        h.SetMarkerSize(0.9)
        #h.GetXaxis().SetRangeUser(0,1)
        h.SetMinimum(0.)
        h.GetXaxis().SetTitle(leg)
        h.GetYaxis().SetTitle("acc x eff")
        h.SetLineColor(colors[i])
        h.SetLineStyle(1)
        h.SetLineWidth(2)
        if i==1 : h.Draw("HIST")
        else: h.Draw("HISTsame")

    ceff.Update()
    ceff.SaveAs(oDir+"/"+"eff.pdf")
    ceff.SaveAs(oDir+"/"+"eff.png")
    ceff.SaveAs(oDir+"/"+leg+".root")
#------------

def drawChiSquare(hlist, snames, legstack, h2, hsOpt, oDir, xbmin, headerOpt, isMC, labels):  
    gStyle.SetOptStat(False)
    c1 = TCanvas("c1", hsOpt['hname'], 800, 800)       

    x_data =[]
    chisq_ =[]
    n_sample = 0
    nbins = hlist[0].GetNbinsX()
    print nbins

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
            print b, x
            print b, x_data[b-xbmin] 
            if (x+x_data[b-xbmin]) > 0: chisq_[i] += pow((x_data[b-xbmin]-x),2)/(x+x_data[b-xbmin])
            else: 
#                chisq_[i] +=-99.
                print 'sample data bin {} : null bin content?'.format(b)       

    print "\n", n_sample  

    hchi = (TH1D)("h_chisq", "chi square", n_sample, 0., n_sample) 
    hchi.SetMarkerStyle(8)
    hchi.SetMarkerSize(1.)
    hchi.SetMarkerColor(1)

    for i, h in enumerate(hlist):         
        hchi.SetBinContent(i+1, chisq_[i])
#        hchi.SetBinError(j+1, dDiff)
        print hchi.GetBinContent(i+1) #, hvar.GetBinError(j+1)       

    x_axis = hchi.GetXaxis()
    for b_n, label in enumerate(labels): 
       # print b_n+1, label
        x_axis.SetBinLabel(b_n+1, label)        

    hchi.GetYaxis().SetTitle("chi2 data-mixed")
    hchi.Draw("PE")
    drawCMS(35.9, headerOpt)

    c1.Update()    
    c1.SaveAs(oDir+"/"+"bdtChisq_allnn.pdf") #allnn   5nn
    c1.SaveAs(oDir+"/"+"bdtChisq_allnn.png")  #allnn    5nn    
    c1.SaveAs(oDir+"/"+"bdtChisq_allnn.root")  #allnn  5nn      

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
    print nbins
    for k in range(0, nbins):
        x_.append(0.)
        xsq_.append(0.)
        sum_err.append(0.)

    for i, h in enumerate(hlist):         
        n_sample +=1
        #h.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        #h.SetMinimum(0.)
        for b in range(1, nbins+1):
            x = h.GetBinContent(b)
            if b==nbins: print b, x
            x_[b-1] += x
            xsq_[b-1] += pow(x,2)
            #print b, xs_mean[b-1], xs_sqmean[b-1]

    print "\n", n_sample

    for k in range(0, nbins):
        xs_mean.append(x_[k]/n_sample)
        xs_sqmean.append(xsq_[k]/n_sample)
        #print k, xs_mean[k], xs_sqmean[k]

    for i, h in enumerate(hlist):         
        for b in range(1, nbins+1):
            x = h.GetBinContent(b)
#            sum_err[b-1] +=  (pow(x,2) + pow(xs_mean[b-1],2) + 2*x*xs_mean[b-1])*x            
            sum_err[b-1] += pow((2*x - 2*xs_mean[b-1] -1),2)*x

    for j, m in enumerate(xs_mean):
        var = (xs_sqmean[j] - pow(m,2))
        #dVar = math.sqrt(4./pow(n_sample,2)*sum_err[j])
        #dVar = 2*pow(var,2)/(n_sample-1)
        #print var #, dVar
        #print m, math.sqrt(m)

       # if x>0. : 
           # r = var / m
           # dRat = math.sqrt(pow(dVar/var,2)+pow(math.sqrt(m)/m,2))*r
        diff = var - m        
        dDiff = math.sqrt(1./pow(n_sample,2)*sum_err[j])
       # else: 
        #    print 'null mean - r set to zero'
         #   diff = 0.
          #  dDiff = 0.

        hvar.SetBinContent(j+1, diff)
        hvar.SetBinError(j+1, dDiff)
        print hvar.GetBinContent(j+1), hvar.GetBinError(j+1)       

#    h.SetMarkerSize(0.6)
#    h.SetMarkerColor(col)
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
    c1.SaveAs(oDir+"/"+"bdtVar_nn.pdf") #allnn   5nn
    c1.SaveAs(oDir+"/"+"bdtVar_nn.png")  #allnn    5nn    
    c1.SaveAs(oDir+"/"+"bdtVar_nn.root")  #allnn  5nn      

#------------

def drawH1(hlist1, snames1, legstack1, hlist2, snames2, legstack2, hsOpt, residuals, norm, oDir, colors, dofill, rebin, headerOpt, isMC):  
    gStyle.SetOptStat(False)
    gStyle.SetOptTitle(0);

    c1 = TCanvas("c1", hsOpt['hname'], 800, 800)       
   # c1.cd()
    if residuals==-1:
        pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0.03) 
        pad1.Draw()             
        pad1.cd()               

    if rebin > 0: rb = rebin
    else: rb =  hsOpt['rebin']

    isNevts=False
    if hsOpt['hname']=="h_nevts": isNevts=True
    ymax = 0.

    if(norm): scale1 = getScale(hsOpt,hlist1, hlist2)
    else: scale1 = 1.
    print "scale1", scale1
    hs1, herr1, h1 =  getStackH(hlist1, hsOpt, rb, snames1, colors[0], scale1, dofill[0], hlist2)
    if hs1.GetMaximum() > ymax: ymax = hs1.GetMaximum()*1.15
    if herr1.GetMaximum() > ymax: ymax = herr1.GetMaximum()*1.15
    if isNevts:  print "h1Int ",  h1.GetBinContent(1), h1.GetBinError(1)
#    if hsOpt['hname']=="h_jets_n":  print "h1Int ",  h1.Integral(), h1.GetBinError(1)

    if(norm): scale2 = getScale(hsOpt,hlist2, hlist2)
    else: scale2 = 1.
    print "scale2", scale2
    hs2, herr2, h2 =  getStackH(hlist2, hsOpt, rb, snames2, colors[1], scale2, dofill[1], hlist2)
    if hs2.GetMaximum() > ymax: ymax = hs2.GetMaximum()*1.15
    if herr2.GetMaximum() > ymax: ymax = herr2.GetMaximum()*1.15
    if isNevts:  print "h2Int ",  h2.GetBinContent(1), h2.GetBinError(1)
#    if hsOpt['hname']=="h_jets_n":  print "h2Int ",  h2.Integral(), h2.GetBinError(1)

    #print herr1.GetBinContent(1), herr1.GetBinError(1)

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
    if len(legstack1) > 1:
        #for n, leg in enumerate(legstack1):
        #    legend.AddEntry(hlist1[n], leg)
        legend.AddEntry(hlist1[0], legstack1[0])
        legend.AddEntry(hlist1[len(hlist1)-2], legstack1[1]) #trg
        legend.AddEntry(hlist1[len(hlist1)-1], legstack1[2]) #trg
        #legend.AddEntry(hlist1[len(hlist1)-1], legstack1[1]) #debug
    else:
        legend.AddEntry(hlist1[len(hlist1)-1], legstack1[0])
    if len(hlist1)>1: legend.AddEntry(herr1, 'bkg. unc. (stat.only)')
    legend.Draw("same")
    #-------------

    if(ymax > 1000): TGaxis.SetMaxDigits(3)
    hs1.SetMaximum(ymax)
    herr1.SetMaximum(ymax)
    hs1.Draw("HISTsame") #Esame
    if len(hlist1)>1: herr1.Draw("E2same")
    else: herr1.Draw("Esame")
    hs2.SetMaximum(ymax)
    herr2.SetMaximum(ymax)
    if dofill[1]: hs2.Draw("HISTsame")
    else: hs2.Draw("Esame")
    #if len(hlist2)>1: herr2.Draw("E2same")
    herr2.Draw("Esameaxis")
#h1->Draw("sameaxis")

    nev1 = 0
    nev2 = 0
    nev1err = 0
    nev2err = 0
    if isNevts: 
#    if hsOpt['hname']=='h_jets_n': 
        nev1 = herr1.GetBinContent(1)
        nev2 = herr2.GetBinContent(1)
        nev1err = herr1.GetBinError(1)
        nev2err = herr2.GetBinError(1)
    #nev1 = herr1.Integral()
    #nev2 = herr2.Integral()
    #nev1err = herr1.Integral()
    #nev2err = herr2.Integral()

    if not residuals==-1:
        c1.Update()    
        c1.SaveAs(oDir+"/"+hsOpt['hname']+".pdf")
        c1.SaveAs(oDir+"/"+hsOpt['hname']+".png")            
        c1.SaveAs(oDir+"/"+hsOpt['hname']+".root") 
    else:
        herr1.GetXaxis().SetLabelSize(0.)

    if residuals==-1: # simple ratio -- data as h2!!!
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

        #for i in range(1, hrat.GetXaxis().GetNbins()+1):
        #    n1 = h1.GetBinContent(i)
        #    n2 = h2.GetBinContent(i)
        #    if n1 and n2 and e1 and e2 : 
        #        hrat.SetBinContent(i,n2/n1)
        # to chek just histos ratio:

        hrat.SetTitle("")
        hrat.GetXaxis().SetTitleSize(20)
        hrat.GetXaxis().SetTitleFont(43)
        hrat.GetXaxis().SetTitleOffset(4.)
        hrat.GetXaxis().SetLabelFont(43)
        hrat.GetXaxis().SetLabelSize(20)
        hrat.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        ymax_ = 1.9 
        ymin_ = 0.1
        ymax = 1.
        ymin = 1.
        for ibin in range(1, hrat.GetNbinsX()+1):       
            binc = hrat.GetBinContent(ibin) 
            if ymax_ > binc > ymax: ymax = (binc + 0.05)
            if ymin_ < binc < ymin: ymin = (binc - 0.05)
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
        hrat.GetYaxis().SetTitle('Data/MC')
        hrat.Draw("E X0")
        herr.SetFillColor(430)
        herr.Draw("E2same")

        #l0 = TLine(hsOpt['xmin'],1.4,hsOpt['xmax'],1.4);
        l1 = TLine(hsOpt['xmin'],1.1,hsOpt['xmax'],1.1);
        l2 = TLine(hsOpt['xmin'],1.,hsOpt['xmax'],1.);
        l3 = TLine(hsOpt['xmin'],0.9,hsOpt['xmax'],0.9);
        #l4 = TLine(hsOpt['xmin'],0.6,hsOpt['xmax'],0.6);
        #l0.SetLineStyle(3)
        l1.SetLineStyle(3)
        l2.SetLineStyle(3)
        l3.SetLineStyle(3)
        #l4.SetLineStyle(3)
        #l0.Draw("same")
        l1.Draw("same")
        l2.Draw("same")
        l3.Draw("same")
        #l4.Draw("same")

        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.pdf")
        c1.SaveAs(oDir+"/"+hsOpt['hname']+"_rat.png")            
        c1.Clear()

    elif residuals >=1: # data as h2!!!
        if residuals==2:
            c2 = TCanvas("c2", "res"+hsOpt['hname'], 800, 800)    
            c2.Divide(1,2)
            c2.cd(1)
        else:
            c2 = TCanvas("c2", "res"+hsOpt['hname'], 800, 400)
        hres = h2.Clone("h_res")
      #  hres = h2.Reset()
        checkbin = False
        for i in range(1, hres.GetXaxis().GetNbins()+1):
            n1 = h1.GetBinContent(i)
            n2 = h2.GetBinContent(i)
            e1 = h1.GetBinError(i)
            e2 = h2.GetBinError(i)
            print  n1, n2, e1, e2
            if n1 and n2 and e1 and e2 : 
                hres.SetBinContent(i,(n2-n1)/math.sqrt(e1*e1+e2*e2)) #debug v1
            #elif not checkbin: hres.SetBinContent(i,0) #to avoid error
               # hres.SetBinError(i,1)
                err = (pow(n1,3) + 15*pow(n1,2)*n2+15*pow(n2,2)*n1 + pow(n2,3))/(4*pow((n1+n2),3))
                hres.SetBinError(i, err)
     # to chek just histos ratio:
     #  hres = Hist(10,0,1, name="")
     #  hres = hdiff.Clone("h_residual")
     #  hres.Divide(h1)

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
            for i in range(0, hres.GetXaxis().GetNbins()): #improve
                if (not hres.IsBinOverflow(i)) and (math.fabs(hres.GetBinContent(i)) != 0):
                    res_a.append(hres.GetBinContent(i))
            res_pull = Hist(60, -6., 6, name = "h_res_pull")
            res_pull.GetXaxis().SetTitle("residuals (sigma units)")
            for v in res_a: res_pull.fill(v)
            res_pull.Fit("gaus", "ILL")
            res_pull.Draw("E X0")

        c2.Update()    
        c2.SaveAs(oDir+"/"+hsOpt['hname']+"_res.pdf")
        c2.SaveAs(oDir+"/"+hsOpt['hname']+"_res.png")            
        c2.SaveAs(oDir+"/"+hsOpt['hname']+"_res.root")   

    return [nev1,nev1err,nev2,nev2err]
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
        print i
        for h__ in h_:
            h =  h__.Clone('h')
        h.SetTitle("")
        h.Divide(hd) 
        h.GetXaxis().SetTitle('generator p_{T} (GeV)')
        h.GetYaxis().SetTitle('acceptance')
       # h.GetYaxis().SetTitleSize(20)
       # h.GetYaxis().SetTitleFont(43)
        h.GetYaxis().SetLabelOffset(0.010)
        h.GetYaxis().SetTitleOffset(1.45)
        h.GetXaxis().SetLabelOffset(0.010)
        h.GetXaxis().SetTitleOffset(1.45)
       # h.GetYaxis().SetLabelFont(43)
       # h.GetYaxis().SetLabelSize(18)
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
#    for i, h in enumerate(hlist1):         

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
        #print scale
        col = colors[i]
        #print hsOpt['xmin'], hsOpt['xmax']
        h.GetXaxis().SetRangeUser(hsOpt['xmin'],hsOpt['xmax'])
        h.SetMinimum(0.)
        h.Scale(scale)
        h.Rebin(rb)
        h.SetMarkerStyle(8)
        h.SetMarkerSize(0.)
        if dofill[i]: 
            h.SetFillColorAlpha(col,0.1)
            h.SetFillStyle(1) #samOpt['fillstyle']
        else:     
            h.SetLineWidth(1)         
            h.SetLineColorAlpha(col,0.)
    #     if ("csv" in hsOpt['hname'] or "cmva" in hsOpt['hname']): ymax = 0.8
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
#    c1.SaveAs(oDir+"/"+hsOpt['hname']+".root") 

    return [nev1,nev1err,nev2,nev2err]
#------------



def drawH1only(hlist1, snames1, legstack1, hsOpt, oDir, colors, dofill, rebin, headerOpt, isMC):
    gStyle.SetOptStat(False)
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
    if len(legstack1) > 1:
        #for n, leg in enumerate(legstack1):
        #    legend.AddEntry(hlist1[n], leg)
        legend.AddEntry(hlist1[0], legstack1[0])
        legend.AddEntry(hlist1[len(hlist1)-1], legstack1[1])
    else:
        legend.AddEntry(hlist1[len(hlist1)-1], legstack1[0])
    legend.Draw("same")
    #-------------

    hs1.Draw("HISTsame")
    if len(hlist1)>1: herr1.Draw("E2same")
    else: herr1.Draw("Esame")

    nev1 = 0
    nev1err = 0
    if hsOpt['hname']=='h_nevts': 
        nev1 = herr1.GetBinContent(1)
        nev1err = herr1.GetBinError(1)

    c1.Update()    
    c1.SaveAs(oDir+"/"+hsOpt['hname']+".pdf")
    c1.SaveAs(oDir+"/"+hsOpt['hname']+".png")            
#    c1.SaveAs(oDir+"/"+hsOpt['hname']+".root") 

    return [nev1,nev1err]
#-------------------------


def fillH2(h,h1):
    h2 = (TH2D)("h2_bdt","h2_bdt",h.GetNbinsX(),0.,1.,h.GetNbinsX(),0.,1.) 
    for i in range(1,h.GetNbinsX()):
        for j in range(1,h1.GetNbinsX()):
         h2.SetBinContent(i, j, h.GetBinContent(i)+h1.GetBinContent(j))
         #print h2.GetBinContent(i,j)
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
        name = hsOpt['hname'] # + "_" + sname[i]
        print name
        c2 = TCanvas(name, name, 800, 800)       #debug
        print h2.GetBinContent(10)
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
        #palette.GetAxis().SetLabelSize(0.005); #debug
        #palette.GetAxis().SetTitleSize(0.005);
       # drawCMS(12.9, "")
        c2.Update() 

        c2.SaveAs(oDir+"/"+name+".pdf")
        c2.SaveAs(oDir+"/"+name+".png")            
        c2.SaveAs(oDir+"/"+name+".root")  

    return True
#------------


##
# DRAWING TOOLS
#################
def drawCMS(lumi, text, onTop=False ):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextAlign(33)     
    if (type(lumi) is float or type(lumi) is int) and float(lumi) > 0: latex.DrawLatex(0.90, 0.94, "%.1f fb^{-1}  (13 TeV)" % (float(lumi)))
    else: latex.DrawLatex(0.90, 0.94, "simulation (13 TeV)")
    if not onTop: latex.SetTextAlign(11)
    latex.SetTextFont(62)
    latex.SetTextSize(0.03 if len(text)>0 else 0.035)
    if not onTop: 
        latex.DrawLatex(0.15, 0.855, "CMS   "+text)
        latex.SetTextFont(52)
        latex.SetTextSize(0.02)
        latex.DrawLatex(0.15, 0.830, "preliminary")
    else: 
        latex.DrawLatex(0.15, 0.88, "CMS   "+text)
        latex.SetTextFont(52)
        latex.SetTextSize(0.02)
        latex.DrawLatex(0.15, 0.850, "preliminary")
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
        latex.DrawLatex(0.15, 0.820, "Simulation")
        latex.SetTextFont(52);
        latex.SetTextSize(0.030);
        latex.DrawLatex(0.15, 0.780, "pp #rightarrow HH #rightarrow b#bar{b}b#bar{b}");

def setLegend(doRight, doTop):
    leg = TLegend(0.55,0.70,0.90,0.90)
    leg.SetTextSize(0.032)
    return leg
#------------

