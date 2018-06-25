#bms=("20171120-160644-bm0" "20171120-161211-bm1" "20171120-161730-bm2" "20171120-162255-bm3" "20171120-162812-bm4" "20171120-163332-bm5" "20171120-163852-bm6" "20171120-164402-bm7" "20171120-164916-bm8" "20171120-165436-bm9" "20171120-165954-bm10" "20171120-170526-bm11" "20171120-171044-bm12" "20171120-171608-bm13")
bm_string="20171120-160644-bm0"
bms=(0) # 1 2 3 4 5 6 7 8 9 10 11 12 13)
i=0
#for bm in ${bms[@]}
#do
#    #python scripts/drawcomp_afterBDT.py -b $bm -l 0 -w 6 --res 2 -n
#	#python scripts/drawcomp_afterBDT.py -b $bm -l 0 -w 5 -c --res -1
#    #python scripts/drawcomp_postFit.py -b $bm -l 1 -d reports_with_bias_weight
#	python scripts/drawcomp_afterBDT.py -b $bm -l 0 -w 2072 --res 10 -d ttbar_bias_check/BM$i
#	((i++))
#done
for bm in ${bms[@]}
do
    #python scripts/drawcomp_afterBDT.py -b $bm -l 0 -w 6 --res 2 -n
	python scripts/drawcomp_afterBDT.py -b $bm_string -l 0 -w 51 --res 5  -d classifier_reports/reports_SM_no_bias_corr_mixing_fix/BM$bm    #Figure 2
    #python scripts/drawcomp_afterBDT.py -b $bm_string -l 0 -w 52 --res 5  -d classifier_reports/reports_SM_no_bias_corr_btagside/BM$bm    #Figure 2 btagside
    #python scripts/drawcomp_afterBDT.py -b $bm_string -l 0 -w 5 --res 2  -d classifier_reports/reports_sm/BM$bm

    #python scripts/drawcomp_postFit.py -b $bm_string --bm $bm -l 0 -d classifier_reports/reports_no_bias_corr_SM_mixing_fix_0_2/BM$bm --res -12
    
    #python scripts/drawcomp_postFit.py -b $bm_string --bm $bm -l 1 -d classifier_reports/reports_SM_mixing_fix/BM$bm   #All vars after final fit
    python scripts/drawcomp_postFit.py -b $bm_string --bm $bm -l 0 -d classifier_reports/reports_SM_mixing_fix/BM$bm #Final BDT
    #python scripts/drawcomp_afterBDT.py -b $bm -l 0 -w 2072 --res 10 -d ttbar_bias_check/BM$bm
	
	#Requested by ARC
	#python scripts/drawcomp_afterBDT.py -b $bm_string -w 5 -c --res -1 -d classifier_reports/reports_SM_mixing_fix_nocut_MS_bias_corr/BM$bm
	#python scripts/drawcomp_afterBDT.py -b $bm_string -w 5 -c --res -1 -d classifier_reports/reports_SM_mixing_fix_nocut_MS_bias_corr/BM$bm -r 8
	#python scripts/drawcomp_afterBDT.py -b $bm_string -w 6 -c --res -1 -d classifier_reports/reports_SM_btagside/BM$bm
	#python scripts/drawcomp_afterBDT.py -b $bm_string -w 6 -c --res -1 -d classifier_reports/reports_SM_btagside/BM$bm -r 8

    #CR plots:
    #python scripts/drawcomp_afterBDT.py -b $bm_string -w 6 --res -1 -d classifier_reports/reports_SM_btagside/BM$bm
    #python scripts/drawcomp_afterBDT.py -b $bm_string -w 5 --res -1 -d classifier_reports/reports_SM_nocut/BM$bm
done
