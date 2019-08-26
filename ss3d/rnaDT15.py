class DT15(object):
    # U_BL
    BL_PS = 23.0
    BL_SB = 10.0
    BL_SP = 64.0

    # U_BA
    BA_PSB =  5.0
    BA_PSP = 20.0
    BA_BSP =  5.0
    BA_SPS = 20.0

    # U_EV
    #exv_dist = 3.2
    # Other distances are defined in another place. See <<<< DT15_exv_param 
    #exv_coef = 1.0
    #exv_adjust = 1.5852
    #n_sep_nlocal_P = 3
    #n_sep_nlocal_S = 3
    #n_sep_nlocal_B = 2
    #exv_inf = 0.2

    # U_ST (consecutive)
    ST_DIST = 1.4
    ST_DIH = 4.0
    #*** st_u0 is defined in another place. See <<<< DT15_stack_param

    # U_ST (non-consecutive)
    TST_DIST = 5.00
    TST_ANGL = 1.50
    TST_DIH = 0.15
    TST_U0 = -6.50

    # U_HB
    HB_DIST = 5.0
    HB_ANGL = 1.5
    HB_DIH_HBOND = 0.15
    HB_DIH_CHAIN = 0.15
    HB_U0 = -2.93168175
    #hb_cutoff_dist = 2.0

#<<<< DT15_stack_param  
#**** 
#**** U0 = - h + kB ( T - Tm ) * s
#**** 
#**    h        s       Tm
#AA  3.99578  -0.319  299.058
#AC  3.96314  -0.319  299.058
#AG  4.75714   5.301  341.349
#AU  3.96314  -0.319  299.058
#CA  3.92173  -0.319  299.058
#CC  3.65726  -1.567  285.968
#CG  4.23498   0.774  315.673
#CU  3.62896  -1.567  285.968
#GA  4.71668   5.301  341.349
#GC  4.71997   4.370  343.363
#GG  5.19674   7.346  366.523
#GU  4.61901   2.924  338.329
#UA  3.92173  -0.319  299.058
#UC  3.62234  -1.567  285.968
#UG  4.67613   2.924  338.329
#UU  2.99569  -3.563  251.733
#>>>>
#
#<<<< DT15_exv_param 
#**    R    epsilon
#P    2.1     0.2
#S    2.9     0.2
#A    2.8     0.2
#G    3.0     0.2
#C    2.7     0.2
#U    2.7     0.2
#Mg2  0.7926  0.894700
#Ca2  1.7131  0.459789
#Cl   1.9480  0.265000
#K    2.6580  0.000328
#Na   1.8680  0.002770
#X1   2.1     0.2
#>>>>
