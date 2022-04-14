
# Part 1: write the model equations ---------------------------------------

#### Write the equations of your model, followng RxODE format


model_RxODE <-RxODE({
  
  
  ####### Fixed values #######
  kdeg_Bcl2 <- 0.139
  kprod_Bcl2 <- Bcl20 * kdeg_Bcl2
  kdeg_Bclxl <- 0.139
  kprod_Bclxl <- Bclxl0 * kdeg_Bclxl
  kdeg_Mcl1 <- 0.925
  kprod_Mcl1 <- Mcl10 * kdeg_Mcl1
  kdeg_BIM <- 0.173
  kdeg_tBID <- 0.554
  kdeg_PUMA <- 0.204
  kdeg_NOXA <- 0.695
  kdeg_Bcl2_BIM <- 0.554
  kforward_Bcl2_BIM <- 0.108
  kbackward_Bcl2_BIM <- 0.504
  kdeg_Bclxl_BIM <- 0.554
  kforward_Bclxl_BIM <- 1.98
  kbackward_Bclxl_BIM <- 1.58
  kdeg_Mcl1_BIM <- 0.277
  kforward_Mcl1_BIM <- 4.68
  kbackward_Mcl1_BIM <- 0.936
  kdeg_Bcl2_tBID <- 0.554
  kforward_Bcl2_tBID <- 0.126
  kbackward_Bcl2_tBID <- 0.504
  kdeg_Bclxl_tBID <- 0.554
  kforward_Bclxl_tBID <- 0.0583
  kbackward_Bclxl_tBID <- 1.58
  kdeg_Mcl1_tBID <- 0.554
  kforward_Mcl1_tBID <- 0.0947
  kbackward_Mcl1_tBID <- 0.936
  kdeg_Bcl2_PUMA <- 0.554
  kforward_Bcl2_PUMA <- 0.028
  kbackward_Bcl2_PUMA <- 0.504
  kdeg_Bclxl_PUMA <- 0.554
  kforward_Bclxl_PUMA <- 0.311
  kbackward_Bclxl_PUMA <- 1.58
  kdeg_Mcl1_PUMA <- 0.277
  kforward_Mcl1_PUMA <- 0.493
  kbackward_Mcl1_PUMA <- 0.936
  kdeg_Bcl2_NOXA <- 0.554
  kforward_Bcl2_NOXA <- 0.00262
  kbackward_Bcl2_NOXA <- 0.504
  kdeg_Bclxl_NOXA <- 0.554
  kforward_Bclxl_NOXA <- 0.000158
  kbackward_Bclxl_NOXA <- 1.58
  kdeg_Mcl1_NOXA <- 0.925
  kforward_Mcl1_NOXA <- 0.0237
  kbackward_Mcl1_NOXA <- 0.936
  kforward_Bcl2_BAKa <- 5.04e-05
  kbackward_Bcl2_BAKa <- 0.504
  kforward_Bcl2_BAXma <- 0.0336
  kbackward_Bcl2_BAXma <- 0.504
  kforward_Bclxl_BAKa <- 0.0198
  kbackward_Bclxl_BAKa <- 1.58
  kforward_Bclxl_BAXma <- 0.00186
  kbackward_Bclxl_BAXma <- 1.58
  kforward_Mcl1_BAKa <- 0.117
  kbackward_Mcl1_BAKa <- 0.936
  kforward_Mcl1_BAXma <- 9.36e-05
  kbackward_Mcl1_BAXma <- 0.936
  kforward_BIM_BAXc <- 0.00925
  kbackward_BIM_BAXc <- 0.925
  kforward_tBID_BAXc <- 0.00925
  kbackward_tBID_BAXc <- 0.925
  kforward_PUMA_BAXc <- 0.00925
  kbackward_PUMA_BAXc <- 0.925
  kforward_BIM_BAK <- 0.00925
  kbackward_BIM_BAK <- 0.925
  kforward_tBID_BAK <- 0.00925
  kbackward_tBID_BAK <- 0.925
  kforward_PUMA_BAK <- 0.00925
  kbackward_PUMA_BAK <- 0.925
  k_BIM_BAK <- 418
  k_tBID_BAK <- 418
  k_PUMA_BAK <- 418
  k_BIM_BAXc <- 418
  k_tBID_BAXc <- 418
  k_PUMA_BAXc <- 418
  k_BAXca <- 418
  kforward_BAK_VDAC2 <- 0.00832
  kbackward_BAK_VDAC2 <- 8.32
  kforward_BAK2 <- 0.0461
  kbackward_BAK2 <- 0.695
  kforward_BAK4 <- 0.0461
  kbackward_BAK4 <- 0.695
  kforward_BAK6 <- 0.0461
  kbackward_BAK6 <- 0.695
  kforward_BAK8 <- 0.0461
  kbackward_BAK8 <- 0.695
  kforward_BAK10 <- 0.0461
  kbackward_BAK10 <- 0.695
  kforward_BAK12 <- 0.0461
  kbackward_BAK12 <- 0.695
  kforward_BAK8 <- 0.0461
  kbackward_BAK8 <- 0.695
  kforward_BAK10 <- 0.0461
  kbackward_BAK10 <- 0.695
  kforward_BAK12 <- 0.0461
  kbackward_BAK12 <- 0.695
  kforward_BAK12 <- 0.0461
  kbackward_BAK12 <- 0.695
  kforward_BAX2 <- 0.0461
  kbackward_BAX2 <- 0.695
  kforward_BAX4 <- 0.0461
  kbackward_BAX4 <- 0.695
  kforward_BAX6 <- 0.0461
  kbackward_BAX6 <- 0.695
  kforward_BAX8 <- 0.0461
  kbackward_BAX8 <- 0.695
  kforward_BAX10 <- 0.0461
  kbackward_BAX10 <- 0.695
  kforward_BAX12 <- 0.0461
  kbackward_BAX12 <- 0.695
  kforward_BAX8 <- 0.0461
  kbackward_BAX8 <- 0.695
  kforward_BAX10 <- 0.0461
  kbackward_BAX10 <- 0.695
  kforward_BAX12 <- 0.0461
  kbackward_BAX12 <- 0.695
  kforward_BAX12 <- 0.0461
  kbackward_BAX12 <- 0.695
  kdeg_ApoG2_Bcl2 <- 0.554
  kforward_ApoG2_Bcl2 <- 0.0198
  kbackward_ApoG2_Bcl2 <- 0.695
  kdeg_ApoG2_Bclxl <- 0.554
  kforward_ApoG2_Bclxl <- 0.00105
  kbackward_ApoG2_Bclxl <- 0.695
  kdeg_ApoG2_Mcl1 <- 0.277
  kforward_ApoG2_Mcl1 <- 0.0277
  kbackward_ApoG2_Mcl1 <- 0.695
  kdeg_ABT737_Bcl2 <- 0.554
  kforward_ABT737_Bcl2 <- 0.698
  kbackward_ABT737_Bcl2 <- 0.695
  kdeg_ABT737_Bclxl <- 0.554
  kforward_ABT737_Bclxl <- 1.39
  kbackward_ABT737_Bclxl <- 0.695
  kdeg_ABT737_Mcl1 <- 0.277
  kforward_ABT737_Mcl1 <- 0.000709
  kbackward_ABT737_Mcl1 <- 0.695
  kdeg_ApoG2 <- 0.0868
  kdeg_ABT737 <- 0.347
  kdeg_Bcl2I_Bcl2 <- 0.554
  kforward_Bcl2I_Bcl2 <- 1.39
  kbackward_Bcl2I_Bcl2 <- 0.695
  kdeg_Bcl2I <- 0.0868
  kdeg_Bclxl_Bclxl <- 0.554
  kforward_Bclxl_Bclxl <- 1.39
  kbackward_Bclxl_Bclxl <- 0.695
  kdeg_BclxlI <- 0.0868
  kdeg_Mcl1I_Mcl1 <- 0.554
  kforward_Mcl1I_Mcl1 <- 1.39
  kbackward_Mcl1I_Mcl1 <- 0.695
  kdeg_Mcl1I <- 0.0868
  ####### End fixed values #######
  
  
  ####### Initial conditions #######
  Bcl2(0) <- Bcl20
  Bclxl(0) <- Bclxl0
  Mcl1(0) <- Mcl10
  BAXc(0) <- BAXc0
  BAK(0) <- BAK0
  VDAC2(0) <- BAK0
  ####### End initial conditions #######
  
  
  
  ####### All rates #######
  R1_deg_Bcl2 <- - kdeg_Bcl2 * Bcl2
  
  R2_prod_Bcl2 <- kprod_Bcl2
  
  R3_deg_Bclxl <- - kdeg_Bclxl * Bclxl
  
  R4_prod_Bclxl <- kprod_Bclxl
  
  R5_deg_Mcl1 <- - kdeg_Mcl1 * Mcl1
  
  R6_prod_Mcl1 <- kprod_Mcl1
  
  R7_deg_BIM <- - kdeg_BIM * BIM
  
  R8_deg_tBID <- - kdeg_tBID * tBID
  
  R9_deg_PUMA <- - kdeg_PUMA * PUMA
  
  R10_deg_NOXA <- - kdeg_NOXA * NOXA
  
  R11_deg_Bcl2_BIM <- - kdeg_Bcl2_BIM * Bcl2_BIM
  
  R12_complex_Bcl2_BIM <- kforward_Bcl2_BIM * Bcl2 * BIM - kbackward_Bcl2_BIM * Bcl2_BIM
  
  R13_deg_Bclxl_BIM <- - kdeg_Bclxl_BIM * Bclxl_BIM
  
  R14_complex_Bclxl_BIM <- kforward_Bclxl_BIM * Bclxl * BIM - kbackward_Bclxl_BIM * Bclxl_BIM
  
  R15_deg_Mcl1_BIM <- - kdeg_Mcl1_BIM * Mcl1_BIM
  
  R16_complex_Mcl1_BIM <- kforward_Mcl1_BIM * Mcl1 * BIM - kbackward_Mcl1_BIM * Mcl1_BIM
  
  R17_deg_Bcl2_tBID <- - kdeg_Bcl2_tBID * Bcl2_tBID
  
  R18_complex_Bcl2_tBID <- kforward_Bcl2_tBID * Bcl2 * tBID - kbackward_Bcl2_tBID * Bcl2_tBID
  
  R19_deg_Bclxl_tBID <- - kdeg_Bclxl_tBID * Bclxl_tBID
  
  R20_complex_Bclxl_tBID <- kforward_Bclxl_tBID * Bclxl * tBID - kbackward_Bclxl_tBID * Bclxl_tBID
  
  R21_deg_Mcl1_tBID <- - kdeg_Mcl1_tBID * Mcl1_tBID
  
  R22_complex_Mcl1_tBID <- kforward_Mcl1_tBID * Mcl1 * tBID - kbackward_Mcl1_tBID * Mcl1_tBID
  
  R23_deg_Bcl2_PUMA <- - kdeg_Bcl2_PUMA * Bcl2_PUMA
  
  R24_complex_Bcl2_PUMA <- kforward_Bcl2_PUMA * Bcl2 * PUMA - kbackward_Bcl2_PUMA * Bcl2_PUMA
  
  R25_deg_Bclxl_PUMA <- - kdeg_Bclxl_PUMA * Bclxl_PUMA
  
  R26_complex_Bclxl_PUMA <- kforward_Bclxl_PUMA * Bclxl * PUMA - kbackward_Bclxl_PUMA * Bclxl_PUMA
  
  R27_deg_Mcl1_PUMA <- - kdeg_Mcl1_PUMA * Mcl1_PUMA
  
  R28_complex_Mcl1_PUMA <- kforward_Mcl1_PUMA * Mcl1 * PUMA - kbackward_Mcl1_PUMA * Mcl1_PUMA
  
  R29_deg_Bcl2_NOXA <- - kdeg_Bcl2_NOXA * Bcl2_NOXA
  
  R30_complex_Bcl2_NOXA <- kforward_Bcl2_NOXA * Bcl2 * NOXA - kbackward_Bcl2_NOXA * Bcl2_NOXA
  
  R31_deg_Bclxl_NOXA <- - kdeg_Bclxl_NOXA * Bclxl_NOXA
  
  R32_complex_Bclxl_NOXA <- kforward_Bclxl_NOXA * Bclxl * NOXA - kbackward_Bclxl_NOXA * Bclxl_NOXA
  
  R33_deg_Mcl1_NOXA <- - kdeg_Mcl1_NOXA * Mcl1_NOXA
  
  R34_complex_Mcl1_NOXA <- kforward_Mcl1_NOXA * Mcl1 * NOXA - kbackward_Mcl1_NOXA * Mcl1_NOXA
  
  R35_complex_Bcl2_BAKa <- kforward_Bcl2_BAKa * Bcl2 * BAKa - kbackward_Bcl2_BAKa * Bcl2_BAKa
  
  R36_complex_Bcl2_BAXma <- kforward_Bcl2_BAXma * Bcl2 * BAXma - kbackward_Bcl2_BAXma * Bcl2_BAXma
  
  R37_complex_Bclxl_BAKa <- kforward_Bclxl_BAKa * Bclxl * BAKa - kbackward_Bclxl_BAKa * Bclxl_BAKa
  
  R38_complex_Bclxl_BAXma <- kforward_Bclxl_BAXma * Bclxl * BAXma - kbackward_Bclxl_BAXma * Bclxl_BAXma
  
  R39_complex_Mcl1_BAKa <- kforward_Mcl1_BAKa * Mcl1 * BAKa - kbackward_Mcl1_BAKa * Mcl1_BAKa
  
  R40_complex_Mcl1_BAXma <- kforward_Mcl1_BAXma * Mcl1 * BAXma - kbackward_Mcl1_BAXma * Mcl1_BAXma
  
  R41_complex_BIM_BAXc <- kforward_BIM_BAXc * BIM * BAXc - kbackward_BIM_BAXc * BIM_BAXc
  
  R42_complex_tBID_BAXc <- kforward_tBID_BAXc * tBID * BAXc - kbackward_tBID_BAXc * tBID_BAXc
  
  R43_complex_PUMA_BAXc <- kforward_PUMA_BAXc * PUMA * BAXc - kbackward_PUMA_BAXc * PUMA_BAXc
  
  R44_complex_BIM_BAK <- kforward_BIM_BAK * BIM * BAK - kbackward_BIM_BAK * BIM_BAK
  
  R45_complex_tBID_BAK <- kforward_tBID_BAK * tBID * BAK - kbackward_tBID_BAK * tBID_BAK
  
  R46_complex_PUMA_BAK <- kforward_PUMA_BAK * PUMA * BAK - kbackward_PUMA_BAK * PUMA_BAK
  
  R47_disso_BIM_BAK <-  - k_BIM_BAK * BIM_BAK
  
  R48_disso_tBID_BAK <-  - k_tBID_BAK * tBID_BAK
  
  R49_disso_PUMA_BAK <-  - k_PUMA_BAK * PUMA_BAK
  
  R50_disso_BIM_BAXc <-  - k_BIM_BAXc * BIM_BAXc
  
  R51_disso_tBID_BAXc <-  - k_tBID_BAXc * tBID_BAXc
  
  R52_disso_PUMA_BAXc <-  - k_PUMA_BAXc * PUMA_BAXc
  
  R53_disso_BAXca <-  - k_BAXca * BAXca
  
  R54_complex_BAK_VDAC2 <- kforward_BAK_VDAC2 * BAK * VDAC2 - kbackward_BAK_VDAC2 * BAK_VDAC2
  
  R55_complex_BAK2 <- kforward_BAK2 * BAKa * BAKa - kbackward_BAK2 * BAK2
  
  R56_complex_BAK4 <- kforward_BAK4 * BAK2 * BAK2 - kbackward_BAK4 * BAK4
  
  R57_complex_BAK6 <- kforward_BAK6 * BAK2 * BAK4 - kbackward_BAK6 * BAK6
  
  R58_complex_BAK8 <- kforward_BAK8 * BAK2 * BAK6 - kbackward_BAK8 * BAK8
  
  R59_complex_BAK10 <- kforward_BAK10 * BAK2 * BAK8 - kbackward_BAK10 * BAK10
  
  R60_complex_BAK12 <- kforward_BAK12 * BAK2 * BAK10 - kbackward_BAK12 * BAK12
  
  R61_complex_BAK8 <- kforward_BAK8 * BAK4 * BAK4 - kbackward_BAK8 * BAK8
  
  R62_complex_BAK10 <- kforward_BAK10 * BAK4 * BAK6 - kbackward_BAK10 * BAK10
  
  R63_complex_BAK12 <- kforward_BAK12 * BAK4 * BAK8 - kbackward_BAK12 * BAK12
  
  R64_complex_BAK12 <- kforward_BAK12 * BAK6 * BAK6 - kbackward_BAK12 * BAK12
  
  R65_complex_BAX2 <- kforward_BAX2 * BAXma * BAXma - kbackward_BAX2 * BAX2
  
  R66_complex_BAX4 <- kforward_BAX4 * BAX2 * BAX2 - kbackward_BAX4 * BAX4
  
  R67_complex_BAX6 <- kforward_BAX6 * BAX2 * BAX4 - kbackward_BAX6 * BAX6
  
  R68_complex_BAX8 <- kforward_BAX8 * BAX2 * BAX6 - kbackward_BAX8 * BAX8
  
  R69_complex_BAX10 <- kforward_BAX10 * BAX2 * BAX8 - kbackward_BAX10 * BAX10
  
  R70_complex_BAX12 <- kforward_BAX12 * BAX2 * BAX10 - kbackward_BAX12 * BAX12
  
  R71_complex_BAX8 <- kforward_BAX8 * BAX4 * BAX4 - kbackward_BAX8 * BAX8
  
  R72_complex_BAX10 <- kforward_BAX10 * BAX4 * BAX6 - kbackward_BAX10 * BAX10
  
  R73_complex_BAX12 <- kforward_BAX12 * BAX4 * BAX8 - kbackward_BAX12 * BAX12
  
  R74_complex_BAX12 <- kforward_BAX12 * BAX6 * BAX6 - kbackward_BAX12 * BAX12
  
  R75_deg_ApoG2_Bcl2 <- - kdeg_ApoG2_Bcl2 * ApoG2_Bcl2
  
  R76_complex_ApoG2_Bcl2 <- kforward_ApoG2_Bcl2 * ApoG2 * Bcl2 - kbackward_ApoG2_Bcl2 * ApoG2_Bcl2
  
  R77_deg_ApoG2_Bclxl <- - kdeg_ApoG2_Bclxl * ApoG2_Bclxl
  
  R78_complex_ApoG2_Bclxl <- kforward_ApoG2_Bclxl * ApoG2 * Bclxl - kbackward_ApoG2_Bclxl * ApoG2_Bclxl
  
  R79_deg_ApoG2_Mcl1 <- - kdeg_ApoG2_Mcl1 * ApoG2_Mcl1
  
  R80_complex_ApoG2_Mcl1 <- kforward_ApoG2_Mcl1 * ApoG2 * Mcl1 - kbackward_ApoG2_Mcl1 * ApoG2_Mcl1
  
  R81_deg_ABT737_Bcl2 <- - kdeg_ABT737_Bcl2 * ABT737_Bcl2
  
  R82_complex_ABT737_Bcl2 <- kforward_ABT737_Bcl2 * ABT737 * Bcl2 - kbackward_ABT737_Bcl2 * ABT737_Bcl2
  
  R83_deg_ABT737_Bclxl <- - kdeg_ABT737_Bclxl * ABT737_Bclxl
  
  R84_complex_ABT737_Bclxl <- kforward_ABT737_Bclxl * ABT737 * Bclxl - kbackward_ABT737_Bclxl * ABT737_Bclxl
  
  R85_deg_ABT737_Mcl1 <- - kdeg_ABT737_Mcl1 * ABT737_Mcl1
  
  R86_complex_ABT737_Mcl1 <- kforward_ABT737_Mcl1 * ABT737 * Mcl1 - kbackward_ABT737_Mcl1 * ABT737_Mcl1
  
  R87_deg_ApoG2 <- - kdeg_ApoG2 * ApoG2
  
  R88_deg_ABT737 <- - kdeg_ABT737 * ABT737
  
  R89_deg_Bcl2I_Bcl2 <- - kdeg_Bcl2I_Bcl2 * Bcl2I_Bcl2
  
  R90_complex_Bcl2I_Bcl2 <- kforward_Bcl2I_Bcl2 * Bcl2I * Bcl2 - kbackward_Bcl2I_Bcl2 * Bcl2I_Bcl2
  
  R91_deg_Bcl2I <- - kdeg_Bcl2I * Bcl2I
  
  R92_deg_Bclxl_Bclxl <- - kdeg_Bclxl_Bclxl * Bclxl_Bclxl
  
  R93_complex_Bclxl_Bclxl <- kforward_Bclxl_Bclxl * Bclxl * Bclxl - kbackward_Bclxl_Bclxl * Bclxl_Bclxl
  
  R94_deg_BclxlI <- - kdeg_BclxlI * BclxlI
  
  R95_deg_Mcl1I_Mcl1 <- - kdeg_Mcl1I_Mcl1 * Mcl1I_Mcl1
  
  R96_complex_Mcl1I_Mcl1 <- kforward_Mcl1I_Mcl1 * Mcl1I * Mcl1 - kbackward_Mcl1I_Mcl1 * Mcl1I_Mcl1
  
  R97_deg_Mcl1I <- - kdeg_Mcl1I * Mcl1I
  ####### End All rates #######
  
  
  ## Brut addition of perfusion##

  if(t > 50 & t <62 ){
    BIMperf <- BIM0 / 12 
    PUMAperf <- PUMA0 / 12 
    NOXAperf <- NOXA0 / 12
    
  }else{
    
    BIMperf <- 0
    PUMAperf <- 0
    NOXAperf <- 0
    
  }
 ## End addition of perfusion 
  
  

  
  ####### All ODES #######
  
  d/dt(Bcl2) <- R1_deg_Bcl2 + R2_prod_Bcl2 - R12_complex_Bcl2_BIM - R18_complex_Bcl2_tBID - R24_complex_Bcl2_PUMA - R30_complex_Bcl2_NOXA - R35_complex_Bcl2_BAKa - R36_complex_Bcl2_BAXma - R76_complex_ApoG2_Bcl2 - R82_complex_ABT737_Bcl2 - R90_complex_Bcl2I_Bcl2
  d/dt(Bclxl) <- R3_deg_Bclxl + R4_prod_Bclxl - R14_complex_Bclxl_BIM - R20_complex_Bclxl_tBID - R26_complex_Bclxl_PUMA - R32_complex_Bclxl_NOXA - R37_complex_Bclxl_BAKa - R38_complex_Bclxl_BAXma - R78_complex_ApoG2_Bclxl - R84_complex_ABT737_Bclxl - R93_complex_Bclxl_Bclxl - R93_complex_Bclxl_Bclxl
  d/dt(Mcl1) <- R5_deg_Mcl1 + R6_prod_Mcl1 - R16_complex_Mcl1_BIM - R22_complex_Mcl1_tBID - R28_complex_Mcl1_PUMA - R34_complex_Mcl1_NOXA - R39_complex_Mcl1_BAKa - R40_complex_Mcl1_BAXma - R80_complex_ApoG2_Mcl1 - R86_complex_ABT737_Mcl1 - R96_complex_Mcl1I_Mcl1
  d/dt(BIM) <- BIMperf + R7_deg_BIM - R12_complex_Bcl2_BIM - R14_complex_Bclxl_BIM - R16_complex_Mcl1_BIM - R41_complex_BIM_BAXc - R44_complex_BIM_BAK - R47_disso_BIM_BAK - R50_disso_BIM_BAXc
  d/dt(tBID) <- R8_deg_tBID - R18_complex_Bcl2_tBID - R20_complex_Bclxl_tBID - R22_complex_Mcl1_tBID - R42_complex_tBID_BAXc - R45_complex_tBID_BAK - R48_disso_tBID_BAK - R51_disso_tBID_BAXc
  d/dt(PUMA) <- PUMAperf +  R9_deg_PUMA - R24_complex_Bcl2_PUMA - R26_complex_Bclxl_PUMA - R28_complex_Mcl1_PUMA - R43_complex_PUMA_BAXc - R46_complex_PUMA_BAK - R49_disso_PUMA_BAK - R52_disso_PUMA_BAXc
  d/dt(NOXA) <- NOXAperf + R10_deg_NOXA - R30_complex_Bcl2_NOXA - R32_complex_Bclxl_NOXA - R34_complex_Mcl1_NOXA
  d/dt(Bcl2_BIM) <- R11_deg_Bcl2_BIM + R12_complex_Bcl2_BIM
  d/dt(Bclxl_BIM) <- R13_deg_Bclxl_BIM + R14_complex_Bclxl_BIM
  d/dt(Mcl1_BIM) <- R15_deg_Mcl1_BIM + R16_complex_Mcl1_BIM
  d/dt(Bcl2_tBID) <- R17_deg_Bcl2_tBID + R18_complex_Bcl2_tBID
  d/dt(Bclxl_tBID) <- R19_deg_Bclxl_tBID + R20_complex_Bclxl_tBID
  d/dt(Mcl1_tBID) <- R21_deg_Mcl1_tBID + R22_complex_Mcl1_tBID
  d/dt(Bcl2_PUMA) <- R23_deg_Bcl2_PUMA + R24_complex_Bcl2_PUMA
  d/dt(Bclxl_PUMA) <- R25_deg_Bclxl_PUMA + R26_complex_Bclxl_PUMA
  d/dt(Mcl1_PUMA) <- R27_deg_Mcl1_PUMA + R28_complex_Mcl1_PUMA
  d/dt(Bcl2_NOXA) <- R29_deg_Bcl2_NOXA + R30_complex_Bcl2_NOXA
  d/dt(Bclxl_NOXA) <- R31_deg_Bclxl_NOXA + R32_complex_Bclxl_NOXA
  d/dt(Mcl1_NOXA) <- R33_deg_Mcl1_NOXA + R34_complex_Mcl1_NOXA
  d/dt(Bcl2_BAKa) <- R35_complex_Bcl2_BAKa
  d/dt(BAKa) <- -R35_complex_Bcl2_BAKa - R37_complex_Bclxl_BAKa - R39_complex_Mcl1_BAKa - R47_disso_BIM_BAK - R48_disso_tBID_BAK - R49_disso_PUMA_BAK - R55_complex_BAK2 - R55_complex_BAK2
  d/dt(Bcl2_BAXma) <- R36_complex_Bcl2_BAXma
  d/dt(BAXma) <- -R36_complex_Bcl2_BAXma - R38_complex_Bclxl_BAXma - R40_complex_Mcl1_BAXma - R53_disso_BAXca - R65_complex_BAX2 - R65_complex_BAX2
  d/dt(Bclxl_BAKa) <- R37_complex_Bclxl_BAKa
  d/dt(Bclxl_BAXma) <- R38_complex_Bclxl_BAXma
  d/dt(Mcl1_BAKa) <- R39_complex_Mcl1_BAKa
  d/dt(Mcl1_BAXma) <- R40_complex_Mcl1_BAXma
  d/dt(BIM_BAXc) <- R41_complex_BIM_BAXc + R50_disso_BIM_BAXc
  d/dt(BAXc) <- -R41_complex_BIM_BAXc - R42_complex_tBID_BAXc - R43_complex_PUMA_BAXc
  d/dt(tBID_BAXc) <- R42_complex_tBID_BAXc + R51_disso_tBID_BAXc
  d/dt(PUMA_BAXc) <- R43_complex_PUMA_BAXc + R52_disso_PUMA_BAXc
  d/dt(BIM_BAK) <- R44_complex_BIM_BAK + R47_disso_BIM_BAK
  d/dt(BAK) <- -R44_complex_BIM_BAK - R45_complex_tBID_BAK - R46_complex_PUMA_BAK - R54_complex_BAK_VDAC2
  d/dt(tBID_BAK) <- R45_complex_tBID_BAK + R48_disso_tBID_BAK
  d/dt(PUMA_BAK) <- R46_complex_PUMA_BAK + R49_disso_PUMA_BAK
  d/dt(BAXca) <- -R50_disso_BIM_BAXc - R51_disso_tBID_BAXc - R52_disso_PUMA_BAXc + R53_disso_BAXca
  d/dt(BAK_VDAC2) <- R54_complex_BAK_VDAC2
  d/dt(VDAC2) <- -R54_complex_BAK_VDAC2
  d/dt(BAK2) <- R55_complex_BAK2 - R56_complex_BAK4 - R56_complex_BAK4 - R57_complex_BAK6 - R58_complex_BAK8 - R59_complex_BAK10 - R60_complex_BAK12
  d/dt(BAK4) <- R56_complex_BAK4 - R57_complex_BAK6 - R61_complex_BAK8 - R61_complex_BAK8 - R62_complex_BAK10 - R63_complex_BAK12
  d/dt(BAK6) <- R57_complex_BAK6 - R58_complex_BAK8 - R62_complex_BAK10 - R64_complex_BAK12 - R64_complex_BAK12
  d/dt(BAK8) <- R58_complex_BAK8 - R59_complex_BAK10 + R61_complex_BAK8 - R63_complex_BAK12
  d/dt(BAK10) <- R59_complex_BAK10 - R60_complex_BAK12 + R62_complex_BAK10
  d/dt(BAK12) <- R60_complex_BAK12 + R63_complex_BAK12 + R64_complex_BAK12
  d/dt(BAX2) <- R65_complex_BAX2 - R66_complex_BAX4 - R66_complex_BAX4 - R67_complex_BAX6 - R68_complex_BAX8 - R69_complex_BAX10 - R70_complex_BAX12
  d/dt(BAX4) <- R66_complex_BAX4 - R67_complex_BAX6 - R71_complex_BAX8 - R71_complex_BAX8 - R72_complex_BAX10 - R73_complex_BAX12
  d/dt(BAX6) <- R67_complex_BAX6 - R68_complex_BAX8 - R72_complex_BAX10 - R74_complex_BAX12 - R74_complex_BAX12
  d/dt(BAX8) <- R68_complex_BAX8 - R69_complex_BAX10 + R71_complex_BAX8 - R73_complex_BAX12
  d/dt(BAX10) <- R69_complex_BAX10 - R70_complex_BAX12 + R72_complex_BAX10
  d/dt(BAX12) <- R70_complex_BAX12 + R73_complex_BAX12 + R74_complex_BAX12
  d/dt(ApoG2_Bcl2) <- R75_deg_ApoG2_Bcl2 + R76_complex_ApoG2_Bcl2
  d/dt(ApoG2) <- -R76_complex_ApoG2_Bcl2 - R78_complex_ApoG2_Bclxl - R80_complex_ApoG2_Mcl1 + R87_deg_ApoG2
  d/dt(ApoG2_Bclxl) <- R77_deg_ApoG2_Bclxl + R78_complex_ApoG2_Bclxl
  d/dt(ApoG2_Mcl1) <- R79_deg_ApoG2_Mcl1 + R80_complex_ApoG2_Mcl1
  d/dt(ABT737_Bcl2) <- R81_deg_ABT737_Bcl2 + R82_complex_ABT737_Bcl2
  d/dt(ABT737) <- -R82_complex_ABT737_Bcl2 - R84_complex_ABT737_Bclxl - R86_complex_ABT737_Mcl1 + R88_deg_ABT737
  d/dt(ABT737_Bclxl) <- R83_deg_ABT737_Bclxl + R84_complex_ABT737_Bclxl
  d/dt(ABT737_Mcl1) <- R85_deg_ABT737_Mcl1 + R86_complex_ABT737_Mcl1
  d/dt(Bcl2I_Bcl2) <- R89_deg_Bcl2I_Bcl2 + R90_complex_Bcl2I_Bcl2
  d/dt(Bcl2I) <- -R90_complex_Bcl2I_Bcl2 + R91_deg_Bcl2I
  d/dt(Bclxl_Bclxl) <- R92_deg_Bclxl_Bclxl + R93_complex_Bclxl_Bclxl
  d/dt(BclxlI) <- R94_deg_BclxlI
  d/dt(Mcl1I_Mcl1) <- R95_deg_Mcl1I_Mcl1 + R96_complex_Mcl1I_Mcl1
  d/dt(Mcl1I) <- -R96_complex_Mcl1I_Mcl1 + R97_deg_Mcl1I
  
  
  
  pctBAX <- 100 * (6 * BAX6 + 8 * BAX8 + 10 * BAX10 + 12 * BAX12)/BAXc0
  pctBAK <- 100 * (6 * BAK6 + 8 * BAK8 + 10 * BAK10 + 12 * BAK12)/BAK0
  Pore <- 100 * (6 * BAK6 + 8 * BAK8 + 10 * BAK10 + 12 * BAK12 + 6 * BAX6 + 8 * BAX8 + 10 * BAX10 + 12 * BAX12)/(BAK0 + BAXc0)
  
  
  ## Add death signal
  
  if(Pore > 10){
    
    TimeAboveRate <- 1
  }else{
    
    TimeAboveRate <- 0
    
  }
  
    d/dt(TimeAbove) <- TimeAboveRate
  
  
})



## Verification: perform verification number1


# Part 2: parameters, initial states and time measurement ---------------------------------------

## You need to fill the following section.
## Easiest way is to first eval "model_extract()" and copy paste the output
## To directly have pre-fille the right parameter defaults values and initial states
## You can of course do modifications if needed


parameters_default_values <- c(
  
  kforward_BAK10 = 0.0461
)

# Now the initial values, with can be parameter quoted

initial_cmt_values <- c(
  
  tBID_BAK = 0
)

times <- c(0:100) # times you want to see your observations


protocols <- list( unique = tibble(cmt = "Bcl2", time = 0, amt = 0))

## Verification: perform verification number2

#### HELPER - TO COMMENT AFTER USE ############


# Step1: launch the followind command (tribblecreator) 
#               tribblecreator(model_RxODE)

# Step2 copy paste code printed in the console, then 
# 1) fill the domain table, 
# 2) replace your_output_without_quotes by your outputes (eg. tumVol, Conc, ...), and give a protocol
# name provided in protocols

### replace this line by the code procuded in step1 

# domain <- tribble(~param, ~min, ~ref,  ~max,
#                   "Bcl20",0,500,1000 ,
#                   "Bclxl0",0,500,1000 ,
#                   "Mcl10", 0,50,100,
#                   "BIM0", 0,500,1000,
#                   "PUMA0",0,500,1000 ,
#                   "NOXA0", 0,500,1000,
#                   "BAXc0",0,500,1000 ,
#                   "BAK0", 0,500,1000)
# 
# find_relative(Pore, protocol = "unique", model = model_RxODE, values = domain, sensitivity = 0.1)
# 
# find_relative(Pore, protocol = "unique", model = model_RxODE, values = domain, sensitivity = 0.1, deepAnalysis = T)
# 
# domain <- tribble(~param, ~min, ~ref,  ~max,
#                   "Bcl20",0,500,1000 ,
#                   "Bclxl0",0,500,1000 ,
#                   "Mcl10", 0,500,1000,
#                   "BIM0", 0,500,1000,
#                   "PUMA0",0,500,1000 ,
#                   "NOXA0", 0,500,1000,
#                   "BAXc0",0,500,1000 ,
#                   "BAK0", 0,500,1000)
# 
# find_relative(Pore, protocol = "unique", model = model_RxODE, values = domain, sensitivity = 0.1)
# 
# find_relative(Pore, protocol = "unique", model = model_RxODE, values = domain, time_simul = 1:100, sensitivity = 0.1, deepAnalysis = T)

# Step3: verify the plot, the table summarising the influence, and if you find it
# plausible copy paste the three line (param_reduce, param_increase, param_no_impact)
# in the next bloc. Feel free to modify manually if needed
#### END HELPER - Please comment the full bloc above ############
## Verification: perform verification number3


param_reduce<- list(TimeAbove = c("Bcl20", "Bclxl0", "Mcl10"), Pore = c("Bcl20", "Bclxl0", "Mcl10")) #conc1 not here because....


param_increase <- list(TimeAbove =  c( "BIM0", "PUMA0", "NOXA0" ), Pore =  c( "BIM0", "PUMA0", "NOXA0" ))

param_no_impact <- list(TimeAbove = character())

# Part 4: data and concentration to test -------------------------------------

# data shoud have at least ID, Value and concX columns, X being replace by drug number (one col by drug concentration)
# Avoid any column starting with "conc" if it is not a drug concentration / dose column

data_VT <- read.table("D:/these/Second_project/QSP/VirtualTumor/datademo_rework.csv", sep = ";", header = T) %>%
  as_tibble

