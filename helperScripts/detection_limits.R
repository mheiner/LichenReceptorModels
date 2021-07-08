# from:
# Rachel Buck
# EAL Manager
# Plant & Wildlife Sciences
# Email: rachelbuck@byu.edu
# DOCUMENT: iCAP 7400 Detection limits.pdf (Radial measurements)
# conversion: one migrogram per liter equals 0.001 parts per million

detection_limits_micgrogramperliter = c(Ag=2.46, Al=1.51, As=4.74, B=1.26, Ba=0.17, 
                                        Be=0.07, Ca=0.02, Cd=0.19, Co=1.16, Cr=0.85,
                                        Cu=2.36, Fe=0.80, Hg=1.10, K=5.10, Li=0.83, 
                                        Mg=0.04, Mn=0.21, Mo=1.11, Na=1.80, Ni=2.29,
                                        P=5.66, Pb=4.50, S=2.22, Sb=9.36, Se=7.36, 
                                        Si=7.20, Sn=1.57, Sr=0.04, Ti=0.58, Tl=7.33, 
                                        V=0.80, Zn=0.60)

detection_limits_ppmDilute = detection_limits_micgrogramperliter * 100.0 / 1000.0
