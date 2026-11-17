# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3      1x1        diag1,run2day
smoke          gx3      1x1        diag1,run2day,evp1d     smoke_gx3_1x1_diag1_run2day
smoke          gx3      4x1        diag1,run2day
smoke          gx3      4x1        diag1,run2day,evp1d     smoke_gx3_4x1_diag1_run2day
smoke          gx3      4x2        diag1,run2day
smoke          gx3      4x2        diag1,run2day,evp1d     smoke_gx3_4x2_diag1_run2day
smoke          gx3      1x1        diag1,run2day,nghost2
smoke          gx3      1x1        diag1,run2day,nghost2,evp1d     smoke_gx3_1x1_diag1_nghost2_run2day
smoke          gx1      32x1       diag1,run2day
smoke          gx1      32x1       diag1,run2day,evp1d     smoke_gx1_32x1_diag1_run2day
smoke          tx1m     20x2       diag1,run2day
smoke          tx1m     20x2       diag1,run2day,evp1d     smoke_tx1m_20x2_diag1_run2day
smoke          tx1m     20x2       diag1,run2day,nghost2
smoke          tx1m     20x2       diag1,run2day,nghost2,evp1d     smoke_tx1m_20x2_diag1_nghost2_run2day
smoke          gbox80   1x1        diag1,run1day,bcopen
smoke          gbox80   1x1        diag1,run1day,bcopen,evp1d     smoke_gbox80_1x1_bcopen_diag1_run1day
smoke          gbox80   8x1        diag1,run1day,bcopen
smoke          gbox80   8x1        diag1,run1day,bcopen,evp1d     smoke_gbox80_8x1_bcopen_diag1_run1day
smoke          gbox80   8x1        diag1,run1day,bccyclic
smoke          gbox80   8x1        diag1,run1day,bccyclic,evp1d     smoke_gbox80_8x1_bccyclic_diag1_run1day
smoke          gbox80   8x1        diag1,run1day,boxclosed
smoke          gbox80   8x1        diag1,run1day,boxclosed,evp1d     smoke_gbox80_8x1_boxclosed_diag1_run1day
smoke          gbox80   8x1        diag1,run1day,bcopen,kmtislands
smoke          gbox80   8x1        diag1,run1day,bcopen,kmtislands,evp1d     smoke_gbox80_8x1_bcopen_diag1_kmtislands_run1day
smoke          gbox128  4x2        diag1,run1day,bcopen
smoke          gbox128  4x2        diag1,run1day,bcopen,evp1d     smoke_gbox128_4x2_bcopen_diag1_run1day
smoke          gx3      4x1        debug,diag1,run2day
smoke          gx3      4x1        debug,diag1,run2day,evp1d     smoke_gx3_4x1_debug_diag1_run2day
smoke          tx1m     20x2       debug,diag1,run1day
smoke          tx1m     20x2       debug,diag1,run1day,evp1d     smoke_tx1m_20x2_debug_diag1_run1day
restart        gx3      4x1        diag1
restart        gx3      4x1        diag1,evp1d     restart_gx3_4x1_diag1
restart        tx1m     20x2       diag1
restart        tx1m     20x2       diag1,evp1d     restart_tx1m_20x2_diag1
