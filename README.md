# Advances in Emulation and the TENSE Framework

Supplementary code for MMATH thesis on _Advances in Emulation and the TENSE Framework_.

Code is organised in `/src` with reusable code and functions in the `/src/lib` sub-directory. All figures from the report along with additional ones can be found in the `/figures` directory and miscellaneous data files are stored in `/data`.

We highlight the files `/src/lib/TENSE.R` and `/src/lib/simple_NS_emulator2.R` which implement the TENSE framework in arbitrary dimensions.


## Contents by Chapter

### 2 Introduction to Emulation

- **Simple BL Emulator:** `/src/lib/simple_BL_emulator.R`
- `/src/1d-gp.R`
- `/src/2d-simple-ble.R`

### 3 Covariance Functions

- **Periodic BL Emulator:** `/src/lib/periodic_BL_emulator.R`
- `/src/2d-stationary-cs.R`
- `/src/1d-nscs.R`
- `/src/1d-periodic-cs.R`
- `/src/robot-arm.R` (relies on `src/models/robot.R`)
- `/src/1d-scale-mixture.R`


### 4 Advanced Emulation and Applications

- **Sequential Minimax Algorithm:** `/src/lib/sequential_minimax.R`
- **Maximin 2D LHD Algorithm:** `/src/lib/maximin_lhd_2d.R`
- **Advanced BL Emulator:** `/src/lib/advanced_BL_emulator.R`
- `/src/1d-adv-emulator.R`
- `/src/extrapolation-problems.R`
- `/src/2d-grid-designs.R`
- `/src/exp-designs.R`
- `/src/sequential-minimax.R`
- `/src/1d-implausibility.R`
- `/src/2d-history-match.R`


### 5 Emulation with Partial Known Discontinuities

- **Simple Non-Stationary Emulator 1:** `/src/lib/simple_NS_emulator.R` (For Implementing TENSE)
- **Simple Non-Stationary Emulator 2:** `/src/lib/simple_NS_emulator2.R` (_improved implementation_)
- **Advanced Non-Stationary Emulator:** `/src/lib/advanced_NS_emulator.R`
- **TENSE for Arbitrary Dimension:** `/src/lib/TENSE.R`
- `/src/2d-simple-discontinuity.R`
- `/src/2d-mult-curved-discontinuities.R`
- `/src/dist-warp-fig.R`
- `/src/3d-TENSE.R`

### Misc

- `/src/lib/filled.contour3.R` (modifies contour plotting function)





