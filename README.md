## PROSAIL_D_Forward

This folder is used for forward simulations based on the **PROSAIL model**.  
By varying key vegetation parameters, the directional reflectance of the canopy is analyzed, providing support for subsequent parameter inversion and LUT construction.

---

### PROSAIL_Main_Cab.m
**Purpose:**  To analyze the effect of leaf chlorophyll content (Cab) variation on canopy reflectance.
**Variable parameter setting:**  `Cab`: 10 – 60 μg/cm²
**Output:**  Directional reflectance curves in the principal plane and perpendicular plane at **600 nm** under different Cab conditions.

### PROSAIL_Main_LAI.m
**Purpose:**  To analyze the effect of leaf area index (LAI) variation on canopy reflectance.
**Variable parameter setting:**  `LAI`: 0.5 – 4.0
**Output:**  Directional reflectance curves in the principal plane and perpendicular plane at **800 nm** under different LAI conditions.
