# MultivarHermiteManifoldInterp_SISC
This folder contains code associated with with the paper "MULTIVARIATE HERMITE INTERPOLATION ON RIEMANNIAN MANIFOLDS" by Zimmermann/Bergmann,
to appear in SIAM Journal on Scientific Computing.

The code is published in order to ensure the reproducibility of the numerical experiments that feature in this paper.
It performs mutlivariate interpolation via
* the barycentric   Hermite interpolationm approach (BHI),
* the tangent space Hermite interpolation approach (THI).

As hard-coded examples of manifolds, the unit sphere in 3-space and the 3x3 special orthogonal group SO(3) are considered.

As subroutines, it makes use of
* a simplistic implementation of the gradient enhanced Kriging method, see e.g.,
* an elementary Riemannian gradient descent with fixed step size for computing Riemannian barycenters
Both tasks can be tackled with more sophisticated implementations.

The executable main scripts are:

# Hermite_RiemannBary_Helicoid.m  (for BHI, see Section 5.1)

# Hermite_RiemannTang_Helicoid.m  (for THI, see Section 5.1)

# Hermite_RiemannBary_SO.m        (for BHI, see Section 5.2)

# Hermite_RiemannTang_SO.m        (for THI, see Section 5.2)

# 1.) Numerical experiments corresponding to Section 5.1 of the paper:
The script "Hermite_RiemannBary_Helicoid.m" is associated with the numerical experiments on barycentric interpolation in Section 5.1 of the paper.
This script produces Figure 1 (left, middle) and those results of Table 2, that correspond to BHI.

The script "Hermite_RiemannTang_Helicoid.m" is associated with the numerical experiments on tangent space interpolation in Section 5.1 of the paper.
This script produces Figure 1 (left,right) and those results of Table 2, that correspond to THI.


# 2.) Numerical experiments corresponding to Section 5.2 of the paper:
The script "Hermite_RiemannBary_SO.m" is associated with the numerical experiments on barycentric interpolation in Section 5.2 of the paper.
With the selection "cheby". this script produces the left-hand side of Figure 4 and those results of Table 3, that correspond to BHI.

The script "Hermite_RiemannTang_So.m" is associated with the numerical experiments on tangent space interpolation in Section 5.2 of the paper.
With the sample selection "cheby", this script produces the right-hande side of Figure 4 and those results of Table 3, that correspond to THI.

# 3.) How to reproduce the tea pot pictures in Figure 7 of the paper.

To get the list of reference matrices for the figure, first open the script

    "Hermite_RiemannBary_SO.m",
    
set "n1 = 7", set the boolean variable "do_midpoint_trials = 1" and run it.

This will create the matrix array "ref_mats" in the matlab workspace.

Then go to the folder "aux4test" and run the script "plot_t_pot_interp".

