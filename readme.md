# Fast Macroscopic Forcing Method (Fast MFM)

Code for: __*Bryngelson, S. H., *Schäfer, F., Liu, J., & Mani, A. (2024). Fast Macroscopic Forcing Method. Journal of Computational Physics, 499, 112721. *Authors contributed equally.__

## Abstract

The macroscopic forcing method (MFM) of Mani and Park and similar methods for obtaining
turbulence closure operators, such as the Green’s function-based approach of Hamba, recover
reduced solution operators from repeated direct numerical simulations (DNS). MFM has already
been used to successfully quantify Reynolds-averaged Navier–Stokes (RANS)-like operators for
homogeneous isotropic turbulence and turbulent channel flows. Standard algorithms for MFM force
each coarse-scale degree of freedom (i.e., degree of freedom in the RANS space) and conduct a
corresponding fine-scale simulation (i.e., DNS), which is expensive. We combine this method with
an approach recently proposed by Schäfer and Owhadi to recover elliptic integral operators from
a polylogarithmic number of matrix–vector products. The resulting fast MFM introduced in this
work applies sparse reconstruction to expose local features in the closure operator and reconstructs
this coarse-grained differential operator in only a few matrix–vector products and correspondingly, a
few MFM simulations. For flows with significant nonlocality, the algorithm first "peels" long-range
effects with dense matrix–vector products to expose a more local operator. We demonstrate the
algorithm’s performance for scalar transport in a laminar channel flow and momentum transport in
a turbulent channel flow. For these problems, we recover eddy–diffusivity- and eddy–viscosity-like
operators, respectively, at 1% of the cost of computing the exact operator via a brute-force approach
for the laminar channel flow problem and 13% for the turbulent one. We observe that we can
reconstruct these operators with an increase in accuracy by about a factor of 100 over randomized
low-rank methods. Applying these operators to compute the averaged fields of interest has visually
indistinguishable behavior from the exact solution. Our results show that a similar number of
simulations are required to reconstruct the operators to the same accuracy under grid refinement.
Thus, the accuracy corresponds to the physics of the problem, not the numerics. We glean that
for problems in which the RANS space is reducible to one dimension, eddy diffusivity and eddy
viscosity operators can be reconstructed with reasonable accuracy using only a few simulations,
regardless of simulation resolution or degrees of freedom.

## Citation

```bibtex
@article{bryngelsonfmfm24,
  author = {{*}Bryngelson, S. H. and {*}Sch{\"a}fer, F. and Liu, J. and Mani, A.},
  title = {Fast macroscopic forcing method},
  journal = {Journal of Computational Physics},
  year = {2024},
  volume = {499},
  pages = {112721},
  doi = {10.1016/j.jcp.2023.112721},
  note = {{*}Authors contributed equally}
}
```

## License

MIT
