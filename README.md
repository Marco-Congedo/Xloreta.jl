
![header](Documents/header.png)

---

> [!TIP] 
> 🦅
> This package is part of the [Eegle.jl](https://github.com/Marco-Congedo/Eegle.jl) ecosystem for EEG data analysis and classification.

---

# Xloreta

This is a pure-[julia](https://julialang.org/) 100%-human package for computing, testing and using human EEG 
(Electroencephalography) inverse solutions of the *Minimum Norm* family. Particularly, it implements the following vector-type distributed inverse solutions:
- *weighted minimum norm* — see [^1] for a review,
- *standardized Low-Resolution Electromagnetic Tomography* (sLORETA) [^2],
- *exact Low-Resolution Electromagnetic Tomography* (eLORETA) [^3].

For all of them, the usual *model-driven* and the *data-driven* [^4] versions are provided, with the latter being actually *beamformers* like in [^5] [^6] [^7] and being little known in the current literature.

> All mathematical details can be found in [^1] [^2] [^3] [^4] [^8].
>
> An overview of the formula involved in the implementation is [here](https://drive.google.com/file/d/1QDKkZUuiY1sNz7kw2mqp7itObRYUxk1x/view?usp=drive_web).
>
> Those that are not familiar with the material, may want to start with this [introduction.](https://drive.google.com/file/d/0B_albC6Y6I9KczRoNjlsbWxKZ3c/view?usp=drive_web&resourcekey=0-LJGNC8sOIGlft_FJ565muA)
> 
> An actual leadfield matrix that can be used in this package is available in [Leadfields.jl](https://github.com/Marco-Congedo/Leadfields.jl).


![separator](Documents/separator.png)

## 🧭 Index

- 📦 [Installation](#-installation)
- 🔣 [Problem Statement, Notation and Nomenclature](#-problem-statement-notation-and-nomenclature)
- 🔌 [API](#-api)
- 💡 [Examples](#-examples)
- ✍️ [About the Author](#️-about-the-author)
- 🌱 [Contribute](#-contribute)
- 🎓 [References](#-references)

![separator](Documents/separator.png)

## 📦 Installation

*julia* version 1.10+ is required.

Execute the following command in julia's REPL:

```julia
using Pkg
Pkg.add("Xloreta")
```

To test the package:
```julia
Pkg.test("Xloreta")
```

[▲ index](#-index)

![separator](Documents/separator.png)

## 🔣 Problem Statement, Notation and Nomenclature

We are given an EEG sensor potentials measurement

𝐱(𝑡) ∈ ℝⁿ 

at n electrodes referenced to the common average, in μV units, where 𝑡 is discrete time (samples);

we wish to estimate the current density

𝐣(𝑡) ∈ ℝᵖ 

at p cortical grey matter voxels, in A/m² units, in the three Cartesian spatial directions (x, y, z).

We have therefore:

**Forward equation** — determining the scalp voltage given the current distribution:

𝐱(𝑡) = 𝐊 𝐣(𝑡).

It is unique for a given leadfield matrix 

𝐊 ∈ ℝⁿ×³ᵖ.

Each column of the leadfield is the scalp field for unit-length dipole pointing in one of three orthogonal directions. The leadfield encapsulates a physical head model [^9] [^10].

**Inverse solution** — estimating the current distribution given the scalp voltage:

𝐣(𝑡) = 𝐓 𝐱(𝑡).

It is not unique. Each inverse solution method yields a different transfer matrix

𝐓 ∈ ℝ³ᵖ×ⁿ,

the computation of which is the main purpose of this package. For details on the composition of the 𝐊 and 𝐓 matrices the see [^8].

> [!IMPORTANT] 
> A solution is said *genuine* or to *respect the measurement* if 
> 𝐊 𝐓 = 𝐇ₙ,
> where 𝐇ₙ is the [centeringMatrix](#centeringmatrix) in dimension n (since the leadfield is always to be centered).
> The weighted minimum norm and eLORETA are genuine solutions, while sLORETA is not.
>
> Also, matrix 𝐓 𝐊 is called the resolution matrix [^11]. Its successive groups of three columns, one group per voxel, are called the point-spread functions. 
> They allow to ascertain whether the transfer matrix is capable of correctly localizing a single current dipole, regardless of its position (voxel) and orientation.
>
> This is a minimal localization capability for an inverse solution, as it (unrealistically) assumes the absence of noise in the measurement and the existence of only one active dipole at a time. Nonetheless, it is a minimal requirement. sLORETA and eLORETA possess this property, while the minimum norm does not, like most inverse solution methods found in the literature, and thus should not be used.

> [!TIP] 
> Overall, eLORETA is known to be an excellent choice, as it provides stable results without requiring fine-tuned noise level estimation [^12].

> [!WARNING] 
> Throughout this documentation and in the package it is always assumed both the input data and the leadfield matrix are referenced to the common average — see [centeringMatrix](#centeringmatrix).

[▲ index](#-index)

![separator](Documents/separator.png)

## 🔌 API

The package exports the following functions:

| Function | Description |
|:---------|:---------|
| [centeringMatrix](#centeringmatrix)   | common average reference operator (alias: ℌ)   |
| [c2cd](#cd2sm)                        | compute the squared magnitude of the current density given a current density vector |
| [psfLocError](#psflocerror)           | point spread function localization error   |
| [psfErrors](#psferrors)               | point spread function localization, spread and equalization errors |
| [minnorm](#minnorm)                   | compute minimum norm transfer matrix (model and data-driven) |
| [sLORETA](#sloreta)                   | compute sLORETA transfer matrix (model and data-driven)|
| [eLORETA](#eloreta)                   | compute eLORETA transfer matrix (model and data-driven)|

[▲ index](#-index)

![separator](Documents/separator.png)

### centeringMatrix

```julia
function centeringMatrix(N::Int)
```
The common average reference (CAR) operator for referencing EEG data potentials so that their mean across sensors (space) is zero at all samples.

Let 𝐗 be the s × n EEG recording, where s and n denote the number of samples and channels (sensors), respectively, and let 𝐇ₙ be the n×n centering matrix, then

𝐘 = 𝐗 𝐇ₙ

is the CAR (or centered) data.

𝐇ₙ is named the common average reference operator. It is given by 

𝐇ₙ = 𝐈ₙ − (1/n) (𝟭ₙ 𝟭ₙᵀ),

where 𝐈ₙ is the n-dimensional identity matrix and 𝟭ₙ is the n-dimensional vector of ones — see for example p.67 in [^13].

Alias ℌ (U+210C, with escape sequence "frakH")

Return the n×n centering matrix.

[▲ API index](#-api)

[▲ index](#-index)

---
### cd2sm

```julia
function cd2sm(j::Vector{R}) where R<:Real
```

> 'cd2sm' stands for 'current density to squared magnitude'. 

Return the current density squared magnitude vector comprised of 1/3 of the elements of the input current density vector `j`. The current density vector `j` holds successively the triplets (x, y, z). Return the successive sums (x²+y²+z²) for each triplet.

The input vector `j` may contain any exact multiple of 3 number of elements.

> [!NOTE] 
> Typically, the squared magnitude of the current density is the quantity of interest in neuroimaging studies.

[▲ API index](#-api)

[▲ index](#-index)

---

### psfLocError

```julia
psfLocError(K::Matrix{R}, T::Matrix{R}) where R<:Real
```

> 'psfLocError' stands for 'point spread functions Localization Error'

Given a n × 3p leadfield matrix `K` and an associated 3p × n transfer matrix `T`, where n is the number of electrodes and 3p is the number of p voxels times 3 (the x, y, z source components), return the number of localization errors obtained by point spread functions — see 🔣 [here](#-problem-statement-notation-and-nomenclature).

Any time you create a sLORETA or eLORETA transfer matrix `T` for a given leadfield matrix `K`, you should test it with this function — see the 💡 [examples](#-examples).

[▲ API index](#-api)

[▲ index](#-index)

---

### psfErrors

```julia
function psfErrors(K::Matrix{R}, T::Matrix{R}) where R<:Real
```
> 'psfErrors' stands for 'point spread function Errors'

Given a n × 3p leadfield matrix `K` and an associated 3p × n transfer matrix `T`, where n is the number of electrodes and 3p is the number of p voxels times 3 (the x, y, z source components), return the 3-tuple of vectors holding 3p errors obtained for each component (x, y, z) at each voxel (test locations):

1. *Localization errors* (bool vector),
    true if the maximum current density magnitude is not located in the test component/location, false otherwise.
2. *Spread errors* (Float vector),
    log of the sum of current density squared magnitude in the entire volume divided by the current density squared magnitude in the test component/location.
3. *Equalization errors* (Float vector),
    uncorrected variance of the current density squared magnitude across the entire volume (all locations, i.e., all voxels)

[▲ API index](#-api)

[▲ index](#-index)

---

### minNorm

```julia
function minNorm(K::Matrix{R},
                 α::Real=0.,
                 C::Union{Symbol, Symmetric{R}, Hermitian{R}}=:modelDriven;
                 W::Union{Vector{R}, Nothing}=nothing) where R<:Real
```

Given a n × 3p leadfield matrix, where n is the number of electrodes and 3p is the number of p voxels times 3 (the x, y, z source components),
return the **minimum norm regularized transfer matrix** with regularization `α`.

if `C` is `:modelDriven` (default), compute the model driven solution, otherwise `C` must be the data covariance matrix and in this case compute the
data-driven solution — see [here](https://drive.google.com/file/d/1QDKkZUuiY1sNz7kw2mqp7itObRYUxk1x/view?usp=drive_web). If `C` is a matrix, it must be flagged as `Symmetric` or as `Hermitian`  — see the 💡 [examples](#-examples).

if optional keyword argument `W` is a vector of 3p non-negative weights, compute the **weighted minimum norm solution** instead. In this case `C` must be
equal to `:modelDriven` (default), as a weighted data-driven solution is not defined.

> [!IMPORTANT] 
> If passed as a matrix, `C` must be non-singular. No check is performed.
> 
> The columns of the leadfield matrix must be centered (common average reference).
> 
> A suitable regularization parameter `α`> 0 should be found by cross-validation or any other suitable method. Never assume an arbitrary value is suitable.
> If a pre-whitening by a noise covariance is sought, apply the pre-whitening to the leadfield, then apply the transfer matrix to pre-whitened data.

[▲ API index](#-api)

[▲ index](#-index)

---

### sLORETA

```julia
function sLORETA(K::Matrix{R},
                 α::Real=0.,
                 C::Union{Symbol, Symmetric{R}, Hermitian{R}}=:modelDriven) where R<:Real
```

Given a n × 3p leadfield matrix, where n is the number of electrodes and 3p is the number of p voxels times 3 (the x, y, z source components),
return the **sLORETA transfer matrix** with regularization `α`.

if `C` is `:modelDriven` (default), compute the model driven solution, otherwise `C` must be the data covariance matrix and in this case compute the
data-driven solution, which is similar (actually better) to the linearly constrained minimum variance beamformer — see [here](https://drive.google.com/file/d/1QDKkZUuiY1sNz7kw2mqp7itObRYUxk1x/view?usp=drive_web). If `C` is a matrix, it must be flagged as `Symmetric` or as `Hermitian`  — see the 💡 [examples](#-examples).

> [!IMPORTANT] 
> If passed as a matrix, `C` must be non-singular. No check is performed.
> 
> The columns of the leadfield matrix must be centered (common average reference).
> 
> A suitable regularization parameter `α`> 0 should be found by cross-validation or any other suitable method. Never assume an arbitrary value is suitable.
> If a pre-whitening by a noise covariance is sought, apply the pre-whitening to the leadfield, then apply the transfer matrix to pre-whitened data.

[▲ API index](#-api)

[▲ index](#-index)

---

### eLORETA

```julia
function eLORETA(K::Matrix{R},
                 α::Real=0.,
                 C::Union{Symbol, Symmetric{R}, Hermitian{R}}=:modelDriven,
                 tol::Real=0.0;
            verbose=true) where R<:Real
```

Given a n × 3p leadfield matrix, where n is the number of electrodes and 3p is the number of p voxels times 3 (the x, y, z source components),
return the **eLORETA transfer matrix** with regularization `α`.

if `C` is `:modelDriven` (default), compute the model driven solution, otherwise `C` must be the data covariance matrix and in this case compute the
data-driven solution — see [here](https://drive.google.com/file/d/1QDKkZUuiY1sNz7kw2mqp7itObRYUxk1x/view?usp=drive_web). If `C` is a matrix, it must be flagged as `Symmetric` or as `Hermitian`  — see the 💡 [examples](#-examples).

The model-driven solution is iterative; the convergence at each iteration is printed unless optional keyword argument `verbose` is set to false.

`tol` is the tolerance for establishing convergence; it defaults to the square root of `Base.eps` of the nearest type of the elements of `K`. This corresponds to requiring the average norm of the difference between the 3 × 3 diagonal blocks of the weight matrix in two successive iterations
to vanish for about half the significant digits.

> [!IMPORTANT] 
> If passed as a matrix, `C` must be non-singular. No check is performed.
> 
> The columns of the leadfield matrix must be centered (common average reference).
> 
> A suitable regularization parameter `α`> 0 should be found by cross-validation or any other suitable method. Never assume an arbitrary value is suitable.
> If a pre-whitening by a noise covariance is sought, apply the pre-whitening to the leadfield, then apply the transfer matrix to pre-whitened data.

[▲ API index](#-api)

[▲ index](#-index)

![separator](Documents/separator.png)

## 💡 Examples

To use this package, all you will need is here below, where it is understood that you replace the example data matrix `X` and the leadfield matrix `K` with your own data and leadfield. To access a standard leadfield, see [here](https://github.com/Marco-Congedo/Leadfields.jl):

```julia
using Xloreta, LinearAlgebra

# example number of electrodes, data samples, voxels
n, s, p = 20, 200, 2000

# example random leadfield in common average reference
K = ℌ(n)*randn(n, 3p)

# example random EEG data
X=randn(s, n)

# random weights for weighted minimum norm solutions
weights=abs.(randn(3p))

# - - -

# sample covariance matrix of the random EEG data
C=Symmetric((1/s)*(X'*X))

Tmn1 = minNorm(K, 1)    # unweighted model-driven min norm with α=1
Tmn2 = minNorm(K, 10)   # unweighted model-driven min norm with α=10
Tmn3 = minNorm(K, 1; W=weights) # weighted model-driven min norm with α=1
Tmn4 = minNorm(K, 1, C) # data-driven min norm with α=1

TsLor1 = sLORETA(K, 1)     # model-driven sLORETA with α=1
TsLor2 = sLORETA(K, 10)    # model-driven sLORETA with α=10
TsLor3 = sLORETA(K, 1, C)  # data-driven sLORETA with α=1

TeLor1 = eLORETA(K, 1)     # model-driven eLORETA with α=1
TeLor2 = eLORETA(K, 10)    # model-driven eLORETA with α=10
TeLor3 = eLORETA(K, 1, C)  # data-driven eLORETA with α=1

# test the transfer matrix you creat
psfLocError(K, TeLor1) == 0 ? println("OK") : println("Error")
```

[▲ index](#-index)

![separator](Documents/separator.png)

## ✍️ About the Author

[Marco Congedo](https://sites.google.com/site/marcocongedo) is a Research Director of [CNRS](http://www.cnrs.fr/en) (Centre National de la Recherche Scientifique), working at [UGA](https://www.univ-grenoble-alpes.fr/english/) (University of Grenoble Alpes). **Contact**: first name dot last name at gmail dot com.

[▲ index](#-index)

![separator](Documents/separator.png)

## 🌱 Contribute

Please contact the author if you are interested in contributing.

[▲ index](#-index)

![separator](Documents/separator.png)

## 🎓 References

[^1]: R. D., Pascual-Marqui, “Review of Methods for Solving the EEG Inverse Problem”, Int. J. Bioelectromag., vol. 1, no.1, pp. 75-86; 1999. [pdf](https://www.ijbem.org/volume1/number1/77-90.pdf).

[^2]: R. D. Pascual-Marqui, “Standardized Low Resolution brain electromagnetic Tomography (sLORETA): technical details,” Methods Find. Exp. Clin. Pharmacol., vol 24D, pp. 5-12, 2002. [pdf](https://www.uzh.ch/keyinst/NewLORETA/sLORETA/sLORETA-Math01.pdf).

[^3]: R.D. Pascual-Marqui, "Discrete, 3D distributed, linear imaging methods of electric neuronal activity. Part 1: exact, zero
error localization," arXiv:0710.3341, 2007-October-17. [pdf](http://arxiv.org/pdf/0710.3341)

[^4]: G. Lio, "Identifier l’activité cérébrale responsable de l’état de stress et d’anxiété induit par les acouphènes chroniques. Une étude comparative de l’activité spectrale des aires de Brodmann par tomographie cérébrale électromagnétique chez le sujet sain et pathologique", Mémoire de deuxième année de master en sciences humaines et sociales, Université Grenoble Alpes, 2010. [pdf](https://osf.io/te2j9/files/s5b7e).

[^5]: B. D. van Veen, W. van Drongelen, A. Suzuki, “Localization of Brain Electrical Activity via Linearly Constrained Minimum Variance Spatial
Filter,” IEEE Trans. Biomed. Eng., vol. 44, no. 9, pp. 867-880, 1997.

[^6]: K. Sekihara, M Sahani, S.S: Nagarajan, “Localization Bias and Spatial Resolution of Adaptive and non-Adaptive Spatial Filters for MEG
Source Reconstruction,” Neuroimage, vol. 25, pp. 1056-1067, 2005.

[^7]: R.E. Greenblatt, A. Ossadtchi, M.E. Pflieger, "Local Linear Estimators
for the Bioelectromagnetic Inverse Problem," IEEE Trans. Sig. Pro., vol53, no. 9, pp. 3403-3412, 2005.

[^8]: M. Congedo, "Subspace Projection Filters for Real-Time Brain Electromagnetic Imaging," IEEE Transactions on Biomedical Engineering, 53 (8), pp. 1624-34, 2006. [pdf](https://hal.science/hal-00460510v1/document)

[^9]: F. Lopes da Silva, "Functional Localization of Brain Sources using EEG and/or MEG data: Volume Conductor and Source Models," Magn. Res.
Img., vol. 22, pp. 1533-1538, 2004.

[^10]: J. Sarvas, “Basic Mathematical and Electromagnetic Concepts of the
Biomagnetic Inverse Problem,” Phys. Med. Biol., vol 32, no. 1, pp. 1122, 1987.

[^11]: G. Backus, F. Gilbert, “The resolving power of gross earth data,”Geophys. J. R. Asr. Soc, vol. 16, pp. 169-205, 1968.

[^12]: A Negi, S Haufe, A Gramfort, A Hashemi, "How forgiving are M/EEG inverse solutions to noise level misspecification? An excursion into the BSI-Zoo,"
bioRxiv, 2025.03. 12.642831. [pdf](https://www.biorxiv.org/content/biorxiv/early/2025/03/13/2025.03.12.642831.full.pdf).

[^13]: S.R. Searle, “Matrix Algebra Useful for Statistics,” John Wiley & Sons, New-York, 1982.

[▲ index](#-index)