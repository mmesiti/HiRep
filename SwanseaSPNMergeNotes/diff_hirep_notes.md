# Merge Notes
The following document has been produced looking at the output 
of 
git diff -w COMMIT__BEFORE__FORK LATEST_SPN_COMMIT
where
COMMI_BEFORE_FORK = 6b232fe0146751d9d19e2e9d8f2b06003e9b2fc5 
LATEST_SPN_COMMIT = b34b469b6e0bd7d7925e57a5f792c18069993fa9

## Converter +30 -0

### archive_eolexi.c +16 -0
Changes to archive only the upper half of the matrix (compressed format). Added `#ifdef` 
to change range of the for loop.

### converter.c +14 -0
- trigger compilation error if using quaternions with `GAUGE_SPN`
- added formats "fullSPN" and "compSPN" (compressed)
(see /LibHR +920 -103/IO/spn_archive.c +406 -0)

## GaugeFix/gaugefix.c +4 -0 +4 -0
for loop only up to NG/2 in the SPN case

## Include +101 -5
### TMPL/suN_repr_func.h.tmpl +0 -0
added _`algebra_project_FMAT` macro (see autosun)

### communications.h +5 -0
**For the Clover term**
In the case of SPN in the fundamental representation, use `suNffuld_field` 
(which is the full-matrix type) instead of `suNf_field`, which is the compressed
case, in the calls to `complete_clover_force_sendrecv` and 
`start_clover_force_sendrecv`.

### complex.h +49 -0
Defined a handful more of convenience macros.

### global.h +8 -0
Defining the clover term global variables `cl_term` and `cl_force` as
`suNffull_field` when using GAUGE_SPN in the fundamental representation.
Added global variable for debug.

### io.h +6 -0
Added functions to read full matrices when using GAUGE_SPN

### memory.h +12 -3
`suNffull_field` type used for clover-related allocations and deallocations 
instead of `suNfc_field` when using GAUGE_SPN in the fundamental representation.

### moreio.h +5 -0
Added functions to read full matrices when using GAUGE_SPN

### update.h +3 -0
Added debug function

### utils.h +5 -1
clover term-related logic uses `suNffull_field` when using GAUGE_SPN and 
REPR_FUNDAMENTAL

## LibHR +920 -103
### Geometry/communications.c +12 -1
clover term-related logic uses `suNffull_field` when using GAUGE_SPN and 
REPR_FUNDAMENTAL

### IO +428 -5
#### archive_su2quat.c +18 -1
added functions to translate between quaternion and compressed spn format.

#### logger.c +4 -4
Changes unrelated to SPN (initialised pointers to null).

#### spn_archive.c +406 -0
Totally new file. SPN-full implementations of 
* write_gauge_field
* write_gauge_field_matrix
* read_gauge_field
* read_gauge_field_matrix
with a '_fullSPN' suffix. These functions are only used in converter.c    

### Memory +20 -0
#### amalloc.c +14 -0
Debugging function added - unrelated to SPN

#### field_alloc.c +6 -0
Declaration of memory function for clover terms must use `suNffull_field`.
(clover term-related logic uses `suNffull_field` when using GAUGE_SPN and 
REPR_FUNDAMENTAL)

### Random/random_suNg.c +37 -2
Adaptation to produce SPN random matrices 

### Update +328 -86
#### Dphi.c +9 -9
Changes to `Cphi_`. (_)
We need to use full matrix operations instead of compressed matrix ones.

#### cabmar.c +86 -1
implementation of the Cabibbo-Marinari algorithm for SPN

#### clover_tools.c +16 -7
clover term-related logic uses `suNffull` when using GAUGE_SPN and 
REPR_FUNDAMENTAL (clover loop is fine with compressed matrices, though).


#### force_fermion_core.c +55 -45
* clover term-related logic uses `suNffull` when using GAUGE_SPN and 
REPR_FUNDAMENTAL . Calculations done in compressed formats need to be expanded.
Full matrix operations are needed here in the clover, spn, fundamental case.
* in force_fermion_core, 
  * suNf -> suNf_FMAT
  * _algebra_project -> _algebra_project_FMAT

#### force_hmc_ff.c +0 -0
in force_hmc_ff
* suNf -> suNf_FMAT
* _algebra_project -> _algebra_project_FMAT

#### force_scalar.c +4 -0
in the calculation of the outer product, the loop goes only until NG * NG / 2

#### luscherweisz.c +79 -20
* redefinition of S (shift by a constant) in `test_wilson_action_and_force()`, 
  used only in Tests.
* loop on generators NG * NG - 1 -> NG * (NG+1) /2 for SPN 
(changes to make the tests work again?)

#### random_momenta.c +1 -3
Redefinition of the number of generators as
```
const int ngen=sizeof(suNg_algebra_vector)/sizeof(double);
```
#### representation.c +6 -1
For SPN, the function _group_represent2 is just a thin wrapper around the macro
_group_represent.

#### update_mt.c +72 -0
new function `update_ghmc_stripped` for debug purposes.

### Utils +95 -9
#### HYP_smearing.c +2 -2
Some functions not working for SPN are removed from compilation or throw an 
error at runtime.

#### TMPL/suN_exp.c.tmpl +24 -1
Added taylor exponentiation in the template

#### boundary_conditions.c +25 -1
clover term-related logic uses `suNffull_field` when using GAUGE_SPN and 
REPR_FUNDAMENTAL

#### det_suNg.c +8 -2
In the SPN case the determinant is computed expanding to a full matrix and 
then reusing the existing suN code.

#### inv_suNg.c +9 -2
In the SPN case the inverse is computed expanding to a full matrix and 
then reusing the existing suN code.

#### suN_utils.c +26 -1
rewrite of project_to_suNg for the SPN case.

## Make +1929 -63
#### MkRules +7 -0
Guard against use of SP2

### Utils +1908 -48
#### autosun +538 -30
##### adjoint.h +14 -1
* in init(): case for spn added. dimension of the algebra.
* group_represent(): case for spn added, uses spmatrix

##### antisymmetric.h +41 -0
* the dimension of the representation for SPN is N * (N-1)/2-1, one of 
  the generators must be removed from the usual construction (the omega-like
  one?)
* group_represent(): case for spn added, uses spmatrix

##### fundamental.h +12 -1
* in group_represent, added case for SPN, with compressed version

##### list.h +9 -0
Added header guards and a bunch of includes.

##### matrix.h +81 -0
* Added header guards and a bunch of includes.
* Added string representation for compressed assignment
* Definition of spmatrix

##### polynomial.h +6 -0
Added header guards and a bunch of includes.

##### representation.h +24 -0
Added header guards and a bunch of includes.

##### sparse.h +9 -0
Added header guards and a bunch of includes.

##### sun.h +323 -27
* Added header guards and a bunch of includes.
* from `#ifdef` to `switch()` for group type.
* written luscher exponentiation for SPN 
* Other trivial changes.

##### symmetric.h +10 -1
* group_represent(): case for spn added, uses spmatrix

#### write_suN_headers.pl +1366 -14
Written the code for generating spn-compressed macros.

## RenormalizationFactors/measure_Z_mom.c +2 -0 +2 -0 
just logging in main() for SPN

## Spectrum +6 -0
### measure_formfactor.c +2 -0
just logging in main() for SPN

### measure_spectrum.c +2 -0
just logging in main() for SPN

### mk_mesons_with_z2semwall_new.c +2 -0
just logging in main() for SPN

## TestProgram +2423 -80 
A number of test programs were broken. Some of them were fixed, and some
have been added (see, e.g., SPNtoSUNRegression) with some convenience scripts.

## WilsonFlow +20 -4
### WF_measure.c +2 -0
just logging in main() for SPN

### WF_measure_adaptative.c +2 -0
just logging in main() for SPN

### wilsonflow.c +16 -4
Changes in matrix size NG * NG -> NG * NG / 2 and algebra dimension 
(NG * NG - 1) -> (NG * (NG + 1) / 2) for the SPN case

