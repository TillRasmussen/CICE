# Halo-blocked EVP subcycling (`nghost` > 1)

**Status: experimental. Correct and bit-for-bit at `nghost=1`; NOT yet correct
at `nghost>1`.** See "Known defect" below before enabling it.

## What it does

The EVP solver exchanges `uvel`/`vvel` after every one of `ndte` subcycles.
At `ndte=240` that is 240 latency-bound neighbour exchanges per timestep on
the B grid, and more on C/CD where several fields are exchanged per pass.
Communication, not arithmetic, dominates the solver at scale -- which is the
same observation that motivates the `shared_mem_1d` (evp1d) algorithm.

A subcycle consumes exactly **one ring of halo**. Writing the pass as
`stress` then `stepu`:

```
str (i,j)  = g( uvel(i-1:i  , j-1:j  ) )     ! stress reads W, SW, S
uvel'(i,j) = f( str (i  :i+1, j  :j+1) )     ! stepu  reads E, NE, N
```

so `uvel'` at `(i,j)` depends on `uvel` over `(i-1:i+1, j-1:j+1)` -- a
symmetric one-ring stencil. With `nghost` valid rings after an exchange,
`nghost` passes can therefore run before the next exchange while leaving the
physical domain exact. Exchanges drop from `ndte` to `ceil(ndte/nghost)`.

The trade is redundant computation in the halo instead of communication,
which is a good trade when latency dominates.

## Implementation

Two changes only.

**1. `ice_dyn_shared.F90`, `dyn_prep2`** -- the `icellT` and `icellU` loop
bounds widen by `nghost-1` rings in each direction:

```fortran
do j = jlo-nghost+1, jhi+nghost        ! icellT, was jlo, jhi+1
do j = jlo-nghost+1, jhi+nghost-1      ! icellU, was jlo, jhi
```

Without this the outer rings are never advanced, and a later pass in a block
reads a ring that still holds its exchanged (pass-zero) value.

**2. `ice_dyn_evp.F90`, B-grid subcycle** -- the single `do ksub = 1,ndte`
loop becomes a block/pass pair, with the exchange after the inner loop:

```fortran
ksub = 0
do kouter = 1, (ndte + nghost - 1)/nghost
do kinner = 1, min(nghost, ndte - ksub)
   ksub = ksub + 1
   ... stress, stepu ...
enddo
   ... halo exchange, ice_restoring_halo ...
enddo
```

The kernels need no shrinking loop bounds. The outer rings go stale as a
block proceeds, but that staleness never reaches the physical domain within
`nghost` passes, and the next exchange overwrites it.

### `ndte` not a multiple of `nghost`

Nothing special is required. **A short block is always safe**: running fewer
passes between exchanges consumes *less* halo than is available. Only a block
*longer* than `nghost` is wrong. So `min(nghost, ndte - ksub)` is the entire
remainder handling. With `ndte=137, nghost=2` that is 68 blocks of two passes
plus one block of one pass -- 137 passes, 69 exchanges.

### `nghost=1` is a no-op by construction

`(ndte+0)/1 = ndte` outer blocks of `min(1,·) = 1` pass each, with an
exchange after every pass -- structurally the original loop. The widened
bounds evaluate to `jlo..jhi+1` and `jlo..jhi` exactly. This is verified, not
assumed (test A below).

## Known defect at `nghost>1`

The bit-for-bit argument assumes a halo exchange **copies a neighbour's
computed values**, so recomputing a ghost cell locally reproduces them. That
holds for a ghost facing another block. It **fails at a domain boundary**,
where the exchange *imposes a boundary condition* rather than copying:
`linear_extrap` extrapolates, `closed` zeroes, `cyclic` wraps from the far
edge. Recomputing those cells replaces the boundary condition with a
locally-derived value, and the next pass reads a ghost that no longer
satisfies it.

Evidence: all of `linear_extrap`, `closed` and `cyclic` configurations differ
alike, and adding a full extra ring of margin (`nghost=3` with two passes per
exchange) changes nothing -- the error is at the domain edge, not a depth
off-by-one.

**The fix** is to exclude domain-boundary ghost cells from the widened
`icellT`/`icellU` ranges, so redundant computation happens only where a
genuine neighbour exists. The block structure already carries what is needed
(`this_block%i_glob`, `j_glob` are <= 0 or wrapped at domain edges); no new
machinery is required.

## Test

Three legs, one template, one compiler. Script: `nghost_blocking_test.sh`.

| leg | binary | `nghost` | expectation |
|-----|--------|----------|-------------|
| A | blocked | 1 | **fully bit-for-bit** vs unmodified upstream |
| B | blocked | 2 | the change under test |
| C | unmodified | 2 | control: isolates `nghost` from blocking |

Leg C matters. `nghost` changes `nx_block`/`ny_block`, so the restart file
changes size (`22x22` -> `24x24` blocks, 6.40 MB -> 6.72 MB) and will always
differ across `nghost` values. **Compare diagnostics and history across
`nghost`, never the restart.** Leg C establishes that `nghost=2` alone is
physics-neutral, so any physics difference in leg B is attributable to the
blocking.

### Result as of this commit

```
A  blocked nghost=1  vs upstream  ->  BIT-FOR-BIT IDENTICAL (6 binary files)
C  upstream nghost=2 vs upstream  ->  diagnostics + history IDENTICAL
                                      (restart differs: block shape only)
B  blocked nghost=2  vs upstream  ->  history DIFFERS
```

So the loop structure, the remainder handling and the mask widening are all
sound; the multi-pass is what is not yet correct, for the reason above.

A secondary, benign observation from an isolation run (widened masks but one
pass per exchange): diagnostics and history are identical, but the restart
differs deep in the file. The restart stores the stress arrays *including*
ghost cells, and widening the masks makes `stress` compute at a few more
ghosts. No physical field is affected.
