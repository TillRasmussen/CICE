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

**1. `ice_dyn_shared.F90`** -- `dyn_prep2` widens the `icellT` and `icellU`
bounds, and `init_dyn_shared` widens the `cyp/cxp/cym/cxm` initialisation,
both by `nwide = halo_blocking_width() - 1` rings in each direction:

```fortran
nwide = halo_blocking_width() - 1
do j = jlo-nwide, jhi+1+nwide          ! icellT, was jlo, jhi+1
do j = jlo-nwide, jhi+nwide            ! icellU, was jlo, jhi
```

⚠ Both must take `nwide` from the SAME function the subcycle uses for its
pass count. Widening and pass count disagreeing is a silent wrong answer.

⚠ The widening is currently INERT -- see "Known defect" below.

**2. `ice_dyn_evp.F90`, B-grid subcycle** -- the single `do ksub = 1,ndte`
loop becomes a block/pass pair, with the exchange after the inner loop:

```fortran
ksub = 0
nblk = halo_blocking_width()
do kouter = 1, (ndte + nblk - 1)/nblk
do kinner = 1, min(nblk, ndte - ksub)
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
*longer* than `nghost` is wrong. So `min(nblk, ndte - ksub)` is the entire
remainder handling. With `ndte=137, nghost=2` that is 68 blocks of two passes
plus one block of one pass -- 137 passes, 69 exchanges.

### `nghost=1` is a no-op by construction

With `nblk = 1`: `ndte` outer blocks of `min(1,·) = 1` pass each, with an
exchange after every pass -- structurally the original loop. The widened
bounds evaluate to `jlo..jhi+1` and `jlo..jhi` exactly. This is verified, not
assumed (test A below).

## Known defect at `nghost>1`: halo inputs are never computed

**Root cause: every input the kernels need is only valid in the physical
domain, so the extra rings cannot be advanced at all.**

`grid_average_X2Y` writes only `jlo..jhi, ilo..ihi`, and none of its outputs
-- `umass`, `aiU`, `uocnU`, `vocnU`, `strairxU` -- is halo-updated
afterwards. At any ghost cell `aiU` and `umass` are therefore zero, so
`dyn_prep2` computes `iceUmask = .false.` there, `icellU` excludes the cell,
and **`stepu` never runs at ring 1 regardless of how the loop bounds are
widened**. The widening is inert: a later pass in a block then reads a ring
that nobody advanced.

A second instance of the same class, fixed in this commit because it is a
genuine prerequisite: `init_dyn_shared` initialises `cyp`, `cxp`, `cym`,
`cxm` over a hardcoded `jlo..jhi+1` / `ilo..ihi+1`. `stress` reads them
wherever `icellT` reaches, so at ring 2 they were uninitialised. That fix is
necessary but not sufficient.

### Hypotheses tested and refuted

Recorded so they are not re-tried:

1. **Domain boundary conditions.** That recomputing a ghost cannot reproduce
   an imposed BC (`linear_extrap`, `closed`, `cyclic` wrap, tripole fold).
   Refuted: a **doubly cyclic** domain, where every ghost is backed by real
   data and the local array is physically contiguous, fails identically.
2. **Depth off-by-one.** Refuted: `nghost=3` with two passes per exchange --
   a full extra ring of margin -- fails at the same history byte.
3. **Uninitialised `cyp/cxp/cym/cxm`.** A real bug, fixed here; still fails.
4. **Unexchanged stress tensor.** `stressp_*`, `stressm_*`, `stress12_*` are
   persistent state the subcycle never exchanges, so a stale outer ring
   could corrupt them permanently. Refuted: exchanging all twelve once per
   block changes nothing -- because of the root cause above, the outer rings
   were never being advanced in the first place.

### What it would take

Temporal blocking needs **every** input to `stress` and `stepu` valid
`nghost` rings out, not one ring or zero:

- `grid_average_X2Y` to write the halo, or its outputs to be halo-updated;
- `iceTmask`/`iceUmask` consequently valid over the full halo;
- the per-timestep forcing terms (`waterx`, `forcex`, `umassdti`, ...)
  likewise, which follows once `dyn_prep2`'s inputs are valid;
- and the boundary-condition question of hypothesis 1 revisited *after*
  that, since it was never actually exercised.

That is a change to the whole per-timestep preparation, not to two loop
bounds. It is the reason this branch stops here.

### The gate

`halo_blocking_width()` currently returns `nghost` only for a doubly cyclic
domain and 1 otherwise, so every standard configuration keeps the historical
one-exchange-per-subcycle behaviour and is bit-for-bit. ⚠ The doubly-cyclic
condition was chosen under hypothesis 1 and is **not** the right criterion --
hypothesis 1 is refuted. It is retained only as a conservative switch that
keeps the branch harmless; revisit it once the root cause above is fixed.

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

### Safety: the default path is untouched

`nghost=1` makes `halo_blocking_width()` return 1, every widening evaluate
to its original bounds, and the subcycle degenerate to the original loop.
Verified against a clean `upstream/main` worktree -- not against a stashed
copy of this branch, which is a mistake worth avoiding:

```
                evpB      box2001   rmap      gx3 (tripole)
   gfortran     BFB       BFB       BFB       BFB
   ifx          BFB       BFB       BFB       BFB
```

Four boundary regimes (linear_extrap, closed, cyclic/open, cyclic/tripole),
both compilers, transport active in rmap and gx3.

⚠ Still untested: multi-task MPI.  Every run above is nprocs=1, where a
halo "exchange" is a local copy and no message is ever sent.

### Result as of this commit

```
A  blocked nghost=1  vs upstream  ->  BIT-FOR-BIT IDENTICAL (6 binary files)
C  upstream nghost=2 vs upstream  ->  diagnostics + history IDENTICAL
                                      (restart differs: block shape only)
B  blocked nghost=2  vs upstream  ->  history DIFFERS
```

So the loop structure and the remainder handling are sound, and `nghost=1`
is exact. The mask widening is *inert* for the reason above, and the
multi-pass is therefore not yet correct.

⚠ Note leg B PASSES on any configuration the gate disables (anything not
doubly cyclic) -- that is the gate working, not the blocking working. Run
the cyclic case to exercise it.

A secondary, benign observation from an isolation run (widened masks but one
pass per exchange): diagnostics and history are identical, but the restart
differs deep in the file. The restart stores the stress arrays *including*
ghost cells, and widening the masks makes `stress` compute at a few more
ghosts. No physical field is affected.
