#!/bin/bash
# Three-leg test for halo-blocked EVP subcycling.  See NGHOST_BLOCKING.md.
#
#   nghost_blocking_test.sh <builddir-blocked> <builddir-upstream> <runtemplate>
#
#     A  blocked   nghost=1   must be BIT-FOR-BIT with unmodified upstream
#     B  blocked   nghost=2   the change under test
#     C  upstream  nghost=2   control: isolates nghost from the blocking
#
# ⚠ nghost changes nx_block/ny_block, so the restart file changes SIZE and
#   will always differ across nghost values.  Compare diagnostics and history
#   across nghost; never the restart.  Leg C exists to prove that nghost=2 on
#   its own is physics-neutral, so that any physics difference in leg B is
#   attributable to the blocking and not to the wider halo.
set -eu

BLK=${1:?directory holding the blocked cice binary}
REF=${2:?directory holding the unmodified upstream cice binary}
TPL=${3:?run template directory to copy}
OUT=${OUT:-$PWD/nghost_test}

mkdir -p "$OUT"

run_leg () {   # name binary nghost
   local name=$1 bin=$2 ng=$3 d="$OUT/$1"
   rm -rf "$d"; cp -r "$TPL" "$d"
   # domain_nml may not carry an nghost line; substitute if present, else insert
   if grep -q '^ *nghost *=' "$d/ice_in"; then
      sed -i "s/^\( *nghost *= *\).*/\1$ng/" "$d/ice_in"
   else
      sed -i "s/^\( *max_blocks *=.*\)/\1\n    nghost            = $ng/" "$d/ice_in"
   fi
   cp "$bin" "$d/cice"
   ( cd "$d" && ulimit -s unlimited 2>/dev/null; ./cice > out.log 2>&1 ) || true
   # ⚠ confirm the namelist actually took effect -- a silent sed miss here
   #   makes every leg run at nghost=1 and the test passes meaninglessly
   local got
   got=$(grep -m1 -oE 'nx_block *=? *[0-9]+' "$d/out.log" | grep -oE '[0-9]+' || echo "?")
   printf '  %-10s nghost=%s  nx_block=%s\n' "$name" "$ng" "$got"
}

echo "legs:"
run_leg ref  "$REF/cice" 1
run_leg blk1 "$BLK/cice" 1
run_leg blk2 "$BLK/cice" 2
run_leg ref2 "$REF/cice" 2

cmp_leg () {  # label a b  -- report diagnostics/history, note restart separately
   echo
   echo "######## $1 ########"
   local a="$OUT/$2" b="$OUT/$3" rc=0
   # ⚠ Do NOT raw-diff the log: it carries timestamps and paths and will
   #   always differ.  Filter to the physical diagnostics, as the suite does.
   local fa=$OUT/.diag_a fb=$OUT/.diag_b
   for pair in "$a:$fa" "$b:$fb"; do
      local d=${pair%%:*} o=${pair##*:}
      grep -A 99999999 "total ice area  (km^2)" "$d/out.log" \
        | grep -e istep1: -e = | grep -iv "min, max, sum" \
        | grep -iv "init_vert" | grep -iv "ridge_ice" > "$o" || true
   done
   if [ ! -s "$fa" ] || [ ! -s "$fb" ]; then
      echo "  diagnostics: INCONCLUSIVE (no output)"; rc=1
   elif diff -w -q "$fa" "$fb" >/dev/null; then
      echo "  diagnostics: IDENTICAL"
   else
      echo "  diagnostics: DIFFER"; rc=1
   fi
   for f in "$a"/history/*.da; do
      [ -e "$f" ] || continue
      local n; n=$(basename "$f")
      if cmp -s "$f" "$b/history/$n"; then echo "  history $n: IDENTICAL"
      else echo "  history $n: DIFFER"; rc=1; fi
   done
   for f in "$a"/restart/iced*; do
      [ -e "$f" ] || continue
      local n; n=$(basename "$f")
      if cmp -s "$f" "$b/restart/$n"; then echo "  restart $n: identical"
      else echo "  restart $n: differs (expected across nghost: block shape)"; fi
   done
   return $rc
}

cmp_leg "A  blocked nghost=1 vs upstream   (MUST be bit-for-bit)" ref blk1 \
   && echo "  => A PASS" || echo "  => A FAIL"
cmp_leg "C  upstream nghost=2 vs upstream  (control, physics-neutral)" ref ref2 \
   && echo "  => C PASS" || echo "  => C FAIL"
cmp_leg "B  blocked nghost=2 vs upstream   (change under test)" ref blk2 \
   && echo "  => B PASS" || echo "  => B FAIL (known: domain-boundary ghosts)"
