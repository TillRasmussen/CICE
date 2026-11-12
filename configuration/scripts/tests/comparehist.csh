#!/bin/csh -f

# Compare history data between two directories.
#
#   comparehist.csh <base_dir> <test_dir>
#
# No test has ever compared a history file against a baseline.  The suite
# compares iced restart files (comparerest.csh) and the diagnostic log
# (comparelog.csh), so a divergence confined to history output passes every
# row.  That is a category of output the framework does not inspect, not a
# gap in which configurations run.
#
# Only binary history can be compared byte for byte.  CICE stamps every
# history file with "This dataset was created on <date> at <time>", so two
# identical runs never produce identical files -- except with ICE_IOTYPE
# binary, where the stamp goes to the .hdr and the data to a separate .da.
# Restart files carry no stamp, which is why comparerest.csh can use cmp.
#
# So this compares iceh*.da and nothing else.  With netcdf or pio there are
# no .da files, the script reports that it compared nothing and exits 3, and
# the caller stays quiet.  Covering those needs cprnc or nccmp, which the
# test framework does not depend on today.

set base_dir = $1
set test_dir = $2

if ( $#argv != 2 ) then
  echo "Usage: ${0} <base_dir> <test_dir>"
  echo "   compares iceh*.da history data in the two directories"
  exit 99
endif

set failure = 0
set ncompared = 0

set nonomatch
foreach test_data (${test_dir}/iceh*.da)
  if ( ! -e ${test_data} ) continue
  set test_file = "${test_data:t}"
  set base_data = ${base_dir}/${test_file}
  echo "base_data: ${base_data}"
  echo "test_data: ${test_data}"
  if ( -e ${base_data} ) then
    @ ncompared = ${ncompared} + 1
    if ( { cmp -s ${base_data} ${test_data} } ) then
      echo "  compare OK"
    else
      set failure = 1
      echo "  compare FAIL"
    endif
  else
    if (${failure} == 0) set failure = 2
    echo "  missing data"
  endif
end
unset nonomatch

if (${ncompared} == 0 && ${failure} == 0) then
  echo "  no binary history data to compare"
  exit 3
endif

exit ${failure}
