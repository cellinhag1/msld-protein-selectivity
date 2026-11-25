#!/bin/bash
# Orchestrates the full FFIX pipeline for ONE system, entirely on ONE
# partition (this cluster's sbatch does not accept a combined
# "-p ada5000,ada6000" request, so assign whole systems to whichever
# partition you want them on):
#   1. FFIX production (N_WINDOWS-window array)
#   2. energy-matrix evaluation (N_WINDOWS x CHUNKS array), dependency on 1
#   3. gather + FastMBAR, dependency on 2
#
# Usage: bash submit_ffix_chain.sh <sysdir> <partition> [N_WINDOWS] [CHUNKS] [NEQUIL] [NPROD] [NSAVC]
#   sysdir    - absolute path to apo_patch / apo_6site / apo_7site
#   partition - e.g. ada5000 or ada6000
#   N_WINDOWS - default 20
#   CHUNKS    - default 5 (frame-chunks per window for the energy pass)
#   NEQUIL    - default 150000 steps (300 ps @ dt=0.002)
#   NPROD     - default 5000000 steps (10 ns @ dt=0.002)
#   NSAVC     - default 5000 steps/frame (10 ps/frame -> 1000 frames/window)

set -euo pipefail

SYSDIR="$1"
PART="$2"
N_WINDOWS="${3:-20}"
CHUNKS="${4:-5}"
NEQUIL="${5:-150000}"
NPROD="${6:-5000000}"
NSAVC="${7:-5000}"

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LASTW=$(( N_WINDOWS - 1 ))
LASTT=$(( N_WINDOWS * CHUNKS - 1 ))

cd "$SYSDIR"
mkdir -p ffix/out ffix/res ffix/dcd ffix/output

echo "== [$SYSDIR] 1/3 FFIX production: ${N_WINDOWS} windows on ${PART} =="
PROD_JOBID=$(sbatch --parsable --partition=${PART} --array=0-${LASTW} "${SCRIPTDIR}/run_ffix.sb" ${NEQUIL} ${NPROD} ${NSAVC})
echo "   prod jobid: $PROD_JOBID"

echo "== [$SYSDIR] 2/3 energy matrix: ${N_WINDOWS} windows x ${CHUNKS} chunks on ${PART} =="
ENERGY_JOBID=$(sbatch --parsable --partition=${PART} --array=0-${LASTT} --dependency=afterok:${PROD_JOBID} \
  "${SCRIPTDIR}/compute_energies_ffix.sb" ${N_WINDOWS} ${CHUNKS})
echo "   energy jobid: $ENERGY_JOBID"

echo "== [$SYSDIR] 3/3 gather + FastMBAR =="
MBAR_JOBID=$(sbatch --parsable --dependency=afterok:${ENERGY_JOBID} \
  "${SCRIPTDIR}/gather_and_mbar.sb" ${N_WINDOWS} ${CHUNKS})
echo "   mbar jobid: $MBAR_JOBID"

echo "Done. Final result: ${SYSDIR}/ffix/output/mbar_results.npz"
