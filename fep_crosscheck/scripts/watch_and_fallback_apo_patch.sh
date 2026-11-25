#!/bin/bash
# Watches apo_patch's FFIX production array (job 7575173, on ada5000). If it
# hasn't started running within 3 hours, cancels the whole apo_patch FFIX
# chain (prod 7575173 / energy 7575174 / mbar 7575175) and resubmits it fresh
# on ada6000 instead.
set -u
PROD_JOBID=7575173
ENERGY_JOBID=7575174
MBAR_JOBID=7575175
SYSDIR=${SYSDIR:?set SYSDIR to your system root}
SCRIPTDIR=${SCRIPTDIR:?set SCRIPTDIR to this repo's scripts/ dir}
DEADLINE=$(( $(date +%s) + 3*3600 ))

while [ $(date +%s) -lt $DEADLINE ]; do
    if squeue -j ${PROD_JOBID} -h -o "%T" 2>/dev/null | grep -q RUNNING; then
        echo "[$(date)] apo_patch production (job ${PROD_JOBID}) has started on ada5000 -- no action needed."
        exit 0
    fi
    if ! squeue -j ${PROD_JOBID} -h 2>/dev/null | grep -q .; then
        echo "[$(date)] job ${PROD_JOBID} is no longer in the queue (completed/failed) -- no action needed."
        exit 0
    fi
    sleep 3600
done

echo "[$(date)] 3 hours elapsed, job ${PROD_JOBID} still not running on ada5000 -- falling back to ada6000."
scancel ${PROD_JOBID} ${ENERGY_JOBID} ${MBAR_JOBID}
sleep 3
rm -rf "${SYSDIR}/ffix/dcd"/* "${SYSDIR}/ffix/res"/* "${SYSDIR}/ffix/output"/* "${SYSDIR}/ffix/out"/*
bash "${SCRIPTDIR}/submit_ffix_chain.sh" "${SYSDIR}" ada6000 20 5 150000 5000000 5000
echo "[$(date)] Resubmitted apo_patch FFIX chain on ada6000."
