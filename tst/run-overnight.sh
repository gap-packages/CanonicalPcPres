#!/usr/bin/env bash
#
# Launch CanonicalPcPres overnight stress tests, one GAP process per shard.
#
# Edit the settings below before running on a server, or override them with
# environment variables on the command line.
#

set -u

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PACKAGE_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
TIMESTAMP=$(date +%Y%m%d-%H%M%S)

GAP_BIN=${GAP_BIN:-gap}
SHARDS=${SHARDS:-127}
START=${START:-1}
FINISH=${FINISH:-2000}
CODE_REPEATS=${CODE_REPEATS:-100}
ISO_REPEATS=${ISO_REPEATS:-5}
INCLUDE_NILPOTENT=${INCLUDE_NILPOTENT:-false}
SKIP_MULTIPLES_OF_64=${SKIP_MULTIPLES_OF_64:-true}
LIMIT=${LIMIT:-0}
BASE_SEED=${BASE_SEED:-$(date +%s)}
LOG_ROOT=${LOG_ROOT:-"$PACKAGE_DIR/tst/overnight-logs/$TIMESTAMP"}
PROGRESS_ROOT=${PROGRESS_ROOT:-"$LOG_ROOT/progress"}

mkdir -p "$LOG_ROOT" "$PROGRESS_ROOT/groups" "$PROGRESS_ROOT/current" || exit 1

echo "CanonicalPcPres overnight run"
echo "  package: $PACKAGE_DIR"
echo "  logs: $LOG_ROOT"
echo "  progress: $PROGRESS_ROOT"
echo "  gap: $GAP_BIN"
echo "  shards: $SHARDS"
echo "  orders: $START..$FINISH"
echo "  code repeats: $CODE_REPEATS"
echo "  iso repeats: $ISO_REPEATS"
echo "  include nilpotent groups: $INCLUDE_NILPOTENT"
echo "  skip multiples of 64: $SKIP_MULTIPLES_OF_64"
echo "  limit per shard: $LIMIT"
echo "  base seed: $BASE_SEED"

declare -a PIDS
declare -a LOGS

for shard in $(seq 1 "$SHARDS"); do
  log="$LOG_ROOT/gap-$shard.log"
  seed=$((BASE_SEED + shard))
  LOGS[$shard]="$log"
  printf 'started\t%s\n' "$seed" >"$PROGRESS_ROOT/shard-$shard.started"

  (
    CANFORM_OVERNIGHT_START="$START" \
    CANFORM_OVERNIGHT_FINISH="$FINISH" \
    CANFORM_OVERNIGHT_CODE_REPEATS="$CODE_REPEATS" \
    CANFORM_OVERNIGHT_ISO_REPEATS="$ISO_REPEATS" \
    CANFORM_OVERNIGHT_SEED="$seed" \
    CANFORM_OVERNIGHT_INCLUDE_NILPOTENT="$INCLUDE_NILPOTENT" \
    CANFORM_OVERNIGHT_SKIP_MULTIPLES_OF_64="$SKIP_MULTIPLES_OF_64" \
    CANFORM_OVERNIGHT_SHARDS="$SHARDS" \
    CANFORM_OVERNIGHT_SHARD="$shard" \
    CANFORM_OVERNIGHT_LIMIT="$LIMIT" \
    CANFORM_OVERNIGHT_PROGRESS_DIR="$PROGRESS_ROOT" \
      "$GAP_BIN" -q --quitonbreak \
      --packagedirs "$PACKAGE_DIR" \
      "$PACKAGE_DIR/tst/overnight.g"
  ) >"$log" 2>&1 &

  PIDS[$shard]=$!
  echo "started GAP shard $shard out of $SHARDS, log: $log"
done

status=0

for shard in $(seq 1 "$SHARDS"); do
  pid=${PIDS[$shard]}
  log=${LOGS[$shard]}
  if wait "$pid"; then
    printf 'done\n' >"$PROGRESS_ROOT/shard-$shard.done"
    echo "gap shard $shard out of $SHARDS completed, see $log"
  else
    rc=$?
    printf 'failed\t%s\n' "$rc" >"$PROGRESS_ROOT/shard-$shard.failed"
    echo "gap shard $shard out of $SHARDS failed with exit code $rc,"
    echo "  see $log"
    status=1
  fi
done

if [ "$status" -eq 0 ]; then
  echo "all GAP shards completed successfully"
else
  echo "one or more GAP shards failed"
fi

completed_groups=$(find "$PROGRESS_ROOT/groups" -type f -name 'grp-*.done' | wc -l | tr -d ' ')
echo "completed group markers: $completed_groups"
echo "progress directory: $PROGRESS_ROOT"

exit "$status"
