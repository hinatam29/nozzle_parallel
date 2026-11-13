#!/bin/bash
# ==== すべての計算フォルダで、二酸化炭素吸収量の図を川谷さんと同じ気相の C_total で描き直す（2026-11-11）====
# 使い方:  bash tools/viz_gas_all.sh          … 対象の一覧を表示するだけ
#          nohup bash tools/viz_gas_all.sh run > tools/viz_gas_all.log 2>&1 &   … 実行（数十分）
# - 各フォルダに viz_gas/（86〜100 mm、対向は 30〜70 mm）を作る。並列は viz_gas_full/（0〜100 mm）も
# - flow.dat がフォルダにないもの（10/29 に D:\kawatani_data へ移したもの）は D: のものを使う
# - 以前の液滴側の吸収の図（viz_out/absorb_*）は viz_out/old_droplet_absorb/ に移す（消さない）
# - 除外: 3本の筋になった失敗版（fieldFIX3_multi5）、計算中のフォルダ（go.out.log が 15 分以内に更新）
cd "$(dirname "$0")/.."
ROOT=$(pwd)
mode=${1:-list}
now=$(date +%s)
list=()
for d in 1_single/*/* 2_parallel/*/* 3_opposed/*/* parallel_pitch*; do
  [ -f "$d/main.f" ] || continue
  case $d in 2_parallel/fieldFIX3_multi5/*|2_parallel/fieldFIX3_multi5_notrandom_failed/*) continue;; esac
  if [ -f "$d/go.out.log" ] && [ $(( (now - $(date -r "$d/go.out.log" +%s)) / 60 )) -le 15 ]; then
    echo "  （計算中なので除外）$d"; continue
  fi
  if [ -f "$d/flow.dat" ]; then src="flow.dat"
  elif [ -f "/mnt/d/kawatani_data/$d/flow.dat" ]; then src="D:"
  else echo "  （flow.dat がない）$d"; continue; fi
  list+=("$d"); echo "  $d   [$src]"
done
echo "対象 ${#list[@]} 個"
[ "$mode" = run ] || { echo "（一覧だけ。実行は: nohup bash tools/viz_gas_all.sh run > tools/viz_gas_all.log 2>&1 &）"; exit 0; }
i=0
for d in "${list[@]}"; do
  i=$((i+1)); echo "=== [$i/${#list[@]}] $d  $(date +%H:%M)"
  ( cd "$d" || exit
    if ls viz_out/absorb_* >/dev/null 2>&1; then
      mkdir -p viz_out/old_droplet_absorb && mv viz_out/absorb_* viz_out/old_droplet_absorb/
    fi
    opt=""; case $d in 2_parallel/*|parallel_pitch*) opt="--full";; esac
    python3 "$ROOT/tools/viz_gas.py" --gif $opt 2>&1 | grep -v "^  ->\|^  flow.dat: step" )
done
echo "すべて完了 $(date +%H:%M)"
