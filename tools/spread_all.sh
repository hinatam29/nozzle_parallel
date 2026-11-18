#!/bin/bash
# ==== すべての計算フォルダの液滴の広がりの図を作り直す（2026-11-18、描く順番をランダムにした後）====
# 使い方:  bash tools/spread_all.sh          … 対象の一覧を表示するだけ
#          nohup bash tools/spread_all.sh run > tools/spread_all.log 2>&1 &   … 実行
# - viz_out/spread_*（いつもの高さの範囲、静止画と GIF）を作り直す
# - 並列（2_parallel、直下の parallel_pitch*）は viz_out_full/spread_*（高さ 0〜100 mm）も作り直す
# - 22 ms まで計算したもの（ftime = 22.e-3）は静止画を 12 枚（2 ms ごと）
# - 計算中のフォルダ（go.out.log が 15 分以内に更新）は除外
cd "$(dirname "$0")/.."
mode=${1:-list}
now=$(date +%s)
list=()
for d in 1_single/*/* 2_parallel/*/* 3_opposed/*/* parallel_pitch*; do
  [ -f "$d/visualize.py" ] && [ -f "$d/spray.dat" ] || continue
  if [ -f "$d/go.out.log" ] && [ $(( (now - $(date -r "$d/go.out.log" +%s)) / 60 )) -le 15 ]; then
    echo "  （計算中なので除外）$d"; continue
  fi
  list+=("$d"); echo "  $d"
done
echo "対象 ${#list[@]} 個"
[ "$mode" = run ] || { echo "（一覧だけ。実行は: nohup bash tools/spread_all.sh run > tools/spread_all.log 2>&1 &）"; exit 0; }
i=0
for d in "${list[@]}"; do
  i=$((i+1)); echo "=== [$i/${#list[@]}] $d  $(date +%H:%M)"
  ( cd "$d" || exit
    ns=""; grep -q "^ *ftime = 22" main.f 2>/dev/null && ns="--nsnap 12"
    python3 visualize.py --kind spread --gif $ns 2>&1 | grep -i "保存.*gif\|error\|Traceback"
    case $d in 2_parallel/*|parallel_pitch*)
      python3 visualize.py --kind spread --zmin 0 --zmax 100 --outdir viz_out_full --gif $ns 2>&1 | grep -i "保存.*gif\|error\|Traceback";;
    esac )
done
echo "すべて完了 $(date +%H:%M)"
