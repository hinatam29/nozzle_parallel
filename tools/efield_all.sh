#!/bin/bash
# ==== すべての計算フォルダの電場の図（|E| と電気力線）を描き直す（2026-11-16）====
# 使い方:  bash tools/efield_all.sh          … 対象の一覧を表示するだけ
#          nohup bash tools/efield_all.sh run > tools/efield_all.log 2>&1 &   … 実行
# - viz_out/efield_field.png（いつもの高さの範囲）を作り直す
# - 並列（2_parallel、直下の parallel_pitch*）は viz_efield_wide/efield_field.png（高さ 70〜100 mm）も作る
# - flow.dat がフォルダにないもの（10/29 に D:\kawatani_data へ移したもの）は D: のものを使う
# - 計算中のフォルダ（go.out.log が 15 分以内に更新）は除外
cd "$(dirname "$0")/.."
mode=${1:-list}
now=$(date +%s)
list=()
for d in 1_single/*/* 2_parallel/*/* 3_opposed/*/* parallel_pitch*; do
  [ -f "$d/visualize.py" ] || continue
  if [ -f "$d/go.out.log" ] && [ $(( (now - $(date -r "$d/go.out.log" +%s)) / 60 )) -le 15 ]; then
    echo "  （計算中なので除外）$d"; continue
  fi
  if [ -f "$d/flow.dat" ]; then fl="flow.dat"
  elif [ -f "/mnt/d/kawatani_data/$d/flow.dat" ]; then fl="/mnt/d/kawatani_data/$d/flow.dat"
  else echo "  （flow.dat がない）$d"; continue; fi
  list+=("$d|$fl"); echo "  $d   [$fl]"
done
echo "対象 ${#list[@]} 個"
[ "$mode" = run ] || { echo "（一覧だけ。実行は: nohup bash tools/efield_all.sh run > tools/efield_all.log 2>&1 &）"; exit 0; }
i=0
for e in "${list[@]}"; do
  d=${e%%|*}; fl=${e##*|}; i=$((i+1))
  echo "=== [$i/${#list[@]}] $d  $(date +%H:%M)"
  ( cd "$d" || exit
    python3 visualize.py --kind efield --flow "$fl" 2>&1 | grep -i "保存\|描けません\|error\|スキップ"
    case $d in 2_parallel/*|parallel_pitch*)
      python3 visualize.py --kind efield --flow "$fl" --zmin 70 --zmax 100 --outdir viz_efield_wide 2>&1 | grep -i "保存\|描けません\|error\|スキップ";;
    esac )
done
echo "すべて完了 $(date +%H:%M)"
