#!/bin/bash
# ==== 計算フォルダの名前の整理（2026-11-10）====
# 方針: フォルダ名の最後に「電圧_計算時間」をつける（例 _7kV_11ms）。
#       液滴をランダムに出す版は multi5r → multi5random にまとめる（3kV・22ms 版も同じ場所へ）。
#       ランダムにせず 3 本の筋になった失敗版は fieldFIX3_multi5_notrandom_failed にする。
# 使い方: bash tools/rename_folders.sh        … 変更の一覧を表示するだけ（何も変えない）
#         bash tools/rename_folders.sh run    … 実行する（計算中のフォルダは含まない）
#         bash tools/rename_folders.sh later  … 計算中だった 32mm 22ms の分を、終わってから実行する
cd "$(dirname "$0")/.."
MAIN=(
"1_single/fieldFIX/single_7kV_d11nm_k1e-4|1_single/fieldFIX/single_d11nm_k1e-4_7kV_11ms"
"1_single/fieldFIX/single_7kV_ring_d11nm_k1e-4_fine|1_single/fieldFIX/single_ring_d11nm_k1e-4_7kV_11ms_fine"
"1_single/fieldFIX/single_7kV_ring_d11nm_k1e-4|1_single/fieldFIX/single_ring_d11nm_k1e-4_7kV_11ms"
"1_single/fieldOLD/repro_6kV|1_single/fieldOLD/repro_d500nm_6kV_11ms"
"1_single/fieldOLD/repro_7kV|1_single/fieldOLD/repro_d500nm_7kV_11ms"
"1_single/onefactor/repro_7kV_fieldFIX_d11nm|1_single/onefactor/repro_fieldFIX_d11nm_7kV_11ms"
"1_single/onefactor/repro_7kV_fieldFIX|1_single/onefactor/repro_fieldFIX_d500nm_7kV_11ms"
"1_single/xy/single_7kV_d11nm_k1e-4_xy|1_single/xy/single_xy_d11nm_k1e-4_7kV_11ms"
"parallel_pitch11mm_hole5mm_7kV_q100|2_parallel/fieldOLD/pitch11mm_hole5mm_d500nm_q100_7kV_11ms"
"parallel_pitch5mm_7kV_q100|2_parallel/fieldOLD/pitch5mm_d500nm_q100_7kV_11ms"
"3_opposed/fieldFIX/opp13mm_ring_7kV_d11nm_k1e-4_fine|3_opposed/fieldFIX/opp13mm_ring_d11nm_k1e-4_7kV_11ms_fine"
"3_opposed/fieldFIX/opp13mm_ring_7kV_d11nm_k1e-4|3_opposed/fieldFIX/opp13mm_ring_d11nm_k1e-4_7kV_11ms"
"3_opposed/fieldFIX/opp33mm_ring_7kV_d11nm_k1e-4_fine|3_opposed/fieldFIX/opp33mm_ring_d11nm_k1e-4_7kV_11ms_fine"
"3_opposed/fieldFIX/opp33mm_ring_7kV_d11nm_k1e-4|3_opposed/fieldFIX/opp33mm_ring_d11nm_k1e-4_7kV_11ms"
)
for g in fieldOLD fieldFIX fieldFIX2; do
  for d in 2_parallel/$g/pitch*; do
    case $d in *_7kV_11ms) continue;; esac
    [ -d "$d" ] && MAIN+=("$d|${d}_7kV_11ms")
  done
done
for d in 2_parallel/fieldFIX3_multi5/pitch*; do
  [ -d "$d" ] && MAIN+=("$d|2_parallel/fieldFIX3_multi5_notrandom_failed/$(basename $d)_7kV_11ms")
done
for d in 2_parallel/fieldFIX3_multi5r/pitch*; do
  [ -d "$d" ] && MAIN+=("$d|2_parallel/fieldFIX3_multi5random/$(basename $d)_7kV_11ms")
done
for d in 2_parallel/fieldFIX3_multi5r_3kV/pitch*; do
  b=$(basename $d); b=${b%_3kV}
  [ -d "$d" ] && MAIN+=("$d|2_parallel/fieldFIX3_multi5random/${b}_3kV_11ms")
done
MAIN+=("2_parallel/fieldFIX3_multi5r_long/pitch11mm_hole5mm_d11nm_k1e-4_q100_22ms|2_parallel/fieldFIX3_multi5random/pitch11mm_hole5mm_d11nm_k1e-4_q100_7kV_22ms")
LATER=("2_parallel/fieldFIX3_multi5r_long/pitch32mm_hole5mm_d11nm_k1e-4_q100_22ms|2_parallel/fieldFIX3_multi5random/pitch32mm_hole5mm_d11nm_k1e-4_q100_7kV_22ms")

mode=${1:-list}
if [ "$mode" = later ]; then LIST=("${LATER[@]}"); else LIST=("${MAIN[@]}"); fi
n=0; bad=0
for e in "${LIST[@]}"; do
  s=${e%%|*}; t=${e##*|}
  if [ ! -d "$s" ]; then echo "  （なし）$s"; continue; fi
  if [ -e "$t" ]; then echo "  ×すでにある: $t"; bad=1; continue; fi
  printf "  %s\n    → %s\n" "$s" "$t"; n=$((n+1))
done
echo "対象 $n 個"
[ "$mode" = list ] && { echo "（一覧だけ。実行は: bash tools/rename_folders.sh run）"; echo "後で（32mm 22ms が終わってから）: ${LATER[0]}"; exit 0; }
[ $bad = 1 ] && { echo "同じ名前のフォルダがあるので中止しました"; exit 1; }
for e in "${LIST[@]}"; do
  s=${e%%|*}; t=${e##*|}
  [ -d "$s" ] || continue
  mkdir -p "$(dirname "$t")"
  mv "$s" "$t" || { echo "移動できませんでした: $s（エクスプローラーや別のタブで開いていないか確認）"; exit 1; }
  # README の中の古い名前も書き換える（後ろに英数字が続かない場合だけ）
  sed -i -E "s#${s}([^A-Za-z0-9_]|\$)#${t}\1#g" README.md
done
for g in fieldFIX3_multi5 fieldFIX3_multi5r fieldFIX3_multi5r_3kV fieldFIX3_multi5r_long; do
  rmdir "2_parallel/$g" 2>/dev/null && echo "空になったので削除: 2_parallel/$g"
done
echo "完了"
