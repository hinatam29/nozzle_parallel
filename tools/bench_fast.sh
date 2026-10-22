#!/bin/bash
# 高速化の効果を、このパソコンで測る（2026-10-22、改訂版）
# 使い方: kawatani フォルダで  bash tools/bench_fast.sh
# 並列 fieldFIX2 11 mm の最後（液滴 約 1.6 万個）から、
#   1 ステップだけ（＝最初の電位の計算の時間）と 101 ステップを計算し、
#   その差から「1 ステップあたりの時間」を出す。
# 作業用のフォルダ _bench を作る（終わったら消してよい）。
set -e
K=$(cd "$(dirname "$0")/.." && pwd)
SRC=$K/2_parallel/fieldFIX2/pitch11mm_hole5mm_d11nm_k1e-4
B=$K/_bench
mkdir -p $B
for v in orig fast; do
  for n in 1 101; do
    d=$B/${v}_$n
    mkdir -p $d
    cp $SRC/*.f $SRC/subcom.inc $SRC/Makefile $SRC/restart.dat $d/
    if [ $n = 1 ]; then ft=11.0e-3; else ft=11.02e-3; fi
    sed -i "s/^      ftime = 11.e-3/      ftime = $ft/" $d/main.f
    [ $v = fast ] && python3 $K/tools/patch_fast.py $d > /dev/null
    (cd $d && rm -f *.o go.out && make -s > /dev/null 2>&1)
  done
done
ulimit -s unlimited 2>/dev/null || true
export OMP_STACKSIZE=512M
echo "（全部で 5〜15 分かかります）"
for t in 8 12 24; do
  export OMP_NUM_THREADS=$t
  for v in orig fast; do
    for n in 1 101; do
      cd $B/${v}_$n
      s=$(date +%s.%N); ./go.out > run_$t.log 2>&1; e=$(date +%s.%N)
      eval "T_${v}_$n=$(awk -v a=$s -v b=$e 'BEGIN{printf "%.2f", b-a}')"
    done
    a=$(eval echo \$T_${v}_1); b=$(eval echo \$T_${v}_101)
    echo "$v ${t}コア: 最初の電位 ${a} 秒, 1ステップ $(awk -v a=$a -v b=$b 'BEGIN{printf "%.3f", (b-a)/100}') 秒"
  done
done
echo "--- 結果が同じかの確認（101 ステップ後の ctotal、24コア）---"
tail -n 1 $B/orig_101/run_24.log
tail -n 1 $B/fast_101/run_24.log
