#!/bin/bash
# ==== go.out 起動スクリプト（OpenMP版）====
# 使い方:  bash run.sh          … 4コアで実行
#          bash run.sh 8        … 8コアで実行
#          bash run.sh 1        … 1コア（従来と同じ動作の確認用）

# 大きな作業配列をスタックに置くため、スタック上限を撤廃（segfault対策）
ulimit -s unlimited 2>/dev/null || ulimit -s 1048576
export OMP_STACKSIZE=512M

# 使用コア数（引数で指定。無指定なら4）
export OMP_NUM_THREADS=${1:-4}

echo "OMP_NUM_THREADS=$OMP_NUM_THREADS で実行します"
./go.out
