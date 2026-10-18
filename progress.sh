#!/bin/bash
# ==== 計算の進み具合をまとめて表示 ====
# 使い方:  bash progress.sh            … 1回だけ表示
#          watch -n 60 bash progress.sh … 60秒ごとに自動更新（Ctrl+C で終了）
# 対象: 1_single/*/*、2_parallel/*/pitch*、3_opposed/*/*、直下に残っている parallel_pitch*
cd "$(dirname "$0")"
now=$(date +%s)
printf "%-58s %8s %7s %7s %6s  %s\n" "フォルダ" "ステップ" "進捗" "液滴数" "BADP" "状態（ログ更新）"
for d in 1_single/*/* 2_parallel/*/pitch* 3_opposed/*/* parallel_pitch*; do
  log=$d/go.out.log
  [ -f "$log" ] || continue
  # 最終ステップ数 = ftime/dt（main.f のコメントでない行から読む）
  iend=$(awk '/^ +dt *=/{gsub(/[^0-9.eE+-]/,"",$3);dt=$3} /^ +ftime *=/{gsub(/[^0-9.eE+-]/,"",$3);ft=$3} END{if(dt>0) printf "%d", ft/dt+0.5}' $d/main.f)
  last=$(grep -E '^ +[0-9]+ ' "$log" | tail -n 1)
  stp=$(echo "$last" | awk '{print $1}'); np=$(echo "$last" | awk '{print $3}')
  [ -z "$stp" ] && stp=0
  pct=$(awk -v s="$stp" -v e="$iend" 'BEGIN{ if(e>0) printf "%.1f%%", 100*s/e; else print "-"}')
  age=$(( (now - $(date -r "$log" +%s)) / 60 ))
  if [ "$stp" -ge "$iend" ] 2>/dev/null; then st="完了"
  elif [ $age -le 15 ]; then st="計算中（${age}分前）"
  else st="止まっている？（${age}分前）"; fi
  printf "%-58s %8s %7s %7s %6s  %s\n" "$d" "$stp" "$pct" "$np" "$(grep -c BADP $log)" "$st"
done
echo "現在時刻: $(date +%m/%d\ %H:%M)    動いている go.out の数: $(pgrep -c -x go.out)"
