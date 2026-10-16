# kawatani フォルダの案内

静電噴霧による二酸化炭素吸収シミュレーション（川谷さんのコードを元に改良）。
各計算フォルダは独立していて、その中で `make` → `nohup bash run.sh 8 > go.out.log 2>&1 &` で計算できる。
進み具合は、このフォルダで `bash progress.sh`。

## フォルダの分け方

    kawatani/
    ├─ README.md            … このファイル
    ├─ progress.sh          … 計算の進み具合をまとめて表示
    ├─ CHANGES.md           … 1本→並列にしたときの変更履歴
    ├─ docs/                … 説明スライドなど
    ├─ 0_original/          … 川谷さんの元のコード（1本ノズル）
    ├─ 1_single/            … 1本ノズル・円筒（軸対称）
    │   └─ fieldOLD/        … 電場計算に誤りあり（下記）
    ├─ 2_parallel/          … 並列ノズル（平面2次元・左右が鏡面）
    │   ├─ fieldOLD/        … 電場計算に誤りあり
    │   └─ fieldFIX/        … 電場計算を修正した版（2026-10-16〜）
    └─ 3_opposed/           … 対向ノズル（これから）

各計算フォルダの名前：`pitch〇〇mm`（ノズル間隔）`_hole〇mm`（対向電極の穴の半幅）`_d〇〇nm`（液滴の平均径）`_k〇〇`（吸収係数 m/s）。
特に書いていない条件は共通：ノズル +7 kV・対向電極 0 V、限界帯電量 100 %、初速度 −3.0 m/s、計算時間 11 ms（55,000 ステップ）。

## 計算の一覧

| フォルダ | 内容 |
|---|---|
| 0_original/kawatani_original_1nozzle | 川谷さんの元の6ファイル |
| 1_single/fieldOLD/repro_6kV, repro_7kV | 川谷さん条件の再現（円筒、1本、6 kV / 7 kV） |
| 2_parallel/fieldOLD/pitch16mm_hole5mm_d500nm | 並列 16 mm、液滴 500 nm、横格子 32 µm（もとは kawatani 直下にあった計算） |
| 　〃 の中の `*_3kv_oldpolarity.*` | 同じフォルダの古い計算：ノズル 0 V・対向電極 −3 kV（9/29）。図は `viz_out/ノズル0kV・対向電極-3kV` |
| 　〃 の中の `*_3kv_80pct.*` | ノズル 3 kV・対向電極 0 V・帯電量 80 %（10/2）。図は `viz_out/20261004ノズル3kV・対向電極0kV・帯電量80％` |
| 　〃 の中の `*_7kv_pcountfix.*`、`*_bugged.*`、`*_old.*` | 9/24〜9/28 のバグ修正前後の試し計算 |
| 2_parallel/fieldOLD/pitch05mm_hole1.56mm_d500nm | 並列 5 mm、穴を比例縮小（開口率 62.5 %） |
| 2_parallel/fieldOLD/pitch11mm_hole5mm_d500nm | 並列 11 mm、液滴 500 nm |
| 2_parallel/fieldOLD/pitch32mm_hole5mm_d500nm | 並列 32 mm、液滴 500 nm |
| 2_parallel/fieldOLD/pitch{11,16,32}mm_hole5mm_d11nm_k1e-4 | 液滴径を実験に合わせた版（平均 11.4 nm、k = 1e-4 m/s、カニンガム補正、横格子 22 µm） |
| 2_parallel/fieldFIX/pitch{11,16,32}mm_hole5mm_d11nm_k1e-4 | 上と同じ条件で、電場計算を修正 |

※ 2026-10-16 時点で `parallel_pitch5mm_7kV_q100`・`parallel_pitch11mm_hole5mm_7kV_q100`（使用中で移動できず）と、
計算中の `parallel_pitch16mm_…_d11nm…`・`parallel_pitch32mm_…_d11nm…` はまだ kawatani 直下にある。あとで fieldOLD に移す。

## 電場計算の誤り（fieldOLD と fieldFIX の違い）

`potential.f` の内部の点の係数で、格子の間隔の代わりに点の座標を使っていた
（`drs = sdx1(i,k)`、`sdx1 = 0.5*(x(i+1)+x(i-1))` は間隔ではなく座標）。
そのため解いている式がラプラス方程式になっていなかった。川谷さんの元のコードから引き継がれていた。
fieldFIX では `drs = 0.5*(x(i+1,k)-x(i-1,k))`、`dzs = 0.5*(y(i,k+1)-y(i,k-1))` に直した。
あわせて、電位は液滴の電荷を含まず毎ステップ同じ答えになるので、最初の1回だけ厳しく収束させて解くようにした（計算時間が短くなる）。
