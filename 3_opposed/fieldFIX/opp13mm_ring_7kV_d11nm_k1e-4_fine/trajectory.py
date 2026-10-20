# -*- coding: utf-8 -*-
"""
液滴の軌跡の重ね合わせ図（高速度カメラの多数フレーム重ね合わせと同じ見方）
==========================================================================

spray_fine.dat（10 µs ごとの液滴の位置）を読み、同じ液滴の位置を線でつないで
軌跡を描く。黒い背景に白い線で、実験の重ね合わせ画像と同じ体裁にする。
円筒（軸対称）なので、r を左右に折り返して（-r と +r）噴霧全体の形で描く。

使い方（計算フォルダの中で）
  python3 trajectory.py                   # 最初の 3 ms と最後の 1 ms の2枚
  python3 trajectory.py --tmin 0 --tmax 1 # 好きな時間の範囲 [ms]
  python3 trajectory.py --color           # 線を液滴径で色分け
  python3 trajectory.py --gif             # 最初の 3 ms の動き（10 µs ごと）を GIF に

出力: viz_traj/ フォルダ
必要ライブラリ: numpy, matplotlib, pillow（GIF のとき）
"""
import argparse
import os
import re
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from matplotlib.colors import BoundaryNorm

DT = 2.0e-7          # 時間刻み [s]（main.f と同じ）
D_LEVELS = np.arange(4.0, 30.0, 2.0)   # 液滴径の色 [nm]（visualize.py と同じ）


def to_float(t):
    # Fortran の3桁指数（E の省略、例 0.123-105）にも対応
    try:
        return float(t)
    except ValueError:
        return float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", t))


def read_frames(path, tmin, tmax):
    """[tmin, tmax] ms の範囲のフレームだけ読む。(時刻 ms, 配列) のリスト"""
    frames = []
    cur_t, rows, keep = None, [], False
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if "zone" in line:
                if keep and rows:
                    frames.append((cur_t, np.array(rows)))
                m = re.search(r"n=\s*(\d+)", line)
                istp = int(m.group(1)) if m else 0
                cur_t = (istp - 1) * DT * 1e3
                keep = tmin - 1e-9 <= cur_t <= tmax + 1e-9
                rows = []
            elif keep:
                p = line.split()
                if len(p) == 5:
                    rows.append([to_float(v) for v in p])
    if keep and rows:
        frames.append((cur_t, np.array(rows)))
    return frames


def read_geom():
    g = dict(nopp=1, xlen=0.1, ylen=0.1, xnoz=0.4e-3, xhol=5e-3, xring=6e-3,
             yr1=92.1e-3, yr2=93.1e-3, yt1=99e-3, yt2=0.0)
    try:
        v = open("geom.txt").read().split()
        g.update(nopp=int(v[0]), xlen=float(v[1]), ylen=float(v[2]),
                 xnoz=float(v[3]), xhol=float(v[4]), xring=float(v[5]),
                 yr1=float(v[6]), yr2=float(v[7]), yt1=float(v[8]),
                 yt2=float(v[9]))
    except (OSError, ValueError, IndexError):
        print("  geom.txt がないので既定の形で描きます", file=sys.stderr)
    return g


def draw_geom(ax, g, rmax):
    def rect(r0, r1, z0, z1):
        for s in (1, -1):   # 左右に折り返す
            a, b = sorted((s * r0, s * r1))
            ax.add_patch(Rectangle((a, z0), b - a, z1 - z0, facecolor="0.55",
                                   edgecolor="0.75", lw=0.6, zorder=5))
    m = 1e3
    rect(0, g["xnoz"] * m, g["yt1"] * m, g["ylen"] * m)          # 上のノズル
    rect(g["xhol"] * m, min(g["xring"] * m, rmax), g["yr1"] * m, g["yr2"] * m)
    if g["nopp"] == 2:
        rect(0, g["xnoz"] * m, 0, g["yt2"] * m)                  # 下のノズル


def tracks(frames):
    """液滴番号ごとに (r, z, d, 飛行中か) の列をつくる"""
    n = max(len(a) for _, a in frames)
    T = len(frames)
    R = np.full((T, n), np.nan)
    Z = np.full((T, n), np.nan)
    Dd = np.full(n, np.nan)
    for it, (_, a) in enumerate(frames):
        k = len(a)
        fly = a[:, 4] == 0
        R[it, :k] = np.where(fly, a[:, 0], np.nan)
        Z[it, :k] = np.where(fly, a[:, 1], np.nan)
        Dd[:k] = a[:, 2]
    return R * 1e3, Z * 1e3, Dd * 1e9


def plot(frames, g, out, title, color, zlim, rmax):
    R, Z, Dd = tracks(frames)
    segs, cols = [], []
    for s in (1, -1):
        r0, r1 = s * R[:-1], s * R[1:]
        z0, z1 = Z[:-1], Z[1:]
        ok = np.isfinite(r0) & np.isfinite(r1)
        it, ik = np.nonzero(ok)
        segs.append(np.stack([np.stack([r0[it, ik], z0[it, ik]], -1),
                              np.stack([r1[it, ik], z1[it, ik]], -1)], 1))
        cols.append(Dd[ik])
    segs = np.concatenate(segs)
    cols = np.concatenate(cols)
    fig, ax = plt.subplots(figsize=(6.0, 6.8))
    ax.set_facecolor("black")
    if color:
        norm = BoundaryNorm(D_LEVELS, plt.cm.jet.N, extend="both")
        lc = LineCollection(segs, cmap=plt.cm.jet, norm=norm, lw=0.5,
                            alpha=0.8)
        lc.set_array(cols)
        ax.add_collection(lc)
        fig.colorbar(lc, ax=ax, label="d [nm]", fraction=0.05)
    else:
        ax.add_collection(LineCollection(segs, colors="white", lw=0.35,
                                         alpha=0.5))
    draw_geom(ax, g, rmax)
    ax.set_xlim(-rmax, rmax)
    ax.set_ylim(*zlim)
    ax.set_aspect("equal")
    ax.set_xlabel("r [mm]")
    ax.set_ylabel("z [mm]")
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=170, facecolor="white")
    plt.close(fig)
    print("  保存:", out, " 線の数:", len(segs), file=sys.stderr)


def gif(frames, g, out, zlim, rmax, fps):
    from PIL import Image
    imgs = []
    for t, a in frames:
        fig, ax = plt.subplots(figsize=(5.0, 5.6))
        ax.set_facecolor("black")
        fly = a[:, 4] == 0
        r, z = a[fly, 0] * 1e3, a[fly, 1] * 1e3
        ax.scatter(np.r_[r, -r], np.r_[z, z], s=2, c="white", lw=0)
        draw_geom(ax, g, rmax)
        ax.set_xlim(-rmax, rmax)
        ax.set_ylim(*zlim)
        ax.set_aspect("equal")
        ax.set_xlabel("r [mm]")
        ax.set_ylabel("z [mm]")
        ax.set_title("t = %.2f ms" % t, fontsize=10)
        fig.tight_layout()
        fig.canvas.draw()
        imgs.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3]))
        plt.close(fig)
    pal = imgs[-1].convert("P", palette=Image.ADAPTIVE, colors=64)
    q = [im.quantize(palette=pal) for im in imgs]
    q[0].save(out, save_all=True, append_images=q[1:], duration=int(1000 / fps),
              loop=0)
    print("  保存:", out, " コマ数:", len(q), file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description="液滴の軌跡の重ね合わせ図")
    ap.add_argument("--spray", default="spray_fine.dat")
    ap.add_argument("--outdir", default="viz_traj")
    ap.add_argument("--tmin", type=float, default=None, help="開始 [ms]")
    ap.add_argument("--tmax", type=float, default=None, help="終了 [ms]")
    ap.add_argument("--rmax", type=float, default=8.0, help="表示する半径 [mm]")
    ap.add_argument("--zmin", type=float, default=None)
    ap.add_argument("--zmax", type=float, default=None)
    ap.add_argument("--color", action="store_true", help="液滴径で色分け")
    ap.add_argument("--gif", action="store_true", help="最初の 3 ms の動きを GIF に")
    ap.add_argument("--fps", type=int, default=20)
    args = ap.parse_args()

    g = read_geom()
    if args.zmin is None:
        args.zmin = (g["yt2"] * 1e3 - 2) if g["nopp"] == 2 else g["yr1"] * 1e3 - 6
    if args.zmax is None:
        args.zmax = (g["yt1"] * 1e3 + 2) if g["nopp"] == 2 else g["ylen"] * 1e3
    zlim = (args.zmin, args.zmax)
    os.makedirs(args.outdir, exist_ok=True)

    if args.tmin is not None or args.tmax is not None:
        wins = [(args.tmin or 0.0, args.tmax if args.tmax is not None else 1e9)]
    else:
        wins = [(0.0, 3.0), (10.0, 11.0)]   # 最初の 3 ms と最後の 1 ms
    for t0, t1 in wins:
        print(f"読み込み中 {t0}〜{t1} ms ...", file=sys.stderr)
        fr = read_frames(args.spray, t0, t1)
        if len(fr) < 2:
            print("  この範囲のフレームが 2 枚未満なので飛ばします", file=sys.stderr)
            continue
        ta, tb = fr[0][0], fr[-1][0]
        tag = "c" if args.color else "w"
        out = os.path.join(args.outdir, f"traj_{ta:.2f}-{tb:.2f}ms_{tag}.png")
        plot(fr, g, out,
             f"droplet trajectories  {ta:.2f}-{tb:.2f} ms ({len(fr)} frames)",
             args.color, zlim, args.rmax)
        if args.gif and t0 == wins[0][0]:
            gif(fr, g, os.path.join(args.outdir, f"motion_{ta:.2f}-{tb:.2f}ms.gif"),
                zlim, args.rmax, args.fps)
    print("完了。出力先:", os.path.abspath(args.outdir), file=sys.stderr)


if __name__ == "__main__":
    main()
