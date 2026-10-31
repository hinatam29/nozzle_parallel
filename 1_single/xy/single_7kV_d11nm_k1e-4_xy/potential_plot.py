# -*- coding: utf-8 -*-
"""
電位分布 φ の図（川谷さんの修論 図 3.2.5 と同じ体裁）
=====================================================

flow.dat の最初のコマから電位 φ を読み、川谷さんの図と同じ段階の値
（0.001, 0.01, 0.1, 1, 10, 100, 1000, 3000, 5000, 7000, 9000 V）で塗り分ける。
ノズルと対向電極は灰色で描く（geom.txt があればその形、なければ川谷さんの形）。

使い方（計算フォルダの中で）
  python3 potential_plot.py                       # r 0〜2 mm, z 80〜100 mm（川谷さんの図と同じ範囲）
  python3 potential_plot.py --rmax 8              # r 0〜8 mm
  python3 potential_plot.py --outdir viz_out/xxx  # 出力先を変える

出力: <outdir>/potential_r<rmax>mm.png
"""
import argparse
import itertools
import os
import re
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Rectangle

# 川谷さんの図 3.2.5 と同じ区切り [V]
LEVELS = [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 3000, 5000, 7000, 9000]
LABELS = ["0.001", "0.01", "0.1", "1", "10", "100", "1000", "3000", "5000",
          "7000", "9000"]
# 同じ色の並び（濃い青 → 青 → 水色 → 緑 → 黄緑 → 黄 → 橙 → 赤 → 赤紫）
COLORS = ["#1f3fd8", "#1a8be8", "#00d8f0", "#00e6a0", "#00cc44", "#a6ee00",
          "#ffd200", "#ff7a00", "#f0002a", "#ff00c8"]
CMAP = ListedColormap(COLORS).with_extremes(under="#0a0aa0", over="#ff00ff")


def to_float(t):
    try:
        return float(t)
    except ValueError:
        return float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", t))


def read_phi(path):
    with open(path, "r", errors="replace") as f:
        f.readline()
        f.readline()
        h = f.readline()
        ni, nj = map(int, re.findall(r"i=\s*(\d+)\s*,j=\s*(\d+)", h)[0])
        a = np.array([[to_float(v) for v in l.split()[:3]]
                      for l in itertools.islice(f, ni * nj)])
    x, y, p = a.T
    if x[1] == x[0]:          # y が先に変わる並び
        X, Y, P = (v.reshape(ni, nj) for v in (x, y, p))
    else:
        X, Y, P = (v.reshape(nj, ni).T for v in (x, y, p))
    r, z = X[:, 0], Y[0, :]
    # flow.dat の座標はセルの中心（格子点から半格子ずれている）。phi は格子点の値
    # なので、格子点の位置に戻す
    r = r - 0.5 * (r[1] - r[0])
    z = z - 0.5 * (z[1] - z[0])
    return r, z, P          # r[i], z[j], phi[i, j]


def read_geom():
    # 既定は川谷さん条件（1本・大きな板の電極）。geom.txt がない古い版では、
    # geo.f の格子番号の決め方により電極の格子点は r = 4.6 mm〜、z = 91.8〜92.6 mm、
    # ノズルは z = 98.8 mm〜。格子点の範囲に半格子（0.1 mm）を足して描く
    g = dict(nopp=1, xlen=0.1, ylen=0.1, xnoz=0.4e-3, xhol=4.5e-3, xring=0.1,
             yr1=91.7e-3, yr2=92.7e-3, yt1=98.7e-3, yt2=0.0, old=True)
    try:
        v = open("geom.txt").read().split()
        g.update(nopp=int(v[0]), xlen=float(v[1]), ylen=float(v[2]),
                 xnoz=float(v[3]), xhol=float(v[4]), xring=float(v[5]),
                 yr1=float(v[6]), yr2=float(v[7]), yt1=float(v[8]),
                 yt2=float(v[9]), old=False)
        # リング版: 格子点 r = xhol〜xring に半格子（0.1 mm）を足して描く
        g["xhol"] -= 0.1e-3
        g["xring"] += 0.1e-3
    except (OSError, ValueError, IndexError):
        pass
    return g


def main():
    ap = argparse.ArgumentParser(description="電位分布の図（川谷さんの図と同じ体裁）")
    ap.add_argument("--flow", default="flow.dat")
    ap.add_argument("--outdir", default="viz_out")
    ap.add_argument("--rmax", type=float, default=2.0, help="r の範囲 [mm]")
    ap.add_argument("--zmin", type=float, default=None)
    ap.add_argument("--zmax", type=float, default=None)
    args = ap.parse_args()

    g = read_geom()
    if args.zmin is None:
        args.zmin = (g["yt2"] * 1e3 - 2) if g["nopp"] == 2 else 80.0
    if args.zmax is None:
        args.zmax = (g["yt1"] * 1e3 + 2) if g["nopp"] == 2 else g["ylen"] * 1e3

    print("flow.dat の最初のコマを読み込み中 ...", file=sys.stderr)
    r, z, P = read_phi(args.flow)
    # 上端は電場 0 の条件（一番上の格子点と同じ値）なので、上端まで同じ値で延ばす
    if g["nopp"] == 1 and z[-1] < g["ylen"]:
        z = np.r_[z, g["ylen"]]
        P = np.c_[P, P[:, -1]]
    r_mm, z_mm = r * 1e3, z * 1e3
    mi = r_mm <= args.rmax + 0.5
    mj = (z_mm >= args.zmin - 0.5) & (z_mm <= args.zmax + 0.5)
    Pp = np.clip(P[np.ix_(mi, mj)].T, 1e-6, None)

    fig = plt.figure(figsize=(4.6, 5.0))
    ax = fig.add_axes([0.17, 0.12, 0.56, 0.80])
    cax = fig.add_axes([0.78, 0.30, 0.04, 0.45])
    norm = BoundaryNorm(LEVELS, CMAP.N)
    cf = ax.contourf(r_mm[mi], z_mm[mj], Pp, levels=LEVELS, cmap=CMAP,
                     norm=norm, extend="both")

    def rect(r0, r1, z0, z1):
        ax.add_patch(Rectangle((r0, z0), r1 - r0, z1 - z0, facecolor="0.82",
                               edgecolor="0.45", lw=0.6, zorder=4))
    m = 1e3
    rect(0, g["xnoz"] * m, g["yt1"] * m, g["ylen"] * m)
    rect(g["xhol"] * m, min(g["xring"], g["xlen"]) * m, g["yr1"] * m,
         g["yr2"] * m)
    if g["nopp"] == 2:
        rect(0, g["xnoz"] * m, 0, g["yt2"] * m)

    ax.set_xlim(0, args.rmax)
    ax.set_ylim(args.zmin, args.zmax)
    ax.set_xlabel("x [mm]", style="italic")
    ax.set_ylabel("y [mm]", style="italic")
    cb = fig.colorbar(cf, cax=cax, ticks=LEVELS)
    cb.ax.set_yticklabels(LABELS, fontsize=7)
    cb.ax.set_title(r"$\Phi$ [V]", fontsize=8, loc="left")
    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, f"potential_r{args.rmax:g}mm.png")
    fig.savefig(out, dpi=200)
    print("  保存:", out, file=sys.stderr)


if __name__ == "__main__":
    main()
