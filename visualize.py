# -*- coding: utf-8 -*-
"""
静電噴霧シミュレーション結果の可視化（先行研究フォーマット準拠）
==================================================================

平面2次元・並列2ノズルの噴霧粒子を、先行研究と同じ体裁で2通り可視化します。

  (1) 液滴の広がり  : 液滴を「粒径 d [m]」で色分けした散布図（先行研究 図1）
  (2) 二酸化炭素吸収 : 液滴の吸収量 clt を格子に集計した「C_total [mol]」の
                       塗り等高線（先行研究 図2の右半分）

いずれもノズル・対向電極を灰色の四角で描き、軸は r [mm] / z [mm]。
左右対称なので既定では片側ノズル（r=0側）に拡大して先行研究と比較しやすくします。

読み込み: spray.dat（xp, yp, dp, clt, pout）

使い方（例）
  python visualize.py
  python visualize.py --rmax 8            # 表示する r 範囲[mm]（片ノズル拡大）
  python visualize.py --full              # 全幅（両ノズル）表示
  python visualize.py --kind spread       # 広がりだけ
  python visualize.py --kind absorb       # 吸収だけ
  python visualize.py --max-frames 18     # 正常範囲だけ（発散フレーム除外）

必要ライブラリ : numpy, matplotlib
"""

import argparse
import os
import re
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.patches import Rectangle
from matplotlib.colors import BoundaryNorm, LogNorm


# ---- 幾何形状（geo.f と対応、単位 m）----------------------------------
YLEN = 0.1
XLEN = 0.1
XNOZ = 0.4e-3       # ノズル半径方向幅
YNOZ = 1.0e-3       # ノズル軸方向長さ（上端から）
XHOL = 5.0e-3       # 対向電極の穴（r<XHOL, および r>XLEN-XHOL）
YHOL1 = 93.0e-3
YHOL2 = 92.0e-3
DX = XLEN / 500.0   # 格子間隔（nx=500）


def setup_font():
    from matplotlib import font_manager
    for name in ["Yu Gothic", "Meiryo", "MS Gothic", "MS PGothic",
                 "Hiragino Sans", "Noto Sans CJK JP", "Noto Sans JP",
                 "IPAexGothic", "IPAGothic", "TakaoGothic", "VL Gothic"]:
        if name in {f.name for f in font_manager.fontManager.ttflist}:
            matplotlib.rcParams["font.family"] = name
            matplotlib.rcParams["axes.unicode_minus"] = False
            return True
    return False


JP = setup_font()
L = {
    "r": "r [mm]", "z": "z [mm]",
    "d": "d [m]", "ct": "C_total [mol]",
    "step": "ステップ" if JP else "step",
    "t_spread": "液滴の広がり" if JP else "Droplet spread",
    "t_absorb": "二酸化炭素吸収" if JP else "CO2 absorption",
}


# ---- Tecplot パーサ（3桁指数E省略・NaN対応）--------------------------
_exp_fix = re.compile(r"(\d\.\d+)([+-]\d{2,3})(?![\d.])")


def to_float(t):
    try:
        return float(t.replace("D", "E").replace("d", "e"))
    except ValueError:
        return np.nan


def read_spray(path, stride=1, max_frames=None):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    varnames, frames, zidx = None, [], -1
    with open(path, "r", errors="replace") as fh:
        pend, buf, keep = None, [], False

        def flush():
            nonlocal pend, buf
            if pend is None:
                return
            need = pend["n"] * (pend["nc"] or 0)
            if keep and pend["nc"] and len(buf) >= need:
                arr = np.array(buf[:need]).reshape(pend["n"], pend["nc"])
                frames.append(dict(istp=pend["istp"], data=arr))
            pend, buf = None, []

        for line in fh:
            s = line.strip()
            if not s:
                continue
            low = s.lower()
            if low.startswith("title"):
                continue
            if low.startswith("variable"):
                varnames = [n.strip() for n in re.findall(r'"([^"]*)"', s)]
                continue
            if low.startswith("zone"):
                flush()
                zidx += 1
                m = re.search(r"n=\s*(\d+)", s)
                mi = re.search(r"\bi\s*=\s*(\d+)", s, re.I)
                keep = (zidx % stride == 0) and (
                    max_frames is None or len(frames) < max_frames)
                pend = dict(istp=int(m.group(1)) if m else zidx,
                            n=int(mi.group(1)) if mi else 1,
                            nc=len(varnames) if varnames else None)
                buf = []
                continue
            if pend is None or not keep:
                continue
            toks = _exp_fix.sub(r"\1E\2", s).split()
            if pend["nc"] is None:
                pend["nc"] = len(toks)
            buf.extend(to_float(t) for t in toks)
        flush()
    if varnames is None:
        raise ValueError(f"{path}: variables 行が見つかりません")
    print(f"  -> {len(frames)} フレーム読み込み ({path})", file=sys.stderr)
    return varnames, frames


def active(varnames, frame):
    idx = {n: i for i, n in enumerate(varnames)}
    d = frame["data"]
    xp = d[:, idx.get("xp", 0)]
    yp = d[:, idx.get("yp", 1)]
    dp = d[:, idx["dp"]] if "dp" in idx else np.zeros(len(xp))
    clt = d[:, idx["clt"]] if "clt" in idx else np.zeros(len(xp))
    pout = d[:, idx["pout"]] if "pout" in idx else np.zeros(len(xp))
    m = np.isfinite(xp) & np.isfinite(yp) & (pout == 0)
    return xp[m], yp[m], dp[m], clt[m]


# ---- 灰色の四角（ノズル・電極）----------------------------------------
def draw_geometry(ax):
    """全幾何を絶対座標[mm]で描く。各パネルは xlim で切り取る。
       ノズル1(r=0), ノズル2(r=XLEN), 対向電極(r=XHOL..XLEN-XHOL)。"""
    def rect(r0, r1, z0, z1):
        ax.add_patch(Rectangle((r0, z0), r1 - r0, z1 - z0,
                               facecolor="0.8", edgecolor="0.5",
                               linewidth=0.6, zorder=4))
    znoz0, znoz1 = (YLEN - YNOZ) * 1e3, YLEN * 1e3          # 99..100
    zel0, zel1 = YHOL2 * 1e3, YHOL1 * 1e3                    # 92..93
    rect(0.0, XNOZ * 1e3, znoz0, znoz1)                      # ノズル1
    rect((XLEN - XNOZ) * 1e3, XLEN * 1e3, znoz0, znoz1)      # ノズル2
    rect(XHOL * 1e3, (XLEN - XHOL) * 1e3, zel0, zel1)        # 電極(両端に穴)


# ---- 広がり（粒径 d で色分け）-----------------------------------------
D_LEVELS = np.arange(2.0, 8.5, 0.5) * 1e-7   # 2E-7 .. 8E-7 (先行研究準拠)


def plot_spread(ax, varnames, frame):
    xp, yp, dp, clt = active(varnames, frame)
    cmap = plt.cm.jet
    norm = BoundaryNorm(D_LEVELS, cmap.N, extend="both")
    sc = None
    if len(xp):
        sc = ax.scatter(xp * 1e3, yp * 1e3, s=6, c=dp, cmap=cmap, norm=norm,
                        edgecolors="none", zorder=3)
    draw_geometry(ax)
    return sc


# ---- 吸収（C_total 塗り等高線）----------------------------------------
def bin_ctotal(xp, yp, clt, rmax_m, zlim_m, bin_m):
    nx = max(int(rmax_m / bin_m), 4)
    nz = max(int((zlim_m[1] - zlim_m[0]) / bin_m), 4)
    H, xe, ye = np.histogram2d(
        xp, yp, bins=[nx, nz],
        range=[[0, rmax_m], [zlim_m[0], zlim_m[1]]], weights=clt)
    # 軽い平滑化（3x3 平均、numpy のみ）
    Hs = H.copy()
    Hs[1:-1, 1:-1] = (H[:-2, 1:-1] + H[2:, 1:-1] + H[1:-1, :-2] + H[1:-1, 2:]
                      + H[1:-1, 1:-1] + H[:-2, :-2] + H[:-2, 2:]
                      + H[2:, :-2] + H[2:, 2:]) / 9.0
    xc = 0.5 * (xe[:-1] + xe[1:])
    yc = 0.5 * (ye[:-1] + ye[1:])
    return xc, yc, Hs


def plot_absorb(ax, varnames, frame, zlim, cmax):
    xp, yp, dp, clt = active(varnames, frame)
    zlim_m = (zlim[0] * 1e-3, zlim[1] * 1e-3)
    cf = None
    if len(xp) and cmax > 0:
        xc, yc, H = bin_ctotal(xp, yp, clt, XLEN, zlim_m, 2 * DX)
        levels = np.logspace(np.log10(cmax) - 3, np.log10(cmax), 12)
        Hp = np.clip(H.T, levels[0] * 0.5, None)
        cf = ax.contourf(xc * 1e3, yc * 1e3, Hp, levels=levels,
                         cmap="jet", norm=LogNorm(levels[0], levels[-1]),
                         extend="both", zorder=2)
    draw_geometry(ax)
    return cf


# ---- スケール ---------------------------------------------------------
def clt_cmax(varnames, frames, zlim, bin_m):
    zlim_m = (zlim[0] * 1e-3, zlim[1] * 1e-3)
    mx = 0.0
    for fr in frames:
        xp, yp, dp, clt = active(varnames, fr)
        if len(xp):
            _, _, H = bin_ctotal(xp, yp, clt, XLEN, zlim_m, bin_m)
            if H.size and H.max() > 0:
                mx = max(mx, float(H.max()))
    return mx


def build_fig(single):
    if single:
        fig = plt.figure(figsize=(4.8, 6.6))
        axes = [fig.add_axes([0.16, 0.10, 0.60, 0.80])]
        cax = fig.add_axes([0.80, 0.10, 0.03, 0.80])
    else:
        fig = plt.figure(figsize=(7.6, 6.6))
        axes = [fig.add_axes([0.10, 0.10, 0.34, 0.78]),
                fig.add_axes([0.47, 0.10, 0.34, 0.78])]
        cax = fig.add_axes([0.86, 0.10, 0.025, 0.78])
    return fig, axes, cax


def render_into(fig, axes, cax, kind, vns, frame, zlim, rmax_mm, cmax):
    for ax in axes:
        ax.clear()
    cax.cla()
    if len(axes) == 1:
        xlims = [(0, rmax_mm)]
    else:
        xlims = [(0, rmax_mm), (XLEN * 1e3 - rmax_mm, XLEN * 1e3)]
    mapp = None
    for ax, xl in zip(axes, xlims):
        m = (plot_spread(ax, vns, frame) if kind == "spread"
             else plot_absorb(ax, vns, frame, zlim, cmax))
        if m is not None:
            mapp = m
        ax.set_xlim(*xl)
        ax.set_ylim(*zlim)
        ax.set_xlabel(L["r"])
    axes[0].set_ylabel(L["z"])
    if len(axes) == 2:
        axes[1].set_yticklabels([])
        axes[0].set_title("ノズル1" if JP else "nozzle 1", fontsize=10)
        axes[1].set_title("ノズル2" if JP else "nozzle 2", fontsize=10)
    title = L["t_spread"] if kind == "spread" else L["t_absorb"]
    fig.suptitle(f'{title}   ({L["step"]}={frame.get("istp")})', y=0.985)
    if mapp is not None:
        plt.colorbar(mapp, cax=cax,
                     label=(L["d"] if kind == "spread" else L["ct"]))
    else:
        cax.axis("off")


def main():
    ap = argparse.ArgumentParser(description="先行研究フォーマットの噴霧可視化")
    ap.add_argument("--spray", default="spray.dat")
    ap.add_argument("--outdir", default="viz_out")
    ap.add_argument("--kind", default="both",
                    choices=["both", "spread", "absorb"])
    ap.add_argument("--mode", default="both", choices=["both", "anim", "snap"])
    ap.add_argument("--rmax", type=float, default=8.0,
                    help="各ノズルパネルの r 範囲 [mm]")
    ap.add_argument("--single", action="store_true",
                    help="ノズル1側だけ1枚で表示（既定は2本とも）")
    ap.add_argument("--zmin", type=float, default=86.0)
    ap.add_argument("--zmax", type=float, default=100.0)
    ap.add_argument("--nsnap", type=int, default=6)
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--max-frames", type=int, default=None)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print("噴霧データ読み込み中 ...", file=sys.stderr)
    vns, frames = read_spray(args.spray, stride=args.stride,
                             max_frames=args.max_frames)
    if not frames:
        print("フレームがありません。", file=sys.stderr)
        return

    zlim = (args.zmin, args.zmax)
    single = args.single
    cmax = clt_cmax(vns, frames, zlim, 2 * DX)
    print(f"  r範囲/パネル=0..{args.rmax}mm, z={zlim}mm, "
          f"C_total上限={cmax:.2e}, {'1枚' if single else '2ノズル'}",
          file=sys.stderr)

    kinds = ["spread", "absorb"] if args.kind == "both" else [args.kind]
    for kd in kinds:
        print(f"[{kd}] 作成中 ...", file=sys.stderr)
        if args.mode in ("snap", "both"):
            fig, axes, cax = build_fig(single)
            picks = np.unique(np.linspace(0, len(frames) - 1,
                              min(args.nsnap, len(frames))).astype(int))
            for k in picks:
                render_into(fig, axes, cax, kd, vns, frames[k], zlim,
                            args.rmax, cmax)
                p = os.path.join(args.outdir,
                                 f"{kd}_snap_{k:03d}_istp{frames[k]['istp']}.png")
                fig.savefig(p, dpi=150)
                print(f"  保存: {p}", file=sys.stderr)
            plt.close(fig)
        if args.mode in ("anim", "both"):
            fig, axes, cax = build_fig(single)

            def upd(k, kd=kd, fig=fig, axes=axes, cax=cax):
                render_into(fig, axes, cax, kd, vns, frames[k], zlim,
                            args.rmax, cmax)
                return []

            anim = animation.FuncAnimation(fig, upd, frames=len(frames),
                                           blit=False)
            mp4 = os.path.join(args.outdir, f"{kd}_animation.mp4")
            try:
                anim.save(mp4, writer=animation.FFMpegWriter(fps=args.fps,
                          bitrate=2400), dpi=140)
                print(f"  保存: {mp4}", file=sys.stderr)
            except Exception as e:
                gif = os.path.join(args.outdir, f"{kd}_animation.gif")
                print(f"  MP4不可({e})->GIF", file=sys.stderr)
                anim.save(gif, writer=animation.PillowWriter(fps=args.fps),
                          dpi=100)
                print(f"  保存: {gif}", file=sys.stderr)
            plt.close(fig)
    print("完了。出力先: " + os.path.abspath(args.outdir), file=sys.stderr)


if __name__ == "__main__":
    main()
