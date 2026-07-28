# -*- coding: utf-8 -*-
"""
静電噴霧シミュレーション結果の可視化（先行研究フォーマット準拠）
==================================================================

平面2次元・並列2ノズルの噴霧粒子を、先行研究と同じ体裁で2通り可視化します。

  (1) 液滴の広がり  : 液滴を「粒径 d [m]」で色分けした散布図（先行研究 図1）
  (2) 二酸化炭素吸収 : 液滴の吸収量 clt を格子集計した「C_total [mol]」の
                       塗り等高線（先行研究 図2の右半分）

既定は「1枚の絵に両ノズル」。ノズル・対向電極は灰色の四角、軸は r/z [mm]。
ノズル間隔（領域幅 xlen）は spray.dat から自動判定します。

読み込み: spray.dat（xp, yp, dp, clt, pout）

使い方（例）
  python visualize.py
  python visualize.py --dual          # 左右2パネルで各ノズルを拡大
  python visualize.py --kind spread   # 広がりだけ
  python visualize.py --max-frames 18 # 正常範囲だけ（発散フレーム除外）

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
from matplotlib.cm import ScalarMappable


# ---- 幾何形状（geo.f と対応、単位 m）。xlen はデータから自動判定 ------
YLEN = 0.1
XNOZ = 0.4e-3       # ノズル半径方向幅
YNOZ = 1.0e-3       # ノズル軸方向長さ（上端から）
XHOL = 5.0e-3       # 対向電極の穴（各ノズル側 XHOL 幅）
YHOL1 = 93.0e-3
YHOL2 = 92.0e-3


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
    "r": "r [mm]", "z": "z [mm]", "d": "d [m]", "ct": "C_total [mol]",
    "step": "ステップ" if JP else "step",
    "t_spread": "液滴の広がり" if JP else "Droplet spread",
    "t_absorb": "二酸化炭素吸収" if JP else "CO2 absorption",
}
D_LEVELS = np.arange(2.0, 8.5, 0.5) * 1e-7   # 2E-7 .. 8E-7（先行研究準拠）


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


def infer_xlen(varnames, frames):
    mx = 0.0
    for fr in frames:
        xp, _, _, _ = active(varnames, fr)
        if len(xp):
            mx = max(mx, float(xp.max()))
    return mx if mx > 0 else 0.1


# ---- 灰色の四角（ノズル1・ノズル2・対向電極）---------------------------
def draw_geometry(ax, xlen_m):
    xl = xlen_m * 1e3

    def rect(r0, r1, z0, z1):
        ax.add_patch(Rectangle((r0, z0), r1 - r0, z1 - z0, facecolor="0.8",
                               edgecolor="0.5", linewidth=0.6, zorder=4))
    zn0, zn1 = (YLEN - YNOZ) * 1e3, YLEN * 1e3
    ze0, ze1 = YHOL2 * 1e3, YHOL1 * 1e3
    rect(0.0, XNOZ * 1e3, zn0, zn1)                    # ノズル1
    rect(xl - XNOZ * 1e3, xl, zn0, zn1)                # ノズル2
    rect(XHOL * 1e3, xl - XHOL * 1e3, ze0, ze1)        # 電極（両端に穴）


# ---- C_total の格子集計 ----------------------------------------------
def bin_ctotal(xp, yp, clt, xlen_m, zlim_m, bin_m):
    nx = max(int(xlen_m / bin_m), 4)
    nz = max(int((zlim_m[1] - zlim_m[0]) / bin_m), 4)
    H, xe, ye = np.histogram2d(xp, yp, bins=[nx, nz],
                               range=[[0, xlen_m], list(zlim_m)], weights=clt)
    Hs = H.copy()
    Hs[1:-1, 1:-1] = (H[:-2, 1:-1] + H[2:, 1:-1] + H[1:-1, :-2] + H[1:-1, 2:]
                      + H[1:-1, 1:-1] + H[:-2, :-2] + H[:-2, 2:]
                      + H[2:, :-2] + H[2:, 2:]) / 9.0
    xc = 0.5 * (xe[:-1] + xe[1:])
    yc = 0.5 * (ye[:-1] + ye[1:])
    return xc, yc, Hs


# ---- 1パネルの中身を描く（カラーバーは触らない）-----------------------
def plot_spread(ax, vns, frame, xlen_m, norm):
    xp, yp, dp, clt = active(vns, frame)
    if len(xp):
        ax.scatter(xp * 1e3, yp * 1e3, s=6, c=dp, cmap=plt.cm.jet, norm=norm,
                   edgecolors="none", zorder=3)
    draw_geometry(ax, xlen_m)


def plot_absorb(ax, vns, frame, xlen_m, zlim, bin_m, levels, norm):
    xp, yp, dp, clt = active(vns, frame)
    zlim_m = (zlim[0] * 1e-3, zlim[1] * 1e-3)
    if len(xp):
        xc, yc, H = bin_ctotal(xp, yp, clt, xlen_m, zlim_m, bin_m)
        Hp = np.clip(H.T, levels[0] * 0.5, None)
        ax.contourf(xc * 1e3, yc * 1e3, Hp, levels=levels, cmap="jet",
                    norm=norm, extend="both", zorder=2)
    draw_geometry(ax, xlen_m)


def clt_cmax(vns, frames, xlen_m, zlim, bin_m):
    zlim_m = (zlim[0] * 1e-3, zlim[1] * 1e-3)
    mx = 0.0
    for fr in frames:
        xp, yp, dp, clt = active(vns, fr)
        if len(xp):
            _, _, H = bin_ctotal(xp, yp, clt, xlen_m, zlim_m, bin_m)
            if H.size and H.max() > 0:
                mx = max(mx, float(H.max()))
    return mx


# ---- 図の構築（カラーバーは固定スケールで一度だけ作る）----------------
def build(kind, dual, xlen_m, norm, cmap, label):
    if dual:
        fig = plt.figure(figsize=(7.6, 6.6))
        axes = [fig.add_axes([0.10, 0.10, 0.34, 0.78]),
                fig.add_axes([0.47, 0.10, 0.34, 0.78])]
        cax = fig.add_axes([0.86, 0.10, 0.025, 0.78])
    else:
        w = 9.5 if xlen_m > 0.06 else 7.0
        fig = plt.figure(figsize=(w, 4.6))
        axes = [fig.add_axes([0.08, 0.14, 0.80, 0.74])]
        cax = fig.add_axes([0.90, 0.14, 0.02, 0.74])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    fig.colorbar(sm, cax=cax, label=label)   # 固定：以後は再作成しない
    return fig, axes


def render(fig, axes, kind, vns, frame, xlen_m, zlim, rmax_mm,
           norm, bin_m, levels):
    for ax in axes:
        ax.clear()
    if len(axes) == 1:
        xlims = [(0, xlen_m * 1e3)]
    else:
        xlims = [(0, rmax_mm), (xlen_m * 1e3 - rmax_mm, xlen_m * 1e3)]
    for ax, xl in zip(axes, xlims):
        if kind == "spread":
            plot_spread(ax, vns, frame, xlen_m, norm)
        else:
            plot_absorb(ax, vns, frame, xlen_m, zlim, bin_m, levels, norm)
        ax.set_xlim(*xl)
        ax.set_ylim(*zlim)
        ax.set_xlabel(L["r"])
    axes[0].set_ylabel(L["z"])
    if len(axes) == 2:
        axes[1].set_yticklabels([])
        axes[0].set_title("ノズル1" if JP else "nozzle 1", fontsize=10)
        axes[1].set_title("ノズル2" if JP else "nozzle 2", fontsize=10)
    title = L["t_spread"] if kind == "spread" else L["t_absorb"]
    fig.suptitle(f'{title}   ({L["step"]}={frame.get("istp")})', y=0.98)


def read_flow_phi(path):
    """flow.dat の先頭1ゾーンから x, y, phi を読み、2D配列(nj,ni)で返す。
       電位はソース項が無くラプラス方程式（幾何とBCのみで決まる）なので
       時間で変わらず、1ゾーンで電場を代表できる。"""
    if not os.path.exists(path):
        return None
    varnames, ni, nj, buf = None, None, None, []
    with open(path, "r", errors="replace") as fh:
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
                if ni is not None:
                    break               # 2つ目のzoneで打ち切り
                mi = re.search(r"\bi\s*=\s*(\d+)", s, re.I)
                mj = re.search(r"\bj\s*=\s*(\d+)", s, re.I)
                ni = int(mi.group(1)); nj = int(mj.group(1))
                continue
            if ni is None:
                continue
            buf.extend(to_float(t) for t in _exp_fix.sub(r"\1E\2", s).split())
            if len(buf) >= ni * nj * len(varnames):
                break
    if not varnames or ni is None or len(buf) < ni * nj * len(varnames):
        return None
    nc = len(varnames)
    a = np.array(buf[:ni * nj * nc]).reshape(nj, ni, nc)
    idx = {n: i for i, n in enumerate(varnames)}
    X = a[:, :, idx.get("x", 0)]
    Y = a[:, :, idx.get("y", 1)]
    PHI = a[:, :, idx["phi"]] if "phi" in idx else a[:, :, 2]
    return X, Y, PHI


def plot_efield(ax, X, Y, PHI, xlen_m):
    """電位 phi から電場 E=-∇φ を計算し、|E| の等高線＋電気力線を描く。"""
    xc = X[0, :]
    zc = Y[:, 0]
    dpz, dpx = np.gradient(PHI, zc, xc)
    Ex, Ez = -dpx, -dpz
    Emag = np.sqrt(Ex * Ex + Ez * Ez)
    pos = Emag[np.isfinite(Emag) & (Emag > 0)]
    vmax = float(np.percentile(pos, 99)) if pos.size else 1.0
    vmin = vmax / 1e3
    levels = np.logspace(np.log10(vmin), np.log10(vmax), 12)
    Ep = np.clip(np.nan_to_num(Emag, nan=vmin), vmin, vmax)
    cf = ax.contourf(xc * 1e3, zc * 1e3, Ep, levels=levels, cmap="viridis",
                     norm=LogNorm(vmin, vmax), extend="both", zorder=1)
    try:
        ax.streamplot(xc * 1e3, zc * 1e3, Ex, Ez, color="white",
                      density=1.1, linewidth=0.5, arrowsize=0.7)
    except Exception:
        pass
    draw_geometry(ax, xlen_m)
    return cf


EPS0 = 8.8542e-12
FS_SURF = 42.9e-3          # 表面張力（intraction.f と同じ）


def droplet_charge(dp):
    """レイリー限界の電荷 q = pi*sqrt(8*eps0*fs*dp^3)（コードと同じ）。"""
    return np.pi * np.sqrt(8.0 * EPS0 * FS_SURF * np.clip(dp, 0, None) ** 3)


def efield_from_charges(xp, yp, q, Xg, Zg, xlen_m, soft):
    """帯電液滴群が格子点に作る電場 E=(1/4πε0)Σ q(r-r_k)/|r-r_k|^3。
       対称面 x=0, x=xlen に対する鏡像も含める（動力学と整合）。"""
    kc = 1.0 / (4.0 * np.pi * EPS0)
    sx = np.concatenate([xp, -xp, 2.0 * xlen_m - xp])   # 実体＋鏡像2枚
    sy = np.concatenate([yp, yp, yp])
    sq = np.concatenate([q, q, q])
    gx = Xg.ravel()
    gz = Zg.ravel()
    ex = np.zeros(gx.size)
    ez = np.zeros(gx.size)
    s2 = soft * soft
    ch = 200
    for i in range(0, gx.size, ch):
        dx = gx[i:i + ch][:, None] - sx[None, :]
        dz = gz[i:i + ch][:, None] - sy[None, :]
        r2 = dx * dx + dz * dz + s2
        inv = sq[None, :] / (r2 * np.sqrt(r2))
        ex[i:i + ch] = kc * np.sum(inv * dx, axis=1)
        ez[i:i + ch] = kc * np.sum(inv * dz, axis=1)
    return ex.reshape(Xg.shape), ez.reshape(Xg.shape)


def main():
    ap = argparse.ArgumentParser(description="先行研究フォーマットの噴霧可視化")
    ap.add_argument("--spray", default="spray.dat")
    ap.add_argument("--flow", default="flow.dat")
    ap.add_argument("--outdir", default="viz_out")
    ap.add_argument("--kind", default="all",
                    choices=["all", "both", "spread", "absorb", "efield",
                             "efieldt"])
    ap.add_argument("--mode", default="both", choices=["both", "anim", "snap"])
    ap.add_argument("--dual", action="store_true",
                    help="左右2パネルで各ノズルを拡大表示")
    ap.add_argument("--rmax", type=float, default=8.0,
                    help="--dual 時の各パネル r 範囲 [mm]")
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

    xlen_m = infer_xlen(vns, frames)
    zlim = (args.zmin, args.zmax)
    bin_m = xlen_m / 200.0
    cmax = clt_cmax(vns, frames, xlen_m, zlim, bin_m)
    levels = np.logspace(np.log10(cmax) - 3, np.log10(cmax), 12) \
        if cmax > 0 else np.logspace(-24, -21, 12)
    norm_s = BoundaryNorm(D_LEVELS, plt.cm.jet.N, extend="both")
    norm_a = LogNorm(levels[0], levels[-1])
    print(f"  ノズル間隔(xlen)={xlen_m*1e3:.1f}mm, z={zlim}mm, "
          f"C_total上限={cmax:.2e}", file=sys.stderr)

    if args.kind == "all":
        kinds = ["spread", "absorb", "efield"]
    elif args.kind == "both":
        kinds = ["spread", "absorb"]
    else:
        kinds = [args.kind]
    for kd in [k for k in kinds if k in ("spread", "absorb")]:
        norm = norm_s if kd == "spread" else norm_a
        label = L["d"] if kd == "spread" else L["ct"]
        print(f"[{kd}] 作成中 ...", file=sys.stderr)
        if args.mode in ("snap", "both"):
            fig, axes = build(kd, args.dual, xlen_m, norm, plt.cm.jet, label)
            picks = np.unique(np.linspace(0, len(frames) - 1,
                              min(args.nsnap, len(frames))).astype(int))
            for k in picks:
                render(fig, axes, kd, vns, frames[k], xlen_m, zlim,
                       args.rmax, norm, bin_m, levels)
                p = os.path.join(args.outdir,
                                 f"{kd}_snap_{k:03d}_istp{frames[k]['istp']}.png")
                fig.savefig(p, dpi=150)
                print(f"  保存: {p}", file=sys.stderr)
            plt.close(fig)
        if args.mode in ("anim", "both"):
            fig, axes = build(kd, args.dual, xlen_m, norm, plt.cm.jet, label)

            def upd(k, kd=kd, fig=fig, axes=axes, norm=norm):
                render(fig, axes, kd, vns, frames[k], xlen_m, zlim,
                       args.rmax, norm, bin_m, levels)
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

    if "efield" in kinds:
        print("[efield] 電場を作成中 ...", file=sys.stderr)
        res = read_flow_phi(args.flow)
        if res is None:
            print(f"  {args.flow} が読めないため電場はスキップ", file=sys.stderr)
        else:
            X, Y, PHI = res
            fig = plt.figure(figsize=(7.0, 4.6))
            ax = fig.add_axes([0.09, 0.14, 0.78, 0.74])
            cax = fig.add_axes([0.90, 0.14, 0.02, 0.74])
            cf = plot_efield(ax, X, Y, PHI, xlen_m)
            ax.set_xlim(0, xlen_m * 1e3)
            ax.set_ylim(*zlim)
            ax.set_xlabel(L["r"])
            ax.set_ylabel(L["z"])
            fig.suptitle("電場分布 |E| と電気力線" if JP
                         else "Electric field |E| and field lines", y=0.98)
            fig.colorbar(cf, cax=cax, label="|E| [V/m]")
            p = os.path.join(args.outdir, "efield_field.png")
            fig.savefig(p, dpi=150)
            plt.close(fig)
            print(f"  保存: {p}", file=sys.stderr)

    if "efieldt" in kinds:
        print("[efieldt] 時間変化する電場(電極+液滴空間電荷)を作成中 ...",
              file=sys.stderr)
        res = read_flow_phi(args.flow)
        if res is None:
            print(f"  {args.flow} が読めないのでスキップ", file=sys.stderr)
        else:
            Xf, Yf, PHI = res
            s = 10
            Xc, Zc, PHc = Xf[::s, ::s], Yf[::s, ::s], PHI[::s, ::s]
            xc1, zc1 = Xc[0, :], Zc[:, 0]
            dpz, dpx = np.gradient(PHc, zc1, xc1)
            Eex, Eez = -dpx, -dpz               # 電極(静的)場
            soft = (xc1[1] - xc1[0]) * 0.6
            # 色スケールは最終フレーム（液滴が最も多い）の空間電荷場で固定
            xq, yq, dq, cq = active(vns, frames[-1])
            if len(xq):
                e0x, e0z = efield_from_charges(xq, yq, droplet_charge(dq),
                                               Xc, Zc, xlen_m, soft)
                E0 = np.sqrt(e0x ** 2 + e0z ** 2)
                pos = E0[np.isfinite(E0) & (E0 > 0)]
                vmax = float(np.percentile(pos, 99)) if pos.size else 1.e4
            else:
                vmax = 1.e4
            vmin = vmax / 1e3
            levels = np.logspace(np.log10(vmin), np.log10(vmax), 12)
            fig = plt.figure(figsize=(7.0, 4.6))
            ax = fig.add_axes([0.09, 0.14, 0.78, 0.74])
            cax = fig.add_axes([0.90, 0.14, 0.02, 0.74])
            sm = ScalarMappable(norm=LogNorm(vmin, vmax), cmap="plasma")
            sm.set_array([])
            fig.colorbar(sm, cax=cax, label="|E| [V/m]")

            def updE(kk):
                ax.clear()
                xp, yp, dp, clt = active(vns, frames[kk])
                if len(xp):
                    Ex, Ez = efield_from_charges(xp, yp, droplet_charge(dp),
                                                 Xc, Zc, xlen_m, soft)
                else:
                    Ex = np.zeros_like(Xc)
                    Ez = np.zeros_like(Xc)
                Em = np.clip(np.sqrt(Ex ** 2 + Ez ** 2), vmin, vmax)
                ax.contourf(xc1 * 1e3, zc1 * 1e3, Em, levels=levels,
                            cmap="plasma", norm=LogNorm(vmin, vmax),
                            extend="both", zorder=1)
                try:
                    ax.streamplot(xc1 * 1e3, zc1 * 1e3, Ex, Ez, color="white",
                                  density=1.0, linewidth=0.5, arrowsize=0.7)
                except Exception:
                    pass
                if len(xp):
                    ax.scatter(xp * 1e3, yp * 1e3, s=2, c="cyan", alpha=0.4,
                               zorder=3)
                draw_geometry(ax, xlen_m)
                ax.set_xlim(0, xlen_m * 1e3)
                ax.set_ylim(*zlim)
                ax.set_xlabel(L["r"])
                ax.set_ylabel(L["z"])
                t = "液滴(空間電荷)がつくる電場 |E|" if JP \
                    else "|E| from droplet space charge"
                fig.suptitle(f'{t}   ({L["step"]}={frames[kk]["istp"]})',
                             y=0.98)
                return []

            anim = animation.FuncAnimation(fig, updE, frames=len(frames),
                                           blit=False)
            mp4 = os.path.join(args.outdir, "efield_animation.mp4")
            try:
                anim.save(mp4, writer=animation.FFMpegWriter(fps=args.fps,
                          bitrate=2400), dpi=130)
                print(f"  保存: {mp4}", file=sys.stderr)
            except Exception as e:
                gif = os.path.join(args.outdir, "efield_animation.gif")
                print(f"  MP4不可({e})->GIF", file=sys.stderr)
                anim.save(gif, writer=animation.PillowWriter(fps=args.fps),
                          dpi=100)
                print(f"  保存: {gif}", file=sys.stderr)
            plt.close(fig)

    print("完了。出力先: " + os.path.abspath(args.outdir), file=sys.stderr)


if __name__ == "__main__":
    main()
