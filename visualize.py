# -*- coding: utf-8 -*-
"""
静電噴霧シミュレーション結果の可視化プログラム（噴霧粒子のみ）
=================================================================

2 本ノズル並列配置の噴霧粒子を、2 通りの見方で可視化します。

  (1) 液滴の広がり      : 液滴の位置分布（両ノズルから出た噴霧の広がり）
  (2) 液滴の二酸化炭素吸収 : 液滴を「その液滴が吸収した二酸化炭素量 clt」で
                            色分け（対数スケール）。吸収の進み具合が分かる。

読み込むファイル
  spray.dat : 噴霧粒子データ（Tecplot point 形式）
              変数 = xp, yp, dp(粒径), clt(粒子ごとの吸収量), pout(状態)
  （気相の場データ flow.dat は使いません。）

生成物（既定は両方）
  spread_snapshot_*.png / spread_animation.(mp4|gif)   … 広がり
  absorb_snapshot_*.png / absorb_animation.(mp4|gif)   … 吸収

使い方（例）
  python visualize.py
  python visualize.py --kind spread          # 広がりだけ
  python visualize.py --kind absorb          # 吸収だけ
  python visualize.py --max-frames 18        # 正常範囲(istp<=17000)だけ
  python visualize.py --mode anim --fps 8

注意
  粒子は istp=18000 で数値発散(NaN)します。発散フレームは自動的に除外
  されますが、'--max-frames 18' を付けると正常範囲だけを確実に扱えます。

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
from matplotlib.colors import LogNorm


# ----------------------------------------------------------------------
# 日本語フォント
# ----------------------------------------------------------------------
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
    "r": "半径方向 r [m]" if JP else "radial  r [m]",
    "z": "軸方向 z [m]" if JP else "axial  z [m]",
    "cbar_clt": "液滴が吸収した二酸化炭素量 clt"
                if JP else "droplet absorbed carbon dioxide  clt",
    "step": "ステップ" if JP else "step",
    "t_spread": "液滴の広がり" if JP else "Droplet spread",
    "t_absorb": "液滴の二酸化炭素吸収"
                if JP else "Droplet carbon dioxide absorption",
    "n1": "ノズル1" if JP else "nozzle 1",
    "n2": "ノズル2" if JP else "nozzle 2",
}


# ----------------------------------------------------------------------
# Tecplot point 形式パーサ（3桁指数の E 省略・NaN に対応）
# ----------------------------------------------------------------------
_exp_fix = re.compile(r"(\d\.\d+)([+-]\d{2,3})(?![\d.])")


def fix_exp(s):
    return _exp_fix.sub(r"\1E\2", s)


def to_float(t):
    try:
        return float(t.replace("D", "E").replace("d", "e"))
    except ValueError:
        return np.nan


def read_tecplot(path, stride=1, max_frames=None):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    varnames = None
    frames = []
    zidx = -1
    with open(path, "r", errors="replace") as fh:
        pend, buf, keep = None, [], False

        def flush():
            nonlocal pend, buf
            if pend is None:
                return
            need = pend["npoints"] * (pend["ncols"] or 0)
            if keep and pend["ncols"] and len(buf) >= need:
                arr = np.array(buf[:need], dtype=np.float64).reshape(
                    pend["npoints"], pend["ncols"])
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
                istp = int(m.group(1)) if m else zidx
                mi = re.search(r"\bi\s*=\s*(\d+)", s, re.I)
                ni = int(mi.group(1)) if mi else 1
                keep = (zidx % stride == 0) and (
                    max_frames is None or len(frames) < max_frames)
                pend = dict(istp=istp, npoints=ni,
                            ncols=len(varnames) if varnames else None)
                buf = []
                continue
            if pend is None or not keep:
                continue
            toks = fix_exp(s).split()
            if pend["ncols"] is None:
                pend["ncols"] = len(toks)
            buf.extend(to_float(t) for t in toks)
        flush()
    if varnames is None:
        raise ValueError(f"{path}: variables 行が見つかりません")
    print(f"  -> {len(frames)} フレーム読み込み ({path})", file=sys.stderr)
    return varnames, frames


def active_droplets(varnames, frame):
    idx = {n: i for i, n in enumerate(varnames)}
    d = frame["data"]
    xp = d[:, idx.get("xp", 0)]
    yp = d[:, idx.get("yp", 1)]
    dp = d[:, idx["dp"]] if "dp" in idx else np.zeros(len(xp))
    clt = d[:, idx["clt"]] if "clt" in idx else np.zeros(len(xp))
    pout = d[:, idx["pout"]] if "pout" in idx else np.zeros(len(xp))
    m = np.isfinite(xp) & np.isfinite(yp) & (pout == 0)
    return xp[m], yp[m], dp[m], clt[m]


# ----------------------------------------------------------------------
# 全フレーム共通のスケール
# ----------------------------------------------------------------------
def compute_scales(vns, frames):
    xs, ys, dps, clts = [], [], [], []
    for fr in frames:
        xp, yp, dp, clt = active_droplets(vns, fr)
        if len(xp):
            xs.append(xp); ys.append(yp); dps.append(dp); clts.append(clt)
    if not xs:
        return dict(xlim=(0, 0.1), ylim=(0, 0.1), ytop=0.1, xlen=0.1,
                    dp=(1e-7, 1e-6), clt=(1e-24, 1e-21))
    x = np.concatenate(xs); y = np.concatenate(ys)
    dp = np.concatenate(dps); clt = np.concatenate(clts)
    xlen = float(x.max())
    ytop = float(y.max())
    ylo, yhi = float(y.min()), float(y.max())
    pad = max(0.06 * (yhi - ylo), 0.002)
    ylim = (ylo - pad, min(yhi + pad, ytop + pad))
    xpad = 0.02 * xlen
    xlim = (min(0.0, float(x.min())) - xpad, xlen + xpad)
    cpos = clt[clt > 0]
    if cpos.size:
        cmin = float(np.percentile(cpos, 5))
        cmax = float(np.percentile(cpos, 99))
        if cmin <= 0 or cmin >= cmax:
            cmin, cmax = cpos.min(), cpos.max()
    else:
        cmin, cmax = 1e-24, 1e-21
    return dict(xlim=xlim, ylim=ylim, ytop=ytop, xlen=xlen,
                dp=(float(dp.min()), float(dp.max())), clt=(cmin, cmax))


def sizes_from_dp(dp, dprange):
    lo, hi = dprange
    if hi <= lo:
        return np.full(np.shape(dp), 26.0)
    s = (np.asarray(dp) - lo) / (hi - lo)
    return 12.0 + 60.0 * np.clip(s, 0, 1)


def draw_nozzles(ax, sc):
    for xpos, name, ha in ((0.0, L["n1"], "left"),
                           (sc["xlen"], L["n2"], "right")):
        ax.scatter([xpos], [sc["ytop"]], marker="v", s=140, c="red",
                   edgecolors="k", linewidths=0.6, zorder=6, clip_on=False)
        ax.annotate(name, xy=(xpos, sc["ytop"]), xytext=(0, 9),
                    textcoords="offset points", ha=ha, va="bottom",
                    color="red", fontsize=9, zorder=6, clip_on=False)


# ----------------------------------------------------------------------
# 描画（kind = "spread" または "absorb"）
# ----------------------------------------------------------------------
def draw(ax, cax, vns, frame, kind, sc, norm):
    ax.clear()
    xp, yp, dp, clt = active_droplets(vns, frame)
    sizes = sizes_from_dp(dp, sc["dp"]) if len(xp) else 26.0
    sm = None
    if kind == "spread":
        if len(xp):
            ax.scatter(xp, yp, s=sizes, c="#1f77b4", alpha=0.6,
                       edgecolors="none", zorder=5)
        title = L["t_spread"]
    else:  # absorb
        if len(xp):
            sm = ax.scatter(xp, yp, s=sizes, c=np.clip(clt, norm.vmin, None),
                            cmap="inferno", norm=norm, alpha=0.9,
                            edgecolors="k", linewidths=0.2, zorder=5)
        title = L["t_absorb"]
    draw_nozzles(ax, sc)
    ax.set_xlim(*sc["xlim"]); ax.set_ylim(*sc["ylim"])
    ax.set_xlabel(L["r"]); ax.set_ylabel(L["z"])
    ax.set_title(f'{title}   ({L["step"]} = {frame.get("istp")})')
    if cax is not None:
        cax.cla()
        if sm is not None:
            plt.colorbar(sm, cax=cax, label=L["cbar_clt"])
        else:
            cax.axis("off")


def make_fig(kind):
    fig = plt.figure(figsize=(10.5, 4.6))
    if kind == "absorb":
        ax = fig.add_axes([0.08, 0.16, 0.80, 0.72])
        cax = fig.add_axes([0.90, 0.16, 0.02, 0.72])
    else:
        ax = fig.add_axes([0.08, 0.16, 0.86, 0.72])
        cax = None
    return fig, ax, cax


def run_kind(kind, vns, frames, sc, outdir, mode, nsnap, fps):
    norm = LogNorm(vmin=sc["clt"][0], vmax=sc["clt"][1]) \
        if kind == "absorb" else None
    if mode in ("snap", "both"):
        picks = np.unique(np.linspace(0, len(frames) - 1,
                                      min(nsnap, len(frames))).astype(int))
        for k in picks:
            fig, ax, cax = make_fig(kind)
            draw(ax, cax, vns, frames[k], kind, sc, norm)
            p = os.path.join(outdir,
                             f"{kind}_snapshot_{k:03d}_istp{frames[k]['istp']}.png")
            fig.savefig(p, dpi=140); plt.close(fig)
            print(f"  保存: {p}", file=sys.stderr)
    if mode in ("anim", "both"):
        fig, ax, cax = make_fig(kind)

        def update(k):
            draw(ax, cax, vns, frames[k], kind, sc, norm)
            return []

        anim = animation.FuncAnimation(fig, update, frames=len(frames),
                                       blit=False)
        mp4 = os.path.join(outdir, f"{kind}_animation.mp4")
        gif = os.path.join(outdir, f"{kind}_animation.gif")
        try:
            anim.save(mp4, writer=animation.FFMpegWriter(fps=fps, bitrate=2400),
                      dpi=130)
            print(f"  保存: {mp4}", file=sys.stderr)
        except Exception as e:
            print(f"  MP4 不可 ({e}) -> GIF", file=sys.stderr)
            try:
                anim.save(gif, writer=animation.PillowWriter(fps=fps), dpi=100)
                print(f"  保存: {gif}", file=sys.stderr)
            except Exception as e2:
                print(f"  GIF も失敗: {e2}", file=sys.stderr)
        plt.close(fig)


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="噴霧粒子の広がり・二酸化炭素吸収の可視化")
    ap.add_argument("--spray", default="spray.dat")
    ap.add_argument("--outdir", default="viz_out")
    ap.add_argument("--kind", default="both",
                    choices=["both", "spread", "absorb"])
    ap.add_argument("--mode", default="both", choices=["both", "anim", "snap"])
    ap.add_argument("--nsnap", type=int, default=6)
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--max-frames", type=int, default=None)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print("噴霧データを読み込み中 ...", file=sys.stderr)
    vns, frames = read_tecplot(args.spray, stride=args.stride,
                               max_frames=args.max_frames)
    if not frames:
        print("フレームがありません。", file=sys.stderr); return
    sc = compute_scales(vns, frames)
    print(f"  表示範囲 x={sc['xlim']} y={sc['ylim']}", file=sys.stderr)
    print(f"  clt色スケール(対数) {sc['clt'][0]:.2e}..{sc['clt'][1]:.2e}",
          file=sys.stderr)

    kinds = ["spread", "absorb"] if args.kind == "both" else [args.kind]
    for kd in kinds:
        print(f"[{kd}] を作成中 ...", file=sys.stderr)
        run_kind(kd, vns, frames, sc, args.outdir, args.mode,
                 args.nsnap, args.fps)
    print("完了。出力先: " + os.path.abspath(args.outdir), file=sys.stderr)


if __name__ == "__main__":
    main()
