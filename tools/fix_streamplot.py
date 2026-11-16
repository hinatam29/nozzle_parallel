"""並列の visualize.py で電気力線が描かれていなかった不具合を直す（2026-11-16）。
原因: 格子の座標に丸め誤差があり、matplotlib の streamplot が「等間隔でない」とエラーを出し、
try/except で黙って飛ばされていた。等間隔の座標に置き直し、表示する高さの範囲だけで描く。
使い方: python3 tools/fix_streamplot.py   （kawatani の下の visualize.py をすべて直す）"""
import glob, os
old = '''    try:
        ax.streamplot(xc * 1e3, zc * 1e3, np.ma.masked_where(weak, Ex),
                      np.ma.masked_where(weak, Ez), color="black",
                      density=1.1, linewidth=0.5, arrowsize=0.8)
    except Exception:
        pass'''
new = '''    # 2026-11-16: 格子の座標の丸め誤差で streamplot が失敗していたので、等間隔の座標に置き直し、
    # 表示する高さ（SL_ZLIM）の範囲だけで描く（範囲全体で描くと表示部分の線が少なくなる）
    xs = np.linspace(xc[0], xc[-1], len(xc))
    zs = np.linspace(zc[0], zc[-1], len(zc))
    zl = globals().get("SL_ZLIM")
    sel = slice(None)
    if zl is not None:
        jj = np.where((zs * 1e3 >= zl[0] - 0.5) & (zs * 1e3 <= zl[1] + 0.5))[0]
        if len(jj) > 2:
            sel = slice(jj[0], jj[-1] + 1)
    try:
        ax.streamplot(xs * 1e3, zs[sel] * 1e3, np.ma.masked_where(weak, Ex)[sel],
                      np.ma.masked_where(weak, Ez)[sel], color="black",
                      density=1.1, linewidth=0.5, arrowsize=0.8, zorder=3)
    except Exception as ex:
        print("  電気力線を描けませんでした:", ex, file=sys.stderr)'''
old2 = "            cf = plot_efield(ax, X, Y, PHI, xlen_m)\n"
new2 = "            globals()[\"SL_ZLIM\"] = zlim\n" + old2
n = 0
for f in glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '**', 'visualize.py'), recursive=True):
    s = open(f).read()
    if old in s and s.count(old2) == 1:
        open(f, 'w').write(s.replace(old, new).replace(old2, new2)); n += 1
print('直した visualize.py:', n)
