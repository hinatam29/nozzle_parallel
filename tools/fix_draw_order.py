"""液滴の図で、点を描く順番をランダムにする（2026-11-18）。
これまでは spray.dat の順（出た順）に描いていたので、重なった所では一番新しい液滴が上に来て見えていた。
新しいのに遠くまで来ている液滴は速い液滴（10 nm 付近では小さい液滴）なので、帯電量を下げて
液滴が混ざらない計算では、実際より小さい液滴ばかりに見えていた。描く順番をランダムにして、
重なった所の見え方を液滴径の割合どおりにする（毎回同じ順になるように乱数の種を固定）。
使い方: python3 tools/fix_draw_order.py"""
import glob, os
old = "    return xp[m], yp[m], dp[m], clt[m]"
new = """    # 2026-11-18: 重なった所の見え方が偏らないように、描く順番をランダムにする（種は固定）
    o = np.random.default_rng(12345).permutation(int(m.sum()))
    return xp[m][o], yp[m][o], dp[m][o], clt[m][o]"""
n = 0
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
for f in glob.glob(os.path.join(root, '**', 'visualize.py'), recursive=True):
    s = open(f).read()
    if old in s and "permutation(int(m.sum()))" not in s:
        open(f, 'w').write(s.replace(old, new)); n += 1
print("直した visualize.py:", n)
