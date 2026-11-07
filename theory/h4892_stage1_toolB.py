"""H(4892) slice, step 1(a): the complete stage-1 solution list of method B of theory/open_screen_stage1_exact.py
(b_solve_np, --strong pruning), written in full (the tool itself stores only 20).  Method B returns each solution
once up to block order and block sign (blocks with y[0] < 0 discarded, quadruples sorted).
usage: python3 theory/h4892_stage1_toolB.py OUT.json [KEY]   (default KEY: the 4892 slice)"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import open_screen_stage1_exact as X

KEY = ["b1", 1223, [], [47, 2], 13, 4, [[0, 23]]]
if __name__ == "__main__":
    t0 = time.time()
    key = json.loads(sys.argv[2]) if len(sys.argv) > 2 else KEY
    st = X.b_setup(key)
    sols, nb, stats = X.b_solve_np(st, st['T1'])
    json.dump(dict(key=key, method="open_screen_stage1_exact.b_solve_np (strong)", blocks=nb, join=stats,
                   solutions=len(sols), secs=round(time.time() - t0, 1), list=[list(map(list, Y)) for Y in sols]),
              open(sys.argv[1], 'w'))
    print(len(sols), nb, stats, round(time.time() - t0, 1))
