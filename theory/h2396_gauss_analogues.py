# Gate A2: the repository estimator (theory/open_screen_deficit.one_slice) on the q = 13, s = 2, U = 1 analogues of the 599 slice.
# Output: logs/h2396/gauss_analogues.jsonl.  usage: python3 theory/h2396_gauss_analogues.py
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import open_screen_deficit as D
keys = [["b1", 79, [], [3, 2], 13, 2, []], ["b1", 105, [], [8], 13, 2, []], ["b1", 157, [], [3, 4], 13, 2, []],
        ["b1", 183, [], [7, 2], 13, 2, []], ["b1", 209, [], [16], 13, 2, []], ["b1", 235, [], [9, 2], 13, 2, []],
        ["b1", 261, [], [5, 4], 13, 2, []], ["b1", 287, [], [11, 2], 13, 2, []], ["b1", 339, [], [13, 2], 13, 2, []]]
for key in keys:
    r = dict(key=json.dumps(key), fam='b1', v=key[1], r=None, sGeff=None)
    t0 = time.time()
    rec = D.one_slice(r, 600, 6000, 3)
    rec['secs'] = round(time.time() - t0, 1)
    print(json.dumps(rec, default=float), flush=True)
